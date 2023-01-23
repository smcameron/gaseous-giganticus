/*
	Copyright (C) 2014 Stephen M. Cameron
	Author: Stephen M. Cameron

	This file is part of Gaseous Giganticus.

	Gaseous Giganticus is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	Gaseous Giganticus is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Gaseous Giganticus; if not, write to the Free Software
	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <getopt.h>
#include <locale.h>
#include <sys/time.h>
#include <fenv.h>

#include "mtwist.h"
#include "mathutils.h"
#include "quat.h"
#include "open-simplex-noise.h"
#include "png_utils.h"
#include "pthread_util.h"

static struct osn_context *ctx;

#define NPARTICLES 8000000
static int particle_count = NPARTICLES;

static int nthreads = 4;
static int thread_iterations = 1;
static int user_threads = -1;
static const int image_threads = 6; /* for 6 faces of cubemap, don't change this */
static char *default_output_file_prefix = "gasgiant-";
static char *default_input_file = "gasgiant-input.png";
static char *vf_dump_file = NULL;
static char *flowmap_dump_file = NULL;
static int flowmap_impressions = 40;
static int flowmap_has_alpha = 0;
static int restore_vf_data = 0;
static char *output_file_prefix;
static char *input_file;
static char *cubemap_prefix = NULL;
static int nofade = 0;
static int stripe = 0;
static int sinusoidal = 1;
static int use_wstep = 0;
static int wstep_period = 10;
static float wstep = 0.0f;
#define FBM_DEFAULT_FALLOFF (0.5)
static float fbm_falloff = FBM_DEFAULT_FALLOFF; 
static float gain = 1.0;
static float ff[] = {
	1.0,
	FBM_DEFAULT_FALLOFF,
	FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF,
	FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF,
	FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF,
	FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF *
			FBM_DEFAULT_FALLOFF,
	FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF *
			FBM_DEFAULT_FALLOFF * FBM_DEFAULT_FALLOFF,
};
static int noise_levels = 4;
static const int max_noise_levels = 7;
static int cloudmode = 0;
static float fade_rate = -1.0;
static int save_texture_sequence = 0;
static int magic_fluid_flow = 0; /* 0 = skip fluid dynamics, 1 = do fluid dynamics (not yet implemented) */

#define DIM 1024 /* dimensions of cube map face images */
#define VFDIM 2048 /* dimension of velocity field. (2 * DIM) is reasonable */
static int vfdim = VFDIM;
#define FDIM ((float) (DIM))
#define XDIM DIM
#define YDIM DIM

static int niterations = 1000;
static const float default_noise_scale = 2.6;
static float noise_scale;
static float speed_multiplier = 1.0;
static float velocity_factor = 1200.0;
static float num_bands = 6.0f;
static float band_speed_factor = 2.9f;
static int band_speed_power = 1;
static float band_speed_powerf = 1.0;
static int vertical_bands = 1;
static float opacity = 1.0;
static float opacity_limit = 0.2;
static float pole_attenuation = 0.5;

static char *start_image;
static unsigned char *cubemap_image[6] = { 0 };
static int start_image_width, start_image_height, start_image_has_alpha, start_image_bytes_per_row;
static unsigned char *output_image[6];
static unsigned char *flowmap_image[6];
static int image_save_period = 25;
static float w_offset = 0.0;
static int random_mode;
static int large_pixels = 0;
static int export_equirect_image = 0;
static int equirect_height = 0;
static float cache_aware = (6.0f * XDIM * XDIM) / (float) NPARTICLES;

struct timing_data {
	struct timeval begin;
	struct timeval end;
	double elapsed;
};

/* velocity field for 6 faces of a cubemap */
struct velocity_field {
	union vec3 v[6][VFDIM][VFDIM];
};

static struct velocity_field *vf, *old_vf; /* old_vf is used in advection and divergence removal steps */

struct color {
	float r, g, b, a;
};

struct color darkest_color;
const struct color hot_pink = { 1.0f, 0.07843f, 0.57647f, 0.01 };
static int use_hot_pink = 0;

/* face, i, j -- coords on a cube map */
struct fij {
	int f, i, j;
};

/* particles have color, and exist on the surface of a sphere. */
static struct particle {
	union vec3 pos;
	struct color c;
	struct fij fij;
} *particle;

static union vec3 *prev_particle_pos; /* previous particle positions used in (optionally) calculating flow map */

struct movement_thread_info {
	int first_particle, last_particle;
	pthread_t thread;
	struct particle *p; /* ptr to particle[], above */
	struct velocity_field *vf; /* ptr to vf, above. */
	int iterations;
};

union cast {
	double d;
	long l;
};

#define MAXVORTICES 200
static struct vortex {
	union vec3 p;
	union quat q, reverseq;
	float r;
	float angular_vel;
} vort[MAXVORTICES];
static int nvortices = 0;
static float vortex_size = 0.04;
static float vortex_size_variance =  0.02;
static float vortex_band_threshold = 0.2;

static void backspace(int columns)
{
	int i;

	for (i = 0; i < columns; i++)
		printf("\033[D");
	fflush(stdout);
	for (i = 0; i < columns; i++)
		printf(" ");
	fflush(stdout);
	for (i = 0; i < columns; i++)
		printf("\033[D");
	fflush(stdout);
}

static double timeval_difference(struct timeval t1, struct timeval t2)
{
	double elapsed_time;
	/* compute and return the elapsed time in millisec */
	elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000.0;      /* sec to ms */
	elapsed_time += (t2.tv_usec - t1.tv_usec) / 1000.0;   /* us to ms */
	return elapsed_time;
}

static inline int float2int(double d)
{
	volatile union cast c;
	c.d = d + 6755399441055744.0;
	return c.l;
}

static float alphablendcolor(float underchannel, float underalpha, float overchannel, float overalpha)
{
	return overchannel * overalpha + underchannel * underalpha * (1.0 - overalpha);
}

static struct color combine_color(struct color *oc, struct color *c)
{
	struct color nc;

	nc.a = (c->a + oc->a * (1.0f - c->a));
	nc.r = alphablendcolor(oc->r, oc->a, c->r, c->a) / nc.a;
	nc.b = alphablendcolor(oc->b, oc->a, c->b, c->a) / nc.a;
	nc.g = alphablendcolor(oc->g, oc->a, c->g, c->a) / nc.a;
	return nc;
}

static void fade_out_background(int f, struct color *c)
{
	int p, i, j;
	struct color oc, nc;
	unsigned char *pixel;

	for (j = 0; j < DIM; j++) {
		for (i = 0; i < DIM; i++) {
			p = j * DIM + i;
			pixel = &output_image[f][p * 4];
			if (!cloudmode) {
				oc.r = (float) pixel[0] / 255.0f;
				oc.g = (float) pixel[1] / 255.0f;
				oc.b = (float) pixel[2] / 255.0f;
				oc.a = 1.0;
				nc = combine_color(&oc, c);
				pixel[0] = float2int(255.0f * nc.r) & 0xff;
				pixel[1] = float2int(255.0f * nc.g) & 0xff;
				pixel[2] = float2int(255.0f * nc.b) & 0xff;
			} else {
				pixel[3] = float2int((float) pixel[3] * 0.9);
			}
		}
	}
}

static union vec3 fij_to_xyz(int f, int i, int j, const int dim);

static float fbmnoise4_1oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w);
}

static float fbmnoise4_2oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w);
}

static float fbmnoise4_3oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w) +
		ff[2] * open_simplex_noise4(ctx, 4.0f * x, 4.0f * y, 4.0f * z, 4.0f * w);
}

static float fbmnoise4_4oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w) +
		ff[2] * open_simplex_noise4(ctx, 4.0f * x, 4.0f * y, 4.0f * z, 4.0f * w) +
		ff[3] * open_simplex_noise4(ctx, 8.0f * x, 8.0f * y, 8.0f * z, 8.0f * w);
}

static float fbmnoise4_5oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w) +
		ff[2] * open_simplex_noise4(ctx, 4.0f * x, 4.0f * y, 4.0f * z, 4.0f * w) +
		ff[3] * open_simplex_noise4(ctx, 8.0f * x, 8.0f * y, 8.0f * z, 8.0f * w) +
		ff[4] * open_simplex_noise4(ctx, 16.0f * x, 16.0f * y, 16.0f * z, 16.0f * w);
}

static float fbmnoise4_6oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w) +
		ff[2] * open_simplex_noise4(ctx, 4.0f * x, 4.0f * y, 4.0f * z, 4.0f * w) +
		ff[3] * open_simplex_noise4(ctx, 8.0f * x, 8.0f * y, 8.0f * z, 8.0f * w) +
		ff[4] * open_simplex_noise4(ctx, 16.0f * x, 16.0f * y, 16.0f * z, 16.0f * w) +
		ff[5] * open_simplex_noise4(ctx, 32.0f * x, 32.0f * y, 32.0f * z, 32.0f * w);
}

static float fbmnoise4_7oct(float x, float y, float z, float w)
{
	return	ff[0] * open_simplex_noise4(ctx, x, y, z, w) +
		ff[1] * open_simplex_noise4(ctx, 2.0f * x, 2.0f * y, 2.0f * z, 2.0f * w) +
		ff[2] * open_simplex_noise4(ctx, 4.0f * x, 4.0f * y, 4.0f * z, 4.0f * w) +
		ff[3] * open_simplex_noise4(ctx, 8.0f * x, 8.0f * y, 8.0f * z, 8.0f * w) +
		ff[4] * open_simplex_noise4(ctx, 16.0f * x, 16.0f * y, 16.0f * z, 16.0f * w) +
		ff[5] * open_simplex_noise4(ctx, 32.0f * x, 32.0f * y, 32.0f * z, 32.0f * w) +
		ff[6] * open_simplex_noise4(ctx, 64.0f * x, 64.0f * y, 64.0f * z, 64.0f * w);
}

static float (*const fbmnoise4[])(float x, float y, float z, float w) = {
	fbmnoise4_1oct,
	fbmnoise4_2oct,
	fbmnoise4_3oct,
	fbmnoise4_4oct,
	fbmnoise4_5oct,
	fbmnoise4_6oct,
	fbmnoise4_7oct,
};

static void paint_particle(int face, int i, int j, struct color *c, const int octaves)
{
	unsigned char *pixel;
	int p;
	struct color oc, nc;
	float (*const noise)(float x, float y, float z, float w) = fbmnoise4[octaves];

	if (i < 0 || i > DIM - 1 || j < 0 || j > DIM - 1) {
		/* FIXME: We get a handful of these, don't know why.
		 * Some are 1024, some seem to be -MAXINT
		 */
		/* printf("i, j = %d, %d!!!!!1\n", i, j);  */
		return;
	}
	p = j * DIM + i;
	pixel = &output_image[face][p * 4];
#if 0
	pixel[0] = (unsigned char) (255.0f * c->r);
	pixel[1] = (unsigned char) (255.0f * c->g);
	pixel[2] = (unsigned char) (255.0f * c->b);

	pixel[3] = 255;
	return;
#else
	/* FIXME, this is inefficient */
	oc.r = (float) pixel[0] / 255.0f;
	oc.g = (float) pixel[1] / 255.0f;
	oc.b = (float) pixel[2] / 255.0f;
	oc.a = 1.0;
	nc = combine_color(&oc, c);
	if (!cloudmode) {
		pixel[0] = float2int(255.0f * nc.r) & 0xff;
		pixel[1] = float2int(255.0f * nc.g) & 0xff;
		pixel[2] = float2int(255.0f * nc.b) & 0xff;
		pixel[3] = 255;
	} else {
		union vec3 v;
		float n, m;
		v = fij_to_xyz(face, i, j, DIM);
		vec3_normalize_self(&v);
		vec3_mul_self(&v, 3.6 * noise_scale);
		n = noise(v.v.x, v.v.y, v.v.z, (w_offset + 10.0) * 3.33f);
		if (n > 0.5f)
			n = n * (1.0 + n - 0.5);
		if (n < 0.0f)
			n = n * (1.0 + n + 0.25);
		if (n > 0.666f)
			n = 0.666f;
		else if (n < -0.333f)
			n = -0.333f;

		// else if (n < -1.0f)
			//n = -1.0f;
		// n = 0.5f * (n + 1.0f);
		n = n + 0.334f;
		m = (nc.r + nc.g + nc.b) / 3.0f;
#if 0
		pixel[0] = float2int(255.0f * m) & 0xff;
		//pixel[0] = 255;
		pixel[1] = pixel[0];
		pixel[2] = pixel[0];
		pixel[3] = float2int(n * (float) pixel[0]) & 0xff;
#endif
		pixel[0] = 255;
		pixel[1] = 255;
		pixel[2] = 255;
		pixel[3] = float2int(255.0f * m * n) & 0xff;
	}
#endif
}

/* convert from cubemap coords to cartesian coords on surface of sphere */
static union vec3 fij_to_xyz(int f, int i, int j, const int dim)
{
	union vec3 answer;

	switch (f) {
	case 0:
		answer.v.x = (float) (i - dim / 2) / (float) dim;
		answer.v.y = -(float) (j - dim / 2) / (float) dim;
		answer.v.z = 0.5;
		break;
	case 1:
		answer.v.x = 0.5;
		answer.v.y = -(float) (j - dim / 2) / (float) dim;
		answer.v.z = -(float) (i - dim / 2) / (float) dim;
		break;
	case 2:
		answer.v.x = -(float) (i - dim / 2) / (float) dim;
		answer.v.y = -(float) (j - dim / 2) / (float) dim;
		answer.v.z = -0.5;
		break;
	case 3:
		answer.v.x = -0.5;
		answer.v.y = -(float) (j - dim / 2) / (float) dim;
		answer.v.z = (float) (i - dim / 2) / (float) dim;
		break;
	case 4:
		answer.v.x = (float) (i - dim / 2) / (float) dim;
		answer.v.y = 0.5;
		answer.v.z = (float) (j - dim / 2) / (float) dim;
		break;
	case 5:
		answer.v.x = (float) (i - dim / 2) / (float) dim;
		answer.v.y = -0.5;
		answer.v.z = -(float) (j - dim / 2) / (float) dim;
		break;
	}
	vec3_normalize_self(&answer);
	return answer;
}

/* convert from cartesian coords on surface of a sphere to floating point cubemap coords (face encoded in z) */
static union vec3 xyz_to_fij_float(const union vec3 *p, const int dim, int preferred_face)
{
	union vec3 answer;
	union vec3 t;
	float f, i, j;
	float d;
	const float fdim = (float) dim;

	vec3_normalize(&t, p);

	if (fabs(t.v.x) > fabs(t.v.y)) {
		if (fabs(t.v.x) > fabs(t.v.z)) {
			/* x is longest leg */
			d = fabs(t.v.x);
			if (t.v.x < 0) {
				f = 3;
				i = (t.v.z / d) * fdim * 0.5 + 0.5 * (float) fdim;
			} else {
				f = 1;
				i = (-t.v.z / d)  * fdim * 0.5 + 0.5 * fdim;
			}
		} else {
			/* z is longest leg */
			d = fabs(t.v.z);
			if (t.v.z < 0) {
				f = 2;
				i = (-t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			} else {
				f = 0;
				i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			}
		}
		j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
	} else {
		/* x is not longest leg, y or z must be. */
		if (fabs(t.v.y) > fabs(t.v.z)) {
			/* y is longest leg */
			d = fabs(t.v.y);
			if (t.v.y < 0) {
				f = 5;
				j = (-t.v.z / d) * fdim * 0.5 + 0.5 * fdim;
			} else {
				f = 4;
				j = (t.v.z / d) * fdim * 0.5 + 0.5 * fdim;
			}
			i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
		} else {
			/* z is longest leg */
			d = fabs(t.v.z);
			if (t.v.z < 0) {
				f = 2;
				i = (-t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			} else {
				f = 0;
				i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			}
			j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
		}
	}

	if (preferred_face != -1 && f != preferred_face) {
		f = preferred_face;
		switch (preferred_face) {
		case 0:
			d = fabs(t.v.z);
			i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		case 1:
			d = fabs(t.v.x);
			i = (-t.v.z / d)  * fdim * 0.5 + 0.5 * fdim;
			j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		case 2:
			d = fabs(t.v.z);
			i = (-t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		case 3:
			d = fabs(t.v.x);
			i = (t.v.z / d)  * fdim * 0.5 + 0.5 * fdim;
			j = (-t.v.y / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		case 4:
			d = fabs(t.v.y);
			i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			j = (t.v.z / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		case 5:
			d = fabs(t.v.y);
			i = (t.v.x / d) * fdim * 0.5 + 0.5 * fdim;
			j = (-t.v.z / d) * fdim * 0.5 + 0.5 * fdim;
			break;
		}
	}

	answer.v.x = i;
	answer.v.y = j;
	answer.v.z = f;

	return answer;
}

/* convert from cartesian coords on surface of a sphere to cubemap coords */
static struct fij xyz_to_fij(const union vec3 *p, const int dim)
{
	struct fij answer;
	union vec3 t;
	int f, i, j;
	float d;
	const float fdim = (float) dim;

	vec3_normalize(&t, p);

	if (fabs(t.v.x) > fabs(t.v.y)) {
		if (fabs(t.v.x) > fabs(t.v.z)) {
			/* x is longest leg */
			d = fabs(t.v.x);
			if (t.v.x < 0) {
				f = 3;
				i = float2int((t.v.z / d) * fdim * 0.5 + 0.5 * (float) fdim);
			} else {
				f = 1;
				i = float2int((-t.v.z / d)  * fdim * 0.5 + 0.5 * fdim);
			}
		} else {
			/* z is longest leg */
			d = fabs(t.v.z);
			if (t.v.z < 0) {
				f = 2;
				i = float2int((-t.v.x / d) * fdim * 0.5 + 0.5 * fdim);
			} else {
				f = 0;
				i = float2int((t.v.x / d) * fdim * 0.5 + 0.5 * fdim);
#if 0
				/* FIXME: we get this sometimes, not sure why. */
				if (i < 0 || i > 1023)
					printf("i = %d!!!!!!!!!!!!!!!!\n", i);
#endif
			}
		}
		j = float2int((-t.v.y / d) * fdim * 0.5 + 0.5 * fdim);
	} else {
		/* x is not longest leg, y or z must be. */
		if (fabs(t.v.y) > fabs(t.v.z)) {
			/* y is longest leg */
			d = fabs(t.v.y);
			if (t.v.y < 0) {
				f = 5;
				j = float2int((-t.v.z / d) * fdim * 0.5 + 0.5 * fdim);
			} else {
				f = 4;
				j = float2int((t.v.z / d) * fdim * 0.5 + 0.5 * fdim);
			}
			i = float2int((t.v.x / d) * fdim * 0.5 + 0.5 * fdim);
		} else {
			/* z is longest leg */
			d = fabs(t.v.z);
			if (t.v.z < 0) {
				f = 2;
				i = float2int((-t.v.x / d) * fdim * 0.5 + 0.5 * fdim);
			} else {
				f = 0;
				i = float2int((t.v.x / d) * fdim * 0.5 + 0.5 * fdim);
			}
			j = float2int((-t.v.y / d) * fdim * 0.5 + 0.5 * fdim);
		}
	}

	answer.f = f;
	answer.i = i;
	answer.j = j;

	/* FIXME: some other problem makes these checks necessary for now. */
	if (answer.f < 0)
		answer.f = 0;
	else if (answer.f > 5)
		answer.f = 5;
	if (answer.i < 0)
		answer.i = 0;
	else if (answer.i >= dim)
		answer.i = dim - 1;
	if (answer.j < 0)
		answer.j = 0;
	else if (answer.j >= dim)
		answer.j = dim - 1;

	return answer;
}

static const float face_to_xdim_multiplier[] = { 0.25, 0.5, 0.75, 0.0, 0.25, 0.25 };
static const float face_to_ydim_multiplier[] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0,
						0.0, 2.0 / 3.0 };

static void add_jitter(union vec3 *v)
{
	union vec3 jitter;

	jitter.v.x = (((rand() % 10000) / 10000.0) - 0.5) / (2.5f * XDIM);
	jitter.v.y = (((rand() % 10000) / 10000.0) - 0.5) / (2.5f * XDIM);
	jitter.v.z = (((rand() % 10000) / 10000.0) - 0.5) / (2.5f * XDIM);
	vec3_add_self(v, &jitter);
	vec3_normalize_self(v);
}

/* Used if nparticles is >= number of pixels in cubemap */
static void dense_cache_friendly_location(int p, int nparticles, union vec3 *pos)
{
	int f, i, j, n;

	/* choose a face of the cubemap based on particle number */
	f  = p / (nparticles / 6);
	if (f > 5) /* integer math, e.g. 7999997 / (8000000 / 6) == 5, but 7999998 / (8000000 / 6) == 6 */
		f = 5;

	int locations_per_face = XDIM * XDIM;

	n = p % locations_per_face;
	i = n / XDIM;
	j = n % XDIM;
	*pos = fij_to_xyz(f, i, j, XDIM);
	add_jitter(pos); /* gets rid of some weird artifacts */
	vec3_mul_self(pos, (float) XDIM / 2.0f);
}

/* Used if nparticles is < number of pixels in cubemap */
static void sparse_cache_friendly_location(int p, int nparticles, union vec3 *pos)
{
	int f, i, j, n;

	/* choose a face of the cubemap based on particle number */
	f  = p / (nparticles / 6);
	if (f > 5) /* integer math, e.g. 7999997 / (8000000 / 6) == 5, but 7999998 / (8000000 / 6) == 6 */
		f = 5;

	int locations_per_face = XDIM * XDIM;
	float pixels_per_particle = (float) locations_per_face / (float) (nparticles / 6.0f);
	n = (int) ((p % (nparticles / 6)) * pixels_per_particle);
	n = n + (rand() % (int) pixels_per_particle);
	i = n / XDIM;
	j = n % XDIM;
	if (i > XDIM)
		i = i % XDIM;
	*pos = fij_to_xyz(f, i, j, XDIM);
	add_jitter(pos); /* gets rid of some weird artifacts */
	vec3_mul_self(pos, (float) XDIM / 2.0f);
}

/* place particles randomly on the surface of a sphere */
static void init_particles(struct particle **pp, union vec3 **prev_particle_pos, const int nparticles)
{
	float xo, yo;
	const int bytes_per_pixel = start_image_has_alpha ? 4 : 3;
	unsigned char *pixel;
	int pn, px, py;
	struct fij fij = { 0, 0, 0 };
	struct particle *p;
	int cache_aware_limit = (int) (cache_aware * nparticles);

	printf("Initializing %d particles", nparticles); fflush(stdout);
	*pp = malloc(sizeof(**pp) * nparticles);
	p = *pp;

	if (flowmap_dump_file)
		*prev_particle_pos = calloc(nparticles, sizeof(**prev_particle_pos));

	for (int i = 0; i < nparticles; i++) {
		if (i < cache_aware_limit) {
			if (cache_aware_limit >= XDIM * XDIM * 6)
				dense_cache_friendly_location(i, cache_aware_limit, &p[i].pos);
			else
				sparse_cache_friendly_location(i, cache_aware_limit, &p[i].pos);
		} else {
			random_point_on_sphere((float) XDIM / 2.0f, &p[i].pos.v.x, &p[i].pos.v.y, &p[i].pos.v.z);
		}
		if (cubemap_prefix) {
			fij = xyz_to_fij(&p[i].pos, DIM);
			if (fij.i < 0 || fij.i > DIM || fij.j < 0 || fij.j > DIM) {
				printf("BAD fij: %d,%d\n", fij.i, fij.j);
			}
			xo = fij.i * start_image_width / (float) XDIM;
			yo = fij.j * start_image_height / (float) YDIM;
			px = ((int) xo) * bytes_per_pixel;
			py = ((int) yo) * start_image_bytes_per_row;
		} else if (!stripe && !sinusoidal) {
			fij = xyz_to_fij(&p[i].pos, DIM);
			if (fij.i < 0 || fij.i > DIM || fij.j < 0 || fij.j > DIM) {
				printf("BAD fij: %d,%d\n", fij.i, fij.j);
			}
			/* The divide by 4.0 or 3.0 below is there because the faces are laid out like:
			 *   #
			 *  ####   (3 high, and 4 wide.)
			 *   #
			 */
			xo = XDIM * face_to_xdim_multiplier[fij.f] + fij.i / 4.0;
			yo = YDIM * face_to_ydim_multiplier[fij.f] + fij.j / 3.0;
			xo = xo / (float) XDIM;
			yo = yo / (float) YDIM;
			px = ((int) (xo * start_image_width)) * bytes_per_pixel;
			py = (((int) (yo * start_image_height)) * start_image_bytes_per_row);
			//printf("x,y = %d, %d, pn = %d\n", px, py, pn);
		} else if (stripe) {
			if (vertical_bands)
				yo = (p[i].pos.v.z + (float) DIM / 2.0) / (float) DIM;
			else
				yo = (p[i].pos.v.y + (float) DIM / 2.0) / (float) DIM;
			xo = 0.5; 
			px = ((int) (xo * start_image_width)) * bytes_per_pixel;
			py = (((int) (yo * start_image_height)) * start_image_bytes_per_row);
			//printf("xo = %f, yo=%f, px = %d, py = %d\n", xo, yo, px, py);
		} else { /* sinusoidal */
			if (!vertical_bands) {
				float abs_cos_lat, longitude, latitude, x1, x2;

				latitude = asinf(p[i].pos.v.y / ((float) DIM / 2.0f));
				yo = ((latitude + 0.5 * M_PI) / M_PI);
				abs_cos_lat = fabs(cosf(latitude));
				x2 = abs_cos_lat * 0.5 + 0.5;
				x1 = -abs_cos_lat * 0.5 + 0.5;
				longitude = atan2f(p[i].pos.v.z, p[i].pos.v.x);
				xo = (cos(longitude) * abs_cos_lat * 0.5f) * (x2 - x1) + 0.5f;
				px = ((int) (xo * start_image_width)) * bytes_per_pixel;
				py = (((int) (yo * start_image_height)) * start_image_bytes_per_row);
			} else {
				float abs_cos_lat, longitude, latitude, x1, x2;

				latitude = asinf(p[i].pos.v.z / ((float) DIM / 2.0f));
				yo = ((latitude + 0.5 * M_PI) / M_PI);
				abs_cos_lat = fabs(cosf(latitude));
				x2 = abs_cos_lat * 0.5 + 0.5;
				x1 = -abs_cos_lat * 0.5 + 0.5;
				longitude = atan2f(p[i].pos.v.y, p[i].pos.v.x);
				xo = (cos(longitude) * abs_cos_lat * 0.5f) * (x2 - x1) + 0.5f;
				px = ((int) (xo * start_image_width)) * bytes_per_pixel;
				py = (((int) (yo * start_image_height)) * start_image_bytes_per_row);
			}
		}
		pn = py + px;
		//printf("pn = %d\n", pn);
		if (!cubemap_prefix)
			pixel = (unsigned char *) &start_image[pn];
		else
			pixel = (unsigned char *) &cubemap_image[fij.f][pn];
		p[i].c.r = ((float) pixel[0]) / 255.0f;
		p[i].c.g = ((float) pixel[1]) / 255.0f;
		p[i].c.b = ((float) pixel[2]) / 255.0f;
		p[i].c.a = start_image_has_alpha ? (float) pixel[3] / 255.0 : 1.0;
		p[i].c.a = 1.0; //start_image_has_alpha ? (float) pixel[3] / 255.0 : 1.0;
	}
	printf("\n");
}

/* compute the noise gradient at the given point on the surface of a sphere */
static union vec3 noise_gradient(union vec3 position, float w, float noise_scale, const int octaves)
{
	union vec3 g;
	const float dx = noise_scale * (0.05f / (float) DIM);
	const float dy = noise_scale * (0.05f / (float) DIM);
	const float dz = noise_scale * (0.05f / (float) DIM);
	float (*const noise)(float x, float y, float z, float w) = fbmnoise4[octaves];

	g.v.x = noise(position.v.x + dx, position.v.y, position.v.z, w) -
		noise(position.v.x - dx, position.v.y, position.v.z, w);
	g.v.y = noise(position.v.x, position.v.y + dy, position.v.z, w) -
		noise(position.v.x, position.v.y - dy, position.v.z, w);
	g.v.z = noise(position.v.x, position.v.y, position.v.z + dz, w) -
		noise(position.v.x, position.v.y, position.v.z - dz, w);
	return g;
}

static union vec3 curl2(union vec3 pos, union vec3 normalized_pos,
			const float pos_magnitude, union vec3 noise_gradient)
{
	union vec3 pos_plus_noise, proj_ng, rotated_ng;
	union quat rotation;

	/* project noise gradient onto sphere surface */
	vec3_add(&pos_plus_noise, &pos, &noise_gradient);
	vec3_normalize_self(&pos_plus_noise);
	vec3_mul_self(&pos_plus_noise, pos_magnitude);
	vec3_sub(&proj_ng, &pos_plus_noise, &pos);

	/* rotate projected noise gradient 90 degrees about pos. */
	quat_init_axis_v(&rotation, &normalized_pos, M_PI / 2.0);
	quat_rot_vec(&rotated_ng, &proj_ng, &rotation);
	return rotated_ng;
}

static float calculate_band_speed(float latitude)
{
	/* to the 5th power to make regions between bands wider. */
	return ((1 - pole_attenuation) + pole_attenuation *
		cosf(latitude)) * powf(cosf(latitude * num_bands), band_speed_powerf) * band_speed_factor;
}

struct velocity_field_thread_info {
	pthread_t thread;
	int f; /* face */
	float w;
	struct velocity_field *vf;
	double elapsed_time_ms;
};

/* Update 1 face of the 6 velocity maps (1 face for each side of the cubemap) */ 
static void *update_velocity_field_thread_fn(void *info)
{
	struct velocity_field_thread_info *t = info;
	int f = t->f;
	float w = t->w;
	struct velocity_field *vf = t->vf;
	const int octaves = noise_levels - 1;
	struct timeval start_time, end_time;

	int i, j, vortex;
	union vec3 v, c, ng;

	gettimeofday(&start_time, NULL);
	for (i = 0; i < vfdim; i++) {
		for (j = 0; j < vfdim; j++) {
			float band_speed, angle;
			union vec3 ov, bv;

			v = fij_to_xyz(f, i, j, vfdim);
			ov = v;
			vec3_mul_self(&v, noise_scale);
			ng = noise_gradient(v, w * noise_scale, noise_scale, octaves);
			c = curl2(v, ov, noise_scale, ng);
			vec3_mul(&vf->v[f][i][j], &c, velocity_factor);

			/* calculate counter rotating band influence */
			if (num_bands != 0) {
				if (vertical_bands) {
					angle = asinf(ov.v.z);
					band_speed = calculate_band_speed(angle);
					bv.v.x = -ov.v.y;
					bv.v.y = ov.v.x;
					bv.v.z = 0;
				} else {
					angle = asinf(ov.v.y);
					band_speed = calculate_band_speed(angle);
					bv.v.x = ov.v.z;
					bv.v.z = -ov.v.x;
					bv.v.y = 0;
				}
				if (fabsf(vec3_magnitude(&bv)) > 1e-20) { /* watch out for NaNs */
					vec3_normalize_self(&bv);
					vec3_mul_self(&bv, band_speed);
					vec3_add_self(&vf->v[f][i][j], &bv);
				}
			}

			/* Calculate artificial vortex influence. */
			for (vortex = 0; vortex < nvortices; vortex++) {
				union vec3 diff;
				union vec3 p2d, vortex_vel, rotated_vortex_vel;
				float dist;
				struct vortex *vor = &vort[vortex];
				vec3_sub(&diff, &vor->p, &ov);
				dist = vec3_magnitude(&diff);
				if (dist > vor->r)
					continue; /* This vortex's influence is zero. */
				/* rotate this postition to front and center */
				quat_rot_vec(&p2d, &ov, &vor->q);

				/* do a 2d rotation about z axis by at most 5 degrees */
				float angle = (vor->angular_vel  * M_PI / 180.0) * sin(M_PI * dist / vor->r);
				float nx = p2d.v.x * cos(angle) - p2d.v.y * sin(angle);
				float ny = p2d.v.x * sin(angle) + p2d.v.y * cos(angle);
				vortex_vel.v.x = XDIM * (nx - p2d.v.x);
				vortex_vel.v.y = XDIM * (ny - p2d.v.y);
				vortex_vel.v.z = 0.0f;

				/* Now rotate the velocity back from front and center and add contribution */
				quat_rot_vec(&rotated_vortex_vel, &vortex_vel, &vor->reverseq);
				vec3_add_self(&vf->v[f][i][j], &rotated_vortex_vel);
			}
		}
	}
	gettimeofday(&end_time, NULL);
	t->elapsed_time_ms = timeval_difference(start_time, end_time);
	return NULL;
}

/* compute velocity field for all cells in cubemap.  It is scaled curl of gradient of noise field */
static void update_velocity_field(struct velocity_field *vf, float w, int *use_wstep)
{
	struct velocity_field_thread_info t[6];
	void *status;
	int f, rc;
	struct timeval vfbegin, vfend;
	static int first_time = 1;
	char *msg = "Calculating velocity field";

	gettimeofday(&vfbegin, NULL);
	if (first_time)
		printf("Velocity field is 6x %dx%d\n", VFDIM, VFDIM);
	printf("%s", msg); fflush(stdout);
	for (f = 0; f < 6; f++) {
		t[f].f = f;
		t[f].w = w;
		t[f].vf = vf;
		rc = create_thread(&t[f].thread, update_velocity_field_thread_fn, &t[f], "veloc-field", 0);
		if (rc)
			fprintf(stderr, "%s: create_thread failed: %s\n",
					__func__, strerror(errno));
	}
	for (f = 0; f < 6; f++) {
		int rc = pthread_join(t[f].thread, &status);
		if (rc)
			fprintf(stderr, "%s: pthread_join failed: %s\n",
					__func__, strerror(errno));
	}
	gettimeofday(&vfend, NULL);

	if (first_time) {
		printf("\n");
		for (f = 0; f < 6; f++)
			printf("Thread %d required %g ms\n", f, t[f].elapsed_time_ms);
		printf("Velocity field computed in %lu seconds, running simulation\n",
			vfend.tv_sec - vfbegin.tv_sec);
		first_time = 0;
	} else {
		backspace(strlen(msg));
	}
	if (*use_wstep)
		(*use_wstep)++;
}

static void check_vf_dump_file(char *filename)
{
	struct stat statbuf;

	if (!filename)
		return;

	int rc = stat(filename, &statbuf);
	if (rc == 0 && !restore_vf_data) /* file exists... */
		printf("File %s already exists, velocity field will not be dumped.\n",
				filename);
}

static void dump_velocity_field(char *filename_template, struct velocity_field *vf, int use_wstep)
{
	int fd;
	int bytes_left = (int) sizeof(*vf);
	int bytes_written = 0;
	char *filename = NULL;

	if (!filename_template || restore_vf_data)
		return;

	printf("\n"); fflush(stdout);

	if (!use_wstep) {
		filename = strdup(filename_template);
	} else {
		filename = malloc(strlen(filename_template) + 100);
		sprintf(filename, "%s-%d", filename_template, use_wstep);
	}
	fd = open(filename, O_RDWR | O_CREAT | O_EXCL, 0644);
	if (fd < 0) {
		fprintf(stderr, "Cannot create '%s': %s. Velocity field not dumped.\n",
			filename, strerror(errno));
		free(filename);
		return;
	}

	/* FIXME: this is quick and dirty, not endian clean, etc. */
	do {
		unsigned char *x = (unsigned char *) vf;
		int rc = write(fd, &x[bytes_written], bytes_left);
		if (rc < 0 && errno == EINTR)
			continue;
		if (rc < 0) {
			fprintf(stderr, "Error writing to '%s': %s\n", filename, strerror(errno));
			fprintf(stderr, "Velocity field dump failed.\n");
			break;
		}
		bytes_written += rc;
		bytes_left -= rc;
	} while (bytes_left > 0);
	close(fd);
	printf("Velocity field dumped to %s\n", filename);
	free(filename);
}

static int restore_velocity_field(char *filename, struct velocity_field *vf)
{
	int fd;
	int bytes_left = (int) sizeof(*vf);
	int bytes_read = 0;

	if (!filename || !restore_vf_data)
		return -1;

	printf("\n"); fflush(stdout);
	fd = open(filename, O_RDONLY);
	if (fd < 0) {
		fprintf(stderr, "Cannot open '%s' for reading: %s\n",
			filename, strerror(errno));
		return -1;
	}

	/* FIXME: this is quick and dirty, not endian clean, etc. */
	do {
		unsigned char *x = (unsigned char *) vf;
		int rc = read(fd, &x[bytes_read], bytes_left);
		if (rc < 0 && errno == EINTR)
			continue;
		if (rc < 0) {
			fprintf(stderr, "Error reading '%s': %s\n", filename, strerror(errno));
			fprintf(stderr, "Velocity field restoration failed.\n");
			close(fd);
			return -1;
		}
		bytes_read += rc;
		bytes_left -= rc;
	} while (bytes_left > 0);
	close(fd);
	printf("Velocity field restored from %s\n", filename);
	return 0;
}

/* move a particle according to velocity field at its current location */
static void move_particle(struct particle *p, struct velocity_field *vf, union vec3 *prev_particle_pos)
{
	struct fij fij;

	if (prev_particle_pos)
		*prev_particle_pos = p->pos;
	fij = xyz_to_fij(&p->pos, vfdim);
	p->fij = xyz_to_fij(&p->pos, DIM);
	vec3_add_self(&p->pos, &vf->v[fij.f][fij.i][fij.j]);
	vec3_normalize_self(&p->pos);
	vec3_mul_self(&p->pos, (float) vfdim / 2.0f);
}

static void *move_particles_thread_fn(void *info)
{
	struct movement_thread_info *thr = info;

	for (int i = 0; i < thr->iterations; i++)
		for (int j = thr->first_particle; j <= thr->last_particle; j++)
			move_particle(&thr->p[j], thr->vf, prev_particle_pos ? &prev_particle_pos[j] : NULL);
	return NULL;
}

static void create_move_particles_thread(struct particle *p, struct movement_thread_info *thr,
			struct velocity_field *vf, int thread_iterations)
{
	int rc;

	thr->vf = vf;
	thr->p = p;
	thr->iterations = thread_iterations;
	rc = create_thread(&thr->thread, move_particles_thread_fn, thr, "move-particles", 0);
	if (rc)
		fprintf(stderr, "%s: create_thread failed: %s\n",
				__func__, strerror(errno));
}

static void wait_for_movement_threads(struct movement_thread_info ti[], int nthreads)
{
	int i;
	void *status;

	for (i = 0; i < nthreads; i++) {
		int rc = pthread_join(ti[i].thread, &status);
		if (rc)
			fprintf(stderr, "%s: pthread_join failed: %s\n",
					__func__, strerror(errno));
	}
}

static void move_particles(int nthreads, struct particle *p, struct movement_thread_info ti[],
				struct velocity_field *vf, struct timing_data *timing,
				int thread_iterations)
{
	int t;

	gettimeofday(&timing->begin, NULL);
	for (t = 0; t < nthreads; t++)
		create_move_particles_thread(p, &ti[t], vf, thread_iterations);
	wait_for_movement_threads(ti, nthreads);
	gettimeofday(&timing->end, NULL);
	timing->elapsed += timeval_difference(timing->begin, timing->end);
}

struct image_thread_info {
	pthread_t thread;
	struct particle *p;
	int face;
	int nparticles;
};

static const int xo[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
static const int yo[] = { 0, 1, 1, 1, 0, -1, -1, -1 };

static void *update_output_image_thread_fn(void *info)
{
	struct image_thread_info *t = info;
	struct particle *p = t->p;
	int i, j;
	const int octaves = noise_levels - 1;
	/* 0.29289322 == 1 - 1/sqrt(2). */
	const float opacity_filter[] = {0.5, 0.29289322, 0.5, 0.29289322, 0.5, 0.29289322, 0.5, 0.29289322};

	if (!nofade)
		fade_out_background(t->face, &darkest_color);
	for (i = 0; i < t->nparticles; i++) {
		if (p[i].fij.f != t->face)
			continue;
		p[i].c.a = opacity;
		paint_particle(t->face, p[i].fij.i, p[i].fij.j, &p[i].c, octaves);
		if (large_pixels) {
			for (j = 0; j < 8; j++) {
				p[i].c.a = opacity * opacity_filter[j];
				paint_particle(t->face, p[i].fij.i + xo[j], p[i].fij.j + yo[j], &p[i].c, octaves);
			}
		}
	}
	return NULL;
}

static void update_output_images(int image_threads, struct particle p[], const int nparticles,
			struct timing_data *timing)
{
	struct image_thread_info t[image_threads];
	int i, rc;
	void *status;

	gettimeofday(&timing->begin, NULL);
	for (i = 0; i < image_threads; i++) {
		t[i].face = i;
		t[i].p = p;
		t[i].nparticles = nparticles;
		rc = create_thread(&t[i].thread, update_output_image_thread_fn, &t[i], "update-image", 0);
		if (rc)
			fprintf(stderr, "%s: create_thread failed: %s\n",
					__func__, strerror(errno));
	}

	for (i = 0; i < image_threads; i++) {
		int rc = pthread_join(t[i].thread, &status);
		if (rc)
			fprintf(stderr, "%s: pthread_join failed: %s\n",
					__func__, strerror(errno));
	}
	if (opacity > opacity_limit)
		opacity = 0.95 * opacity;
	gettimeofday(&timing->end, NULL);
	timing->elapsed += timeval_difference(timing->begin, timing->end);
}

static void maybe_save_equirectangular_images(char *output_file_prefix, int sequence_number,
			unsigned char *image[6], int has_alpha)
{
	char fname[PATH_MAX];
	int h, w, x, y;
	unsigned char *input_pixel, *output_pixel, *equirect;
	union vec3 p;
	struct fij fij;
	float latitude, longitude;
	int px, py;

	if (!export_equirect_image)
		return;

	if (sequence_number < 0)
		sprintf(fname, "%s-eqr.png", output_file_prefix);
	else
		sprintf(fname, "%s%04d-eqr.png", output_file_prefix, sequence_number);
	h = equirect_height;
	w = 2 * equirect_height;
	equirect = calloc(1, h * w * 4);
	memset(equirect, 255, h * w * 4);
	for (y = 0; y < h; y++) {
		latitude = ((float) y / (float) h) * M_PI - M_PI / 2;
		for (x = 0; x < w; x++) {
			output_pixel = equirect + 4 * (y * w  + x);
			longitude = ((float) x / (float) w) * 2.0 * M_PI;

			if (!vertical_bands) {
				p.v.x = cosf(longitude) * cosf(latitude);
				p.v.z = sinf(longitude) * cosf(latitude);
				p.v.y = sinf(latitude);
			} else {
				p.v.x = cosf(longitude) * cosf(latitude);
				p.v.y = sinf(longitude) * cosf(latitude);
				p.v.z = sinf(latitude);
			}
			fij = xyz_to_fij(&p, DIM);

			px = fij.i * (3 + !!has_alpha);
			py = fij.j * (3 + !!has_alpha) * XDIM;

			input_pixel = &image[fij.f][py + px];
			output_pixel[0] = input_pixel[0];
			output_pixel[1] = input_pixel[1];
			output_pixel[2] = input_pixel[2];
			output_pixel[3] = 255;
		}
	}
	if (png_utils_write_png_image(fname, equirect, w, h, 1, 0))
			fprintf(stderr, "Failed to write %s\n", fname);
}

static void save_output_images(char *output_file_prefix, int sequence_number, unsigned char *image[6], int has_alpha,
			struct timing_data *timing)
{
	int i;
	char fname[PATH_MAX];
	char *msg = "Saving Images";

	gettimeofday(&timing->begin, NULL);
	printf("%s", msg); fflush(stdout);
	for (i = 0; i < 6; i++) {
		if (sequence_number < 0)
			sprintf(fname, "%s%d.png", output_file_prefix, i);
		else
			sprintf(fname, "%s%04d-%d.png", output_file_prefix, i, sequence_number);
		if (png_utils_write_png_image(fname, image[i], DIM, DIM, has_alpha, 0))
			fprintf(stderr, "Failed to write %s\n", fname);
	}
	maybe_save_equirectangular_images(output_file_prefix, sequence_number, image, has_alpha);
	backspace(strlen(msg));
	printf("o");
	fflush(stdout);
	gettimeofday(&timing->end, NULL);
	timing->elapsed += timeval_difference(timing->begin, timing->end);
}

static void find_darkest_pixel(unsigned char *image, int w, int h,
				struct color *darkest_color, int bpr)
{
	int x, y, n;
	float r, g, b, min;

	if (use_hot_pink) {
		*darkest_color = hot_pink;
		return;
	}

	min = 3.0f;
	for (x = 0; x < w; x++) {
		for (y = 0; y < h; y++) {
			n = y * bpr + x * 3;
			r = (float) image[n + 0] / 255.0f;
			g = (float) image[n + 1] / 255.0f;
			b = (float) image[n + 2] / 255.0f;
			if (r + g + b < min) {
				darkest_color->r = r;
				darkest_color->g = g;
				darkest_color->b = b;
				min = r + g + b;
			}
		}
	}
	if (!cloudmode)
		darkest_color->a = 0.01;
	else
		darkest_color->a = 0.5;
	if (fade_rate > 0)
		darkest_color->a = fade_rate;
}

static char *load_image(const char *filename, int *w, int *h, int *a, int *bytes_per_row)
{
	char *i;
	char msg[100];

	i = png_utils_read_png_image(filename, 0, 0, 0, w, h, a, msg, sizeof(msg));
	if (!i) {
		fprintf(stderr, "%s: cannot load image: %s\n", filename, msg);
		exit(1);
	}
	*bytes_per_row = *w * (3 + !!(*a));
	/* align to 4 byte boundary */
	if (*bytes_per_row & 0x03)
		*bytes_per_row += 4 - (*bytes_per_row & 0x03);
	find_darkest_pixel((unsigned char *) i, *w, *h, &darkest_color, *bytes_per_row);
	return i;
}

void allocate_output_images(unsigned char *image[6], int bytes_per_pixel)
{
	int i;

	for (i = 0; i < 6; i++)
		image[i] = calloc(1, bytes_per_pixel * DIM * DIM);
}

static void free_output_images(unsigned char *image[6])
{
	int i;

	for (i = 0; i < 6; i++)
		if (image[i])
			free(image[i]);
}

static void usage(void)
{
	fprintf(stderr, "usage: gaseous-giganticus [-b bands] [-i inputfile] [-o outputfile]\n");
	fprintf(stderr, "       [-w w-offset] [-h] [-n] [-v velocity factor] [-B band-vel-factor]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "   -a, --pole-attenuation: attenuate band velocity near poles.  If\n");
	fprintf(stderr, "                 attenuation is zero, band velocity at poles will be\n");
	fprintf(stderr, "                 relatively high due to short distance around band.\n");
	fprintf(stderr, "                 Default is 0.5. Min is 0.0, max is 1.0.\n");
	fprintf(stderr, "   -b, --bands : Number of counter rotating bands.  Default is 6.0\n");
	fprintf(stderr, "   -B, --band-vel-factor: Multiply band velocity by this number when\n");
	fprintf(stderr, "                   computing velocity field.  Default is 2.9\n");
	fprintf(stderr, "   -c, --count : Number of iterations to run the simulation.\n");
	fprintf(stderr, "                 Default is 1000\n");
	fprintf(stderr, "   -C, --cloudmode: modulate image output by to produce clouds\n");
	fprintf(stderr, "   -d, --dump-velocity-field : dump velocity field data to specified file\n");
	fprintf(stderr, "                               (see -r option, below)\n");
	fprintf(stderr, "   -D, --faderate : Rate at which particles fade, default = 0.01\n");
	fprintf(stderr, "   --dump-flowmap : dump velocity field data to a series of\n");
	fprintf(stderr, "                    six RGB png files encoding x, y in R, G channels\n");
	fprintf(stderr, "   -e, --band-speed-power: Band speed is modulated by\n");
	fprintf(stderr, "	cos(K*latitude) ^ band-speed-power, where band-speed-power is an\n");
	fprintf(stderr, "	odd integer. By default, 1, higher values make the fast moving part\n");
	fprintf(stderr, "	of the bands narrow and the parts between wider.\n");
	fprintf(stderr, "   -E, --equirectangular image-height. Output an equirectangular image in\n");
	fprintf(stderr, "         addition to the usual cubemap images.  The image height must be\n");
	fprintf(stderr, "         an integer power of 2.\n");
	fprintf(stderr, "   -f, --fbm-falloff: Use specified falloff for FBM noise.  Default is 0.5\n");
	fprintf(stderr, "   -F, --vfdim: Set size of velocity field.  Default:2048. Min: 16. Max: 2048\n");
	fprintf(stderr, "   -g, --gain, 2nd and later octaves are multiplied by pow(fbm-falloff, (octave-1)*gain)\n");
	fprintf(stderr, "   -h, --hot-pink: Gradually fade pixels to hot pink.  This will allow\n");
	fprintf(stderr, "                   divergences in the velocity field to be clearly seen,\n");
	fprintf(stderr, "                   as pixels that contain no particles wil not be painted\n");
	fprintf(stderr, "                   and will become hot pink.\n");
	fprintf(stderr, "   -i, --input : Input image filename.  Must be RGB or RGBA png file.\n");
	fprintf(stderr, "   -I, --image-save-period: Interval of simulation iterations after which\n");
	fprintf(stderr, "         to output images.  Default is every 25 iterations\n");
	fprintf(stderr, "   -k, --cubemap: input 6 RGB or RGBA png cubemap images to initialize particle color\n");
	fprintf(stderr, "   -K, --cache-aware fraction:  initialize the first fraction of particles in a cache\n");
	fprintf(stderr, "         aware way (faster) and remaining particles with random locations.\n");
	fprintf(stderr, "         E.g. -K 0.5 will initialize 50%% of particles in a cache aware way and\n");
	fprintf(stderr, "         50%% randomly.\n");
	fprintf(stderr, "   -l, --large-pixels: particles are 3x3 pixels instead of 1 pixel.  Allows\n");
	fprintf(stderr, "                       for fewer particles and thus faster rendering, but this\n");
	fprintf(stderr, "                       WILL leave visible artifacts at cubemap face boundaries.\n");
	fprintf(stderr, "   -L, --noise-levels: Number of fractal noise levels to consider.  Default 4, max %d\n",
			max_noise_levels);
	fprintf(stderr, "   -m, --speed-multiplier:  band_speed_factor and velocity_factor are\n");
	fprintf(stderr, "         by this number.  It is a single option to affect both by a\n");
	fprintf(stderr, "         multiplier which is a bit easier than setting an absolute value\n");
	fprintf(stderr, "         for them individually.\n");
	fprintf(stderr, "   -n, --no-fade:  Do not fade the image at all, divergences will be hidden\n");
	fprintf(stderr, "   -N, --sequence-number:  Do not overwrite output images, instead embed a\n");
	fprintf(stderr, "                           sequence number in the filenames.\n");
	fprintf(stderr, "   -o, --output : Output image filename template.\n");
	fprintf(stderr, "               Example: 'out-' will produces 6 output files\n");
	fprintf(stderr, "               out-0.png, out-1.png, ..., out-5.png\n");
	fprintf(stderr, "   -O, --opacity: Specify minimum opacity of particles between 0 and 1.\n");
	fprintf(stderr, "                  default is 0.2\n");
	fprintf(stderr, "   -p, --particles: Use specified number of particles. Default is 8000000.\n");
	fprintf(stderr, "   -P, --plainmap  Do not use sinusoidal image mapping, instead repeat image\n");
	fprintf(stderr, "                   on six sides of a cubemap.\n");
	fprintf(stderr, "   -r, --restore-velocity-field: Instead of computing the velocity field from\n");
	fprintf(stderr, "                                 scratch, import a previously saved velocity\n");
	fprintf(stderr, "                                 field (see -d option, above)\n");
	fprintf(stderr, "   -R, --random: Random values are used for bands, band-vel-factor, velocity-factor\n"
			"                 noise-scale, and w-offset.  -S and -V options are also set.\n");
	fprintf(stderr, "   -s, --stripe: Begin with stripes from a vertical strip of input image\n");
	fprintf(stderr, "   -S, --sinusoidal: Use sinusoidal projection for input image\n");
	fprintf(stderr, "                 Note: sinusoidal is the default projection.\n");
	fprintf(stderr, "                 Note: --stripe and --sinusoidal are mutually exclusive\n");
	fprintf(stderr, "   -t, --threads: Use the specified number of CPU threads up to the\n");
	fprintf(stderr, "                   number of online CPUs\n");
	fprintf(stderr, "   -T, --thread-iterations: Number of iterations of thread movement particles\n");
	fprintf(stderr, "         before threads join. Default is 1. Higher values will reduce the number\n");
	fprintf(stderr, "         of times threads are created and destroyed.\n");
	fprintf(stderr, "   --trap-nans: trap divide by zero, overflow and invalid floating point\n");
	fprintf(stderr, "                exceptions\n");
	fprintf(stderr, "   -v, --velocity-factor: Multiply velocity field by this number when\n");
	fprintf(stderr, "                   moving particles.  Default is 1200.0\n");
	fprintf(stderr, "   -V, --vertical-bands:  Make bands rotate around Y axis instead of X\n");
	fprintf(stderr, "                          Also affects --stripe option.\n");
	fprintf(stderr, "   --vortex-band-threshold: controls distribution of vortices (see man page)\n");
	fprintf(stderr, "   --vortex-size: radius of vortices as a fraction of planet radius. Default 0.04\n");
	fprintf(stderr, "   --vortex-size-variance: Range of vortex sizes, default is plus or minus 0.02\n");
	fprintf(stderr, "   -w, --w-offset: w dimension offset in 4D open simplex noise field\n");
	fprintf(stderr, "                   Use -w to avoid (or obtain) repetitive results.\n");
	fprintf(stderr, "   -W, --wstep: w coordinate of noise field is incremented by specified\n");
	fprintf(stderr, "                amount periodically and velocity field is recalculated\n");
	fprintf(stderr, "   -x, --vortices: how many artificial circular vortices to add into the v-field\n");
	fprintf(stderr, "   -z, --noise-scale: default is %f\n", default_noise_scale);
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: As noise-scale increases, velocity-factor should generally decrease.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Examples:\n");
        fprintf(stderr, "\n");
	fprintf(stderr, "   ./gaseous-giganticus --noise-scale 2.5 --velocity-factor 1300 --bands 10\\\n"
			"        -i input_image.png -o output_image\n");
	fprintf(stderr, "   ./gaseous-giganticus -R -i input_image.png -o output_image\n");
	fprintf(stderr, "\n");
	exit(1);
}

#define VORTEX_SIZE_OPTION 1000
#define VORTEX_SIZE_VARIANCE_OPTION 1001
#define VORTEX_BAND_THRESHOLD_OPTION 1002
#define DUMP_FLOWMAP_OPTION 1003
#define TRAP_NANS_OPTION 1004

static struct option long_options[] = {
	{ "pole-attenuation", required_argument, NULL, 'a' },
	{ "bands", required_argument, NULL, 'b' },
	{ "count", required_argument, NULL, 'c' },
	{ "cloudmode", required_argument, NULL, 'C' },
	{ "dump-velocity-field", required_argument, NULL, 'd' },
	{ "faderate", required_argument, NULL, 'D' },
	{ "input", required_argument, NULL, 'i' },
	{ "image-save-period", required_argument, NULL, 'I' },
	{ "cubemap", required_argument, NULL, 'k' },
	{ "cache-aware", required_argument, NULL, 'K' },
	{ "large-pixels", no_argument, NULL, 'l' },
	{ "noise-levels", no_argument, NULL, 'L' },
	{ "output", required_argument, NULL, 'o' },
	{ "opacity", required_argument, NULL, 'O' },
	{ "w-offset", required_argument, NULL, 'w' },
	{ "fbm-falloff", required_argument, NULL, 'f' },
	{ "gain", required_argument, NULL, 'g' },
	{ "hot-pink", no_argument, NULL, 'h' },
	{ "help", no_argument, NULL, 'H' },
	{ "no-fade", no_argument, NULL, 'n' },
	{ "sequence-number", no_argument, NULL, 'N' },
	{ "velocity-factor", required_argument, NULL, 'v' },
	{ "vertical-bands", required_argument, NULL, 'V' },
	{ "band-vel-factor", required_argument, NULL, 'B' },
	{ "restore-velocity-field", required_argument, NULL, 'r' },
	{ "random", no_argument, NULL, 'R' },
	{ "stripe", no_argument, NULL, 's' },
	{ "sinusoidal", no_argument, NULL, 'S' },
	{ "threads", required_argument, NULL, 't' },
	{ "thread-iterations", required_argument, NULL, 'T' },
	{ "particles", required_argument, NULL, 'p' },
	{ "plainmap", required_argument, NULL, 'P' },
	{ "speed-multiplier", required_argument, NULL, 'm' },
	{ "magic-fluid-flow", no_argument, NULL, 'M'},
	{ "vfdim", required_argument, NULL, 'F' },
	{ "wstep", required_argument, NULL, 'W' },
	{ "wstep-period", required_argument, NULL, 'q' },
	{ "vortices", required_argument, NULL, 'x' },
	{ "vortex-band-threshold", required_argument, NULL, VORTEX_BAND_THRESHOLD_OPTION },
	{ "vortex-size", required_argument, NULL, VORTEX_SIZE_OPTION },
	{ "vortex-size-variance", required_argument, NULL, VORTEX_SIZE_VARIANCE_OPTION },
	{ "noise-scale", required_argument, NULL, 'z' },
	{ "band-speed-power", required_argument, NULL, 'e' },
	{ "equirectangular", required_argument, NULL, 'E' },
	{ "dump-flowmap", required_argument, NULL, DUMP_FLOWMAP_OPTION },
	{ "trap-nans", no_argument, NULL, TRAP_NANS_OPTION },
	{ 0, 0, 0, 0 },
};

static void process_float_option(char *option_name, char *option_value, float *value)
{
	float tmpf;

	if (sscanf(option_value, "%f", &tmpf) == 1) {
		*value = tmpf;
	} else {
		fprintf(stderr, "Bad %s option '%s'\n", option_name, option_value);
		usage();
	}
}

static void process_int_option(char *option_name, char *option_value, int *value)
{
	int tmp;

	if (sscanf(option_value, "%d", &tmp) == 1) {
		*value = tmp;
	} else {
		fprintf(stderr, "Bad %s option '%s'\n", option_name, option_value);
		usage();
	}
}

static void set_automatic_options(int random_mode)
{
	struct timeval tv;
	int has_vortices;

	if (!random_mode)
		return;

	/* set Sinusoidal mode */
	sinusoidal = 1;
	stripe = 0;
	/* set horizontal bands */
	vertical_bands = 0;
	gettimeofday(&tv, NULL);
	srand(tv.tv_usec);
	has_vortices = (rand() % 1000) < 500;
	noise_scale = ((float) rand() / (float) RAND_MAX) * 3.0 + 1.0;
	velocity_factor = ((float) rand() / (float) RAND_MAX) * 1500.0 + 500.0;
	band_speed_factor = ((float) rand() / (float) RAND_MAX) * 3.0;
	num_bands = ((float) rand() / (float) RAND_MAX) * 10.0 + 5.0;
	w_offset = ((float) rand() / (float) RAND_MAX) * 300.0;

	if (has_vortices) {
		vortex_band_threshold = ((float) rand() / (float) RAND_MAX) * 0.95 + 0.05;
		vortex_size_variance =  ((float) rand() / (float) RAND_MAX) * 0.04;
		vortex_size =  ((float) rand() / (float) RAND_MAX) * 0.04 + 0.02;
		nvortices = ((float) rand() / (float) RAND_MAX) * MAXVORTICES;
	}

	printf("\n");
	printf("Auto parameters:\n");
	printf("  --noise-scale %f\n", noise_scale);
	printf("  --velocity-factor %f\n", velocity_factor);
	printf("  --band-vel-factor %f\n", band_speed_factor);
	printf("  --num_bands %f\n", num_bands);
	printf("  --w-offset %f\n", w_offset);
	if (has_vortices) {
		printf("  --vortices %i\n", nvortices);
		printf("  --vortex-size %f\n", vortex_size);
		printf("  --vortex-size-variance %f\n", vortex_size_variance);
		printf("  --vortex-band-threshold %f\n", vortex_band_threshold);
	}
	printf("\n");
}

static int is_power_of_2(int n)
{
	int setbits = 0;
	int i;

	for (i = 0; (size_t) i < sizeof(int) * 8; i++) {
		if (n & 0x01)
			setbits++;
		n = n >> 1;
	}
	return (setbits == 1);
}

static void process_options(int argc, char *argv[])
{
	int i, c;

	output_file_prefix = default_output_file_prefix;
	input_file = default_input_file;

	while (1) {
		int option_index;
		c = getopt_long(argc, argv, "a:B:b:c:Cd:D:e:E:f:g:F:hHi:k:K:I:lL:nNm:o:O:p:PRr:sSt:T:Vv:w:W:x:z:",
				long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
		case 'a':
			process_float_option("pole-attenuation", optarg, &pole_attenuation);
			if (pole_attenuation < 0.0f)
				pole_attenuation = 0.0f;
			if (pole_attenuation > 1.0f)
				pole_attenuation = 1.0f;
			break;
		case 'B':
			process_float_option("band-vel-factor", optarg, &band_speed_factor);
			break;
		case 'b':
			process_float_option("num-bands", optarg, &num_bands);
			break;
		case 'c':
			process_int_option("count", optarg, &niterations);
			break;
		case 'C':
			cloudmode = 1;
			break;
		case 'd':
			vf_dump_file = optarg;
			restore_vf_data = 0;
			break;
		case 'D':
			process_float_option("faderate", optarg, &fade_rate);
			break;
		case 'e':
			process_int_option("band-speed-power", optarg, &band_speed_power);
			if (band_speed_power < 0 || (band_speed_power % 2) != 1) {
				fprintf(stderr, "band-speed-power must be an odd integer\n");
				usage();
			}
			band_speed_powerf = (float) band_speed_power;
			break;
		case 'E':
			process_int_option("equirectangular", optarg, &equirect_height);
			if (equirect_height < 1 || !is_power_of_2(equirect_height)) {
				fprintf(stderr, "Equirectangular height %d, must be integer power of 2.\n",
					equirect_height);
				usage();
			}
			export_equirect_image = 1;
			break;
		case 'r':
			vf_dump_file = optarg;
			restore_vf_data = 1;
			break;
		case 'R':
			random_mode = 1;
			break;
		case 'f':
			process_float_option("fbm-falloff", optarg, &fbm_falloff);
			ff[0] = 1.0;
			ff[1] = fbm_falloff;
			ff[2] = ff[1] * fbm_falloff;
			ff[3] = ff[2] * fbm_falloff;
			ff[4] = ff[3] * fbm_falloff;
			ff[5] = ff[4] * fbm_falloff;
			ff[6] = ff[5] * fbm_falloff;
			break;
		case 'g':
			process_float_option("gain", optarg, &gain);
			ff[0] = 1.0;
			ff[1] = fbm_falloff; /* same as pow(fbm_fallof, 0.0 * gain); */
			ff[2] = pow(fbm_falloff, 1.0 * gain);
			ff[3] = pow(fbm_falloff, 2.0 * gain);
			ff[4] = pow(fbm_falloff, 3.0 * gain);
			ff[5] = pow(fbm_falloff, 4.0 * gain);
			ff[6] = pow(fbm_falloff, 5.0 * gain);
			break;
		case 'F':
			process_int_option("vfdim", optarg, &vfdim);
			if (vfdim < 16) {
				vfdim = 16;
				fprintf(stderr, "Bad value of vfdim specified, using 16.\n");
			}
			if (vfdim > 2048) {
				vfdim = 2048;
				fprintf(stderr, "Bad value of vfdim specified, using 2048.\n");
			}
			break;
		case 'h':
			use_hot_pink = 1;
			break;
		case 'H':
			usage(); /* does not return */
			break;
		case 'i':
			input_file = optarg;
			break;
		case 'I':
			process_int_option("image-save-period", optarg, &image_save_period);
			break;
		case 'k':
			cubemap_prefix = optarg;
			break;
		case 'K':
			process_float_option("cache-aware", optarg, &cache_aware);
			if (cache_aware < 0.0 || cache_aware > 1.0) {
				fprintf(stderr, "cache_aware (-K) option value must be between 0.0 and 1.0\n");
				usage();
			}
			break;
		case 'm':
			process_float_option("speed-multiplier", optarg, &speed_multiplier);
			break;
		case 'M':
			magic_fluid_flow = 1;
			break;
		case 'n':
			nofade = 1;
			break;
		case 'N':
			save_texture_sequence = 1;
			break;
		case 'l':
			large_pixels = 1;
			break;
		case 'L':
			process_int_option("noise-levels", optarg, &noise_levels);
			if (noise_levels < 1) {
				fprintf(stderr, "Noise levels value of %d is out of range, using 1\n", noise_levels);
				noise_levels = 1;
			} else if (noise_levels > max_noise_levels) {
				fprintf(stderr, "Noise levels value of %d is out of range, using %d\n",
					noise_levels, max_noise_levels);
				noise_levels = max_noise_levels;
			}
			break;
		case 'o':
			output_file_prefix = optarg;
			break;
		case 'O':
			process_float_option("opacity", optarg, &opacity_limit);
			if (opacity_limit < 0.0)
				opacity_limit = 0.0;
			if (opacity_limit > 1.0)
				opacity_limit = 1.0;
			break;
		case 'p':
			process_int_option("particles", optarg, &particle_count);
			break;
		case 'P': /* plain mapping of image to 6 sides of cubemap */
			stripe = 0;
			sinusoidal = 0;
			break;
		case 's':
			stripe = 1;
			sinusoidal = 0;
			break;
		case 'S':
			sinusoidal = 1;
			stripe = 0;
			break;
		case 't':
			process_int_option("threads", optarg, &user_threads);
			break;
		case 'T':
			process_int_option("thread-iterations", optarg, &thread_iterations);
			break;
		case 'w':
			process_float_option("w-offset", optarg, &w_offset);
			break;
		case 'v':
			process_float_option("velocity-factor", optarg, &velocity_factor);
			break;
		case 'V':
			vertical_bands = 0;
			break;
		case 'W':
			process_float_option("wstep", optarg, &wstep);
			use_wstep = 1;
			break;
		case 'q': /* running out of letters, so, 'q'. */
			process_int_option("wstep-period", optarg, &wstep_period);
			break;
		case 'x':
			process_int_option("vortices", optarg, &nvortices);
			if (nvortices > MAXVORTICES) {
				nvortices = MAXVORTICES;
				fprintf(stderr, "Requested %d vortices, maximum is %d\n", nvortices, MAXVORTICES);
			}
			break;
		case 'z':
			process_float_option("noise-scale", optarg, &noise_scale);
			break;
		case VORTEX_SIZE_OPTION:
			process_float_option("vortex-size", optarg, &vortex_size);
			break;
		case VORTEX_SIZE_VARIANCE_OPTION:
			process_float_option("vortex-size", optarg, &vortex_size_variance);
			break;
		case VORTEX_BAND_THRESHOLD_OPTION:
			process_float_option("vortex-band-threshold", optarg, &vortex_band_threshold);
			if (vortex_band_threshold < 0.05)
				vortex_band_threshold = 0.05;
			if (vortex_band_threshold > 1.0)
				vortex_band_threshold = 1.0;
			break;
		case DUMP_FLOWMAP_OPTION:
			flowmap_dump_file = optarg;
			break;
		case TRAP_NANS_OPTION:
			feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
			break;
		default:
			fprintf(stderr, "unknown option '%s'\n",
				option_index > 0 && option_index - 1 < argc &&
				argv[option_index - 1] ? argv[option_index - 1] : "(none)");
			usage();
		}
	}

	/* Scale so that vfdim doesn't change the effect of these */
	band_speed_factor = ((float) vfdim / 2048.0) * band_speed_factor * speed_multiplier;
	velocity_factor = ((float) vfdim / 2048.0) * velocity_factor * speed_multiplier;

	for (i = noise_levels; i < max_noise_levels; i++)
		ff[i] = 0.0;

	return;
}

static void copy_cubemap_face_to_output_face(unsigned char *image, int face, int width, int height, int has_alpha)
{
	int x, y, ix, iy, ip, op;
	unsigned char *input_pixel, *output_pixel;

	for (y = 0; y < DIM; y++) {
		for (x = 0; x < DIM; x++) {
			ix = (x * width) / DIM;
			iy = (y * height) / DIM;
			ip = y * DIM + x;
			op = iy * width + ix;
			input_pixel = &image[ip * (3 + !!has_alpha)];
			output_pixel = &output_image[face][op * 4];
			output_pixel[0] = input_pixel[0];
			output_pixel[1] = input_pixel[1];
			output_pixel[2] = input_pixel[2];
			output_pixel[3] = 255;
		}
	}
}

static void load_cubemap_images(const char *prefix)
{
	char filename[PATH_MAX];
	int i;

	for (i = 0; i < 6; i++) {
		sprintf(filename, "%s%d.png", prefix, i);
		cubemap_image[i] = (unsigned char *) load_image(filename, &start_image_width, &start_image_height,
						&start_image_has_alpha, &start_image_bytes_per_row);
	}
	/* Copy the cubemap images directly to the output images as a starting point. */
	for (i = 0; i < 6; i++)
		copy_cubemap_face_to_output_face(cubemap_image[i], i,
				start_image_width, start_image_height, start_image_has_alpha);
}

/* Return a value between -1.0 and 1.0 with a bias towards zero.
 * Not a normal distribution, but sorta close enough for our purposes.
 */
static float random_squared(void)
{
	float x, y;

	x = (2.0f * (float) rand() / (float) RAND_MAX) - 1.0f;
	y = (2.0f * (float) rand() / (float) RAND_MAX) - 1.0f;
	return x * y;
}

static void create_vortex(int i)
{
	const union vec3 right_at_ya = { { 0.0f, 0.0f, 1.0f } };
	const union vec3 up = { { 0.0f, 1.0f, 0.0f } };
	float angle, band_speed;


	if (num_bands > 0) {
		/* I don't think this is actually correct.  I think
		 * I need to sample the velocity field around the vortex
		 * to determine the already present rotation and make the
		 * vortex rotate in the same direction... this doesn't
		 * quite do that.  And the velocity field is not built yet
		 * when this gets called (that might be change-able though.).
		 *
		 * A better method might be to pre-sample the noise field, then
		 * disturb the sampled noise gradient to make the vortices prior to
		 * constructing the velocity field.
		 */
		do {
			random_point_on_sphere(1.0, &vort[i].p.v.x, &vort[i].p.v.y, &vort[i].p.v.z);
			vort[i].r = vortex_size + random_squared() * vortex_size_variance;
			if (vertical_bands)
				angle = asinf(vort[i].p.v.z);
			else
				angle = asinf(vort[i].p.v.y);
			band_speed = calculate_band_speed(angle);
		} while (fabs(band_speed) > vortex_band_threshold * band_speed_factor &&
			fabs(angle) > 15.0 * M_PI / 180.0); /* exclude vortice within 15 deg of poles */
		if (calculate_band_speed(angle + 0.05) < calculate_band_speed(angle - 0.05))
			vort[i].angular_vel = 2.5;
		else
			vort[i].angular_vel = -2.5;
	} else {
		random_point_on_sphere(1.0, &vort[i].p.v.x, &vort[i].p.v.y, &vort[i].p.v.z);
		vort[i].r = vortex_size + random_squared() * vortex_size_variance;
		/* This will be right half the time :/ */
		if (rand() > (RAND_MAX >> 1))
			vort[i].angular_vel = 2.5;
		else
			vort[i].angular_vel = -2.5;
	}
	quat_from_u2v(&vort[i].q, &vort[i].p, &right_at_ya, &up);
	quat_from_u2v(&vort[i].reverseq, &right_at_ya, &vort[i].p, &up);
}

static void create_vortices()
{
	int i;

	for (i = 0; i < nvortices; i++)
		create_vortex(i);
}

#ifdef _WIN32
#include <windows.h>
static int get_num_cpus()
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
}
#else
#ifdef unix
static int get_num_cpus()
{
#ifdef _SC_NPROCESSORS_ONLN
	return sysconf(_SC_NPROCESSORS_ONLN);
#else
	return 1;
#warning "This variant of unix does not define _SC_NPROCESSORS_ONLN. Assuming 1 CPU. Performance will not be good."
#endif
}
#else
#warning "get_num_cpus() not found. Assuming 1 CPU. Performance will not be good."
static int get_num_cpus()
{
	return 1;
}
#endif
#endif

static double compute_flowmap_scaling_factor(int nparticles)
{
	int i;
	float maxv = 0.0;
	union vec2 p1, p2, pv;
	union vec3 fij1, fij2;
	float vm2;

	/* Project particle movement onto cubemap faces to get 2d velocity vectors
	 * and find the largest such velocity. Return 1.0 / sqrt(largest velocity)
	 * In this way, we can scale all the velocities such that the magnitudes
	 * are between 0 and 1.
	 */
	for (i = 0; i < nparticles; i++) {
		fij1 = xyz_to_fij_float(&prev_particle_pos[i], DIM, -1);
		fij2 = xyz_to_fij_float(&particle[i].pos, DIM, (int) fij1.v.z);
		if (fij1.v.z != fij2.v.z) /* Do not count particles that cross cubemap face boundaries */
			continue;
		p1.v.x = fij1.v.x;
		p1.v.y = fij1.v.y;
		p2.v.x = fij2.v.x;
		p2.v.y = fij2.v.y;
		pv.v.x = p2.v.x - p1.v.x;
		pv.v.y = p2.v.y - p1.v.y;
		vm2 = pv.v.x * pv.v.x + pv.v.y * pv.v.y;
#if 0
		/* Ignore very high velocities */
		if (vm2 > 25.0)
			continue;
#endif
		if (vm2 > maxv)
			maxv = vm2;
	}
	return 1.0 / sqrt(maxv);
}

static void encode_vec2_in_red_green(unsigned char *image, int i, int j, union vec2 *v)
{
	unsigned char *pixel;

	float x = (0.5 * v->v.x + 0.5);
	float y = (0.5 * v->v.y + 0.5);

	int p = j * DIM + i;
	pixel = &image[p * (3 + flowmap_has_alpha)];
	pixel[0] = float2int(255.0f * x) & 0xff;
	pixel[1] = float2int(255.0f * y) & 0xff;
	pixel[2] = 0;
	if (flowmap_has_alpha)
		pixel[3] = 255;
}

/* Compute a flow map in which particle velocity is represented as 6 images of
 * a cubemap with x and y velocities encoded in R and G channels
 */
static void maybe_compute_flowmap(int nparticles)
{
	static float scaling_factor = 1.1e6;
	union vec2 pv;
	union vec3 p1, p2;
	int p, f, i, j;

	if (!flowmap_dump_file)
		return;

	if (scaling_factor > 1e6)
		/* Only need to calculate this once as it is nearly the same every time. */
		scaling_factor = compute_flowmap_scaling_factor(nparticles);

	for (p = 0; p < nparticles; p++) {
		p1 = xyz_to_fij_float(&prev_particle_pos[p], DIM, -1);
		p2 = xyz_to_fij_float(&particle[p].pos, DIM, (int) p1.v.z);
		if (p1.v.z != p2.v.z) /* If the particle crossed cubemap face boundary, skip it */
			continue;
		pv.v.x = p2.v.x - p1.v.x;
		pv.v.y = p2.v.y - p1.v.y;
		pv.v.x *= scaling_factor;
		pv.v.y *= scaling_factor;

		i = float2int(p1.v.x);
		j = float2int(p1.v.y);
		f = float2int(p1.v.z);

		encode_vec2_in_red_green(flowmap_image[f], i, j, &pv);
	}
}

static void advect_velocity_field(void)
{
	/* TODO: simulate fluid momentum via lagrangian advection of velocity field.
	 * See: https://www.karlsims.com/fluid-flow.html
	 * The basic idea is as follows:
	 * For each velocity vector, look backwards (negate the vector) and find
	 * where that spot is.  Then use bilinear interpolation of the surrounding
	 * 4 velocity vectors to figure what the velocity vector at that point was,
	 * and that becomes the new velocity vector.
	 *
	 * There will need to be additional steps of flattening the local region
	 * to a plane, doing the advection step as above, then projecting back onto
	 * the sphere.
	 *
	 * This advection step should be "embarassingly parallel" because the new
	 * velocity field (vf) depends only on the old one (old_vf) which is not
	 * mutated.
	 */
}

static void remove_divergences(void)
{
	/* TODO: figure out how to remove divergences introduced by advection step
	 * from velocity field. See: https://www.karlsims.com/fluid-flow.html
	 *
	 * "... the flow field is repeatedly incremented by the gradient of its divergence.
	 * A gradient is the slope or rate of change of a value across the grid."
	 *
	 * I do not really have a good understanding of how to do this step.
	 */
}

static void do_fluid_dynamics(void)
{
	if (!magic_fluid_flow)
		return;
	memcpy(old_vf, vf, sizeof(struct velocity_field));
	advect_velocity_field();
	remove_divergences();
}

static void print_timing_info(int i, int niterations, struct timing_data *move, struct timing_data *image,
				struct timing_data *png, struct timing_data *prog)
{
	if ((i % 50) == 0) {
		gettimeofday(&prog->end, NULL);
		prog->elapsed = timeval_difference(prog->begin, prog->end);
		if (i != 0)
			printf(" Part mvmt:%g ms  Render:%g ms  PNG:%g ms Elapsed:%g ms\n",
				move->elapsed, image->elapsed, png->elapsed, prog->elapsed);
		else
			printf("\n");
		printf("%5d / %5d ", i, niterations);
	} else {
		printf(".");
	}
	fflush(stdout);
}

int main(int argc, char *argv[])
{
	int i, t;
	struct movement_thread_info *ti;
	int last_imaged_iteration = -1;
	int sequence_number = -1;
	particle_count = NPARTICLES;
	struct timing_data movetime, imagetime, pngtime, progtime;

	movetime.elapsed = 0.0;
	imagetime.elapsed = 0.0;
	pngtime.elapsed = 0;
	progtime.elapsed = 0;
	gettimeofday(&progtime.begin, NULL);

	/* Allocate a ton of memory. We allocate these rather than
	 * declaring them statically because if DIM is too large, the
	 * static data will not fit into the .bss (see https://en.wikipedia.org/wiki/.bss)
	 */
	vf = malloc(sizeof(struct velocity_field));
	if (magic_fluid_flow)
		old_vf = malloc(sizeof(struct velocity_field));
	else
		old_vf = NULL;

	open_simplex_noise(3141592, &ctx);

	setlocale(LC_ALL, "");
	noise_scale = default_noise_scale;

	process_options(argc, argv);
	set_automatic_options(random_mode);
	create_vortices();

	check_vf_dump_file(vf_dump_file);

	int num_online_cpus = get_num_cpus();
	if (num_online_cpus > 0)
		nthreads = num_online_cpus;
	if (user_threads > 0 && user_threads < num_online_cpus)
		nthreads = user_threads;
	printf("Using %d threads for particle motion\n", nthreads);

	int tparticles = particle_count / nthreads;

	ti = malloc(sizeof(*ti) * nthreads);

	for (t = 0; t < nthreads; t++) {
		ti[t].first_particle = t * tparticles;
		ti[t].last_particle = (t + 1) * tparticles - 1;
	}
	/* last thread picks up extra (particle_count - (nthreads * tparticles)) */
	t = nthreads - 1;
	ti[t].last_particle = particle_count - 1;

	printf("Allocating output image space, 6x %dx%d\n", DIM, DIM);
	allocate_output_images(output_image, 4);
	if (flowmap_dump_file)
		allocate_output_images(flowmap_image, 3 + flowmap_has_alpha);
	start_image_has_alpha = 0;
	if (cubemap_prefix) {
		printf("Loading cubemap images\n");
		load_cubemap_images(cubemap_prefix);
	} else {
		printf("Loading image\n");
		start_image = load_image(input_file, &start_image_width, &start_image_height,
						&start_image_has_alpha, &start_image_bytes_per_row);
	}
#if 0	
	write_png_image("blah.png", (unsigned char *) start_image, start_image_width,
			start_image_height, start_image_has_alpha);
#endif
	printf("width, height, bytes per row = %d,%d,%d\n",
			start_image_width, start_image_height, start_image_bytes_per_row);
	init_particles(&particle, &prev_particle_pos, particle_count);
	if (restore_velocity_field(vf_dump_file, vf))
		update_velocity_field(vf, w_offset, &use_wstep);
	dump_velocity_field(vf_dump_file, vf, use_wstep);

	for (i = 0; i < niterations; i++) {
		print_timing_info(i, niterations, &movetime, &imagetime, &pngtime, &progtime);
		move_particles(nthreads, particle, ti, vf, &movetime, thread_iterations);
		do_fluid_dynamics();
		update_output_images(image_threads, particle, particle_count, &imagetime);
		if ((i % image_save_period) == 0) {
			sequence_number += save_texture_sequence;
			save_output_images(output_file_prefix, sequence_number, output_image, 1, &pngtime);
			last_imaged_iteration = i;
		}
		if (use_wstep && (i % wstep_period == 0)) {
			w_offset += wstep;
			update_velocity_field(vf, w_offset, &use_wstep);
			dump_velocity_field(vf_dump_file, vf, use_wstep);
		}
		if (i < flowmap_impressions)
			maybe_compute_flowmap(particle_count);
		if (i == flowmap_impressions && flowmap_dump_file)
			save_output_images(flowmap_dump_file, -1, flowmap_image, flowmap_has_alpha, &pngtime);
	}
	if (last_imaged_iteration != i - 1) {
		sequence_number += save_texture_sequence;
		save_output_images(output_file_prefix, sequence_number, output_image, 1, &pngtime);
	}
	printf("\n%5d / %5d -- done.\n", i, niterations);
	open_simplex_noise_free(ctx);

	/* free a ton of memory */
	free(vf);
	if (old_vf)
		free(old_vf);
	free(particle);
	free(ti);
	free_output_images(output_image);
	free_output_images(flowmap_image);
	if (prev_particle_pos)
		free(prev_particle_pos);

	gettimeofday(&progtime.end, NULL);
	progtime.elapsed = timeval_difference(progtime.begin, progtime.end);
	printf("\nElapsed time: %g ms\n", progtime.elapsed);

	return 0;
}

