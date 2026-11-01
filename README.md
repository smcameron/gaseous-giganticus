# Gaseous Giganticus 

This program is used to create procedurally generated gas giant planet cubemap
textures for the game Space Nerds In Space. Of course it can be used to create
assets for other games,
e.g. [Kerbal Space Program](https://forum.kerbalspaceprogram.com/topic/165285-planetary-texturing-guide-repository/#comment-3735392).

# Sample Output

The output of gaseous-giganticus consists of six square images that can be used as a
cubemap to texture a sphere.  Here the output is viewed with "mesh_viewer", from the
space-nerds-in-space repository.

![sample-gg-output-1.jpg](sample-gg-output-1.jpg?raw=true "gaseous-giganticus output")
![sample-gg-output-2.jpg](sample-gg-output-2.jpg?raw=true "gaseous-giganticus output")

## Getting Started

1. clone the project
2. Type "make"
3. run "./gaseous-giganticus --help" and do what it says.
4. See also: output of "nroff -man < gaseous-giganticus.1 | more"

The program requires one PNG input image, which should ideally be about 200 pixels wide
and about 1200 pixels tall, and blurred.  It outputs 6 square images that can be thought
of as the faces of a cube, which you can imagine being overinflated until it is a sphere.
The layout of those six output images is as follows:

```
        +------+
        |  4   |
        |      |
        +------+------+------+------+
        |  0   |  1   |  2   |  3   |
        |      |      |      |      |
        +------+------+------+------+
        |  5   |
        |      |
        +------+

```

If you wish to watch the progress of the program while it works, you'll need to compile and run
"mesh_viewer", which is part of [Space Nerds In Space](https://github.com/smcameron/space-nerds-in-space).

[Here is a video demoing the gaseous-giganticus and mesh_viewer](https://www.youtube.com/watch?v=8nx5yPpQh2M).
This video is a bit old, and I should probably make a new one, but it will have to suffice for now.

### Prerequisites

1. libpng

### Coding style

1. Run "git diff | ./checkpatch.pl -".  It will tell you what you did wrong.

TLDR: Tabs, not spaces, snake_case, not CamelCase.

## Contributing

To gain an understanding of how the program works, see [this slideshow](https://smcameron.github.io/space-nerds-in-space/gaseous-giganticus-slides/slideshow.html#1).
Use arrow keys to navigate.  That slideshow does not work well on mobile, so find a real computer.

For a long time (from 2014 to 2026) I thought it would be cool to add a fluid dynamics simulation to
gaseous-giganticus, and thought that this would likely produce better results.  When AI coding became
pretty strong, I attempted a few times to get it to add fluid dynamics without success.  But then, I
instead described to Claude Opus 5.5 what gaseous-giganticus does, and asked it to write a new program
that does that, and also fluid dynamics, making the program from scratch.  This worked, and produced
[gaseous-gusty-ball](https://github.com/smcameron/gaseous-gusty-ball).  However, adding fluid dynamics
did not dramatically improve the look of the output.  Feel free to try it out though, I don't know if
it's worse or better or about the same.  In any case, I'm no longer that interested in adding fluid
dynamics to gaseous-giganticus.

Please read CONTRIBUTING.md

## Authors

* **Stephen M. Cameron**
* **Tobias Simon** (initial author of quaternion library)
* **Jeremy Van Grinsven** (many additions to quaternion library)

## License

GNU GPL v. 2

## Notes

On Sunday, October 21, 2018, this code was split off from the Space
Nerds In Space repository at https://github.com/smcameron/space-nerds-in-space
where the code was originally developed.

An effort was made to preserve the history (git log, etc), and this
effort was 99.9% successful, but there are a few small differences from
the original history. Most significantly, the Makefile does not exist for
most of the history in this repository because the Makefile used in
Space Nerds In Space would have been broken anyway and would have
contained tons of irrelevant cruft. The Makefile in this repository was
added only after the entire gaseous-giganticus history was imported.
Less significantly, a few commits were made out-of-order, and although
the original dates were preserved in the log, checking out a particular
sha you may find some small differences if you compare with what was in
the Space Nerds In Space repository at a corresponding time. At the time
the Makefile was added to this repository, all source files present here
were identical to those also present in the Space Nerds In Space
repository, meaning those small differences were resolved eventually.

