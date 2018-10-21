# Gaseous Giganticus 

This program is used to create procedurally generated gas giant planet cubemap
textures for the game Space Nerds In Space.

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

### Prerequisites

1. libpng

### Coding style

1. Run "git diff | ./checkpatch.pl -".  It will tell you what you did wrong.

TLDR: Tabs, not spaces, snake_case, not CamelCase.

## Contributing

Please read CONTRIBUTING.md

## Authors

* **Stephen M. Cameron**
* **Tobias Simon** (initial author of quaternion library)
* **Jeremy Van Grinsven** (many additions to quaternion library)

## License

GNU GPL v. 2

