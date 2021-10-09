# SolidPrep - Experimental method of structural design automation

SolidPrep ("Solid Preprocessor") is an experimental method of automating the 
design of continuum structures, based on the use of continuum beams with 
variable cross-section and taking advantage of topology optimization.

## Installation

Dependencies:
- [OpenCascade](https://git.dev.opencascade.org/gitweb/?p=occt.git;a=summary)
- [LAPACK, LAPACKE, CBLAS](http://www.netlib.org/lapack/)
- [NLopt](https://nlopt.readthedocs.io/en/latest/)
- [Gmsh](https://gmsh.info/)

Installation consists on installing the dependencies somewhere CMake can find 
them (which should be done automatically when building them with CMake or using 
a package manager such as `apt`). After that, the program can be built as usual:

```bash
mkdir build
cd build
cmake ..
make
```

The program has not been tested on Windows nor MacOS. It should work, but it is 
possible that some modifications may be needed.

## Using SolidPrep

To run an optimization problem, just execute the program in a terminal 
environment with the path to the problem's file as argument, for example:

```
./build/SolidPrep ./examples/cantilever.json
```

The JSON files defines the analysis to be done as well as all of the 
characteristics of the problem (except for the geometry of the ground structure, 
which should be supplied as a STEP file).

See the `examples` folder for more information.
