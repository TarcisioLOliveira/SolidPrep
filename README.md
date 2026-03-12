# SolidPrep - Structural design automation and optimization

SolidPrep ("Solid Preprocessor") is a collection of methods developed to 
automate and optimize the design of continuum structures. It is currently 
under development, but is usable.

It serves mainly as research software, and contains experimental methods 
alongside established methods. It is designed to be extensible and flexible, in 
order to allow new methods to be added and compared to established ones.

Currently features:
- Gmsh/OneLab meshing and GUI
- JSON-based problem specification file
- Topology optimization:
    - SIMP
    - Convolution density filtering
    - Compliance and volume minimization
    - MMA and GCMMA
- Beam-based design automation (experimental)
    - Fast design of continuum structures based on methods used for sizing beams
    - Uses OpenCascade for solid geometry

## Installation

Dependencies:
- [OpenCascade](https://git.dev.opencascade.org/gitweb/?p=occt.git;a=summary)
- [LAPACK, LAPACKE, CBLAS](http://www.netlib.org/lapack/)
- [Gmsh](https://gmsh.info/)
- [VTK](https://vtk.org/) (Dependency for OpenCascade, may not be bundled with it)
- [JsonCpp](https://github.com/open-source-parsers/jsoncpp)
- MPI (currently only [MPICH](https://www.mpich.org/))
- [MUMPS](https://graal.ens-lyon.fr/MUMPS/index.php) ([CMake version](https://github.com/scivision/mumps))
- [Eigen](https://eigen.tuxfamily.org/)
- [PETSc](https://petsc.org/release/)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html)
- [TEEM](https://teem.sourceforge.net/download/index.html) (for NRRD files)
- [Catch2 v2](https://github.com/catchorg/Catch2/tree/v2.x) (for tests)

Also includes code from [jdumas/mma](https://github.com/jdumas/mma).

Installation consists on installing the dependencies somewhere CMake can find 
them (which should be done automatically when building them with CMake or using
a Linux package manager). Then, load submodule with:

```bash
git submodule update --init
```

After that, the program can be built as usual:

```bash
mkdir build
cd build
cmake ..
make
```

Currently only works on Linux due to some technical matters. PETSc and MUMPS 
apparently work better on Linux too, so it is the better choice in this case.


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
