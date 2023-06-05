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
- [RapidJSON](https://github.com/Tencent/rapidjson) (included as submodule)
- MPI (currently only [MPICH](https://www.mpich.org/))
- [MUMPS](https://graal.ens-lyon.fr/MUMPS/index.php) ([CMake version](https://github.com/scivision/mumps))
- [Eigen](https://eigen.tuxfamily.org/)
- [PETSc](https://petsc.org/release/)
- [NVIDIA HPC](https://developer.nvidia.com/hpc-sdk) (optional)

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
