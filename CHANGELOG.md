# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

## [v0.2.0]
### Added
- `ElementCommon` family of classes, to enhance code reuse and extension among
  element classes.
- `Geometry` class, replacing `GroundStructure`.
- Support for multiple geometries.
- Support for multiple materials.
- Q4, Q4S, TET4, H8 elements.
- Iterative solvers (gradient descent and PCG) (not recommended).
- MUMPS support.
- PETSc support (CPU + CUDA).
- Eigen3 support.
- Gauss-Legendre integration (using GSL).
- Improved results display using SPView.
- OpenMP support.
- Support for combining different objective and restriction functions.
- Support for different density filters.
- Support for different projection approaches.
- Support for generating fields using the generated meshes.
- Support for non-homogeneous materials.
- Support for different boundary conditions (beyond force and support).
- GlobalStiffnessMatrix base class.
- General solvers.
- Experimental features.

### Fixed
- Many small fixes.
- Improved error messages in ProjectData (more work is still necessary).
- Display error message on missing STEP file.
- Some problems were obsoleted due to rewrites.

### Changed
- Moved project loading from header to source file.
- Remove `results` variable from Node.
- Fully refactor element classes.
- Fully refactor finite element classes.
- Fully refactor meshing classes.
- Improve OpenMP approach to MMA.
- Deprecate `topology_optimization` directory.
- Many optimizations.

### Removed
- `GroundStructure`class.

## [v0.1.2]
### Added
- MMASolver from [jdumas/mma](https://github.com/jdumas/mma).
- Views for stress components for `fea_only` and `beams_only`.
- `get_stress_tensor()` for elements.
- Examples for v0.1.2 (`beams_only`).

### Changed
- Global stiffness matrix now uses minimal K instead of minimal density.
- Organize examples by version.
- Improved resizing method so it works with trussed and indeterminate structures.
- Somewhat improved speed of resizing.
- Make OneLab default to dark mode.
- Make MinimalVolume and MinimalCompliance use MMA.
- Show sizing time when using `beams_only` mode.
- Change font renderer for OneLab.

### Fixed
- Resizing of indeterminate beams.
- Compatibility with GCC 11 and OpenCascade 7.6.1.
- Resizing creating stray circles.
- Make OneLab more responsive after optimization end.
- Make density visualization grayscale.
- Stop flickering when updating multiple views.

### Removed
- Unused/unusable functions from StandardSizing.
- `calculate_stress()` from FiniteElement.
- Unused functions from MeshElement.
- NLopt as dependency

## [v0.1.1]
### Added
- Minimal compliance support.
- GCMMASolver class from [jdumas/mma](https://github.com/jdumas/mma).

### Changed
- Make MinimalVolume work correctly.
- Replace NLopt for GCMMASolver in MinimalVolume.
- Fix force vector generation.
- Make DirectSolver slightly faster.
- Fix StandardSizing generating misplaced circles (in some cases).
- Improve visualization of results from optimization and finite element analysis.

## [v0.1.0]
### Added
- Implementation of pre-sizing algorithm (StandardSizing class)
- Visibility graph pathfinding
- Gmsh-based meshing
- Lapack-based linear FEM direct solver
- Minimal volume topology optimization with global stress based on
LE et al, 2009
- TRI3 elements (implementation not fully tested)
- GT9 elements (YUQIU and YIN, 1994)
- Example problems
- Support for isotropic elastic elements
- Support for orthotropic elastic elements (not fully tested)
- CSG-based model handling (including loading and saving STEP files)
- Onelab-based GUI (Onelab being a part of the Gmsh library)

### Deprecated
- BeamSizing class (not fully implemented)
- MeshlessAStar (Vgraph is more efficient)
- Element-based sizing (present in StandardSizing class, does not give usable
results)
