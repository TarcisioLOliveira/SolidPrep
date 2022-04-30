# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- MMASolver from [jdumas/mma](https://github.com/jdumas/mma).

### Changed
- Global stiffness matrix now uses minimal K instead of minimal density.
- Organize examples by version.
- Improved resizing method so it works with trussed structures.
- Somewhat improved speed of resizing.
- Make OneLab default to dark mode.

### Fixed
- Resizing of indeterminate beams.
- Compatibility with GCC 11 and OpenCascade 7.6.1.
- Resizing creating stray circles.
- Make OneLab more responsive after optimization end.

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
