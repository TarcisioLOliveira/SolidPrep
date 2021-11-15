# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
