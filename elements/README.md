# Element creation scripts

Python scripts designed to aid the creation of element classes. They depend
on `numpy` and `sympy`.

The scripts are used mainly to calculate the element-dependent matrices:
- the stiffness matrix (`k`);
- the `DB` matrix, that is, the matrix resulting from the multiplication of the
  constitutive matrix (D) and the linear elasticity matrix (B);
- the `Nf` matrix, that is, the boundary integral of the interpolation matrix,
  used for calculating the `f` vector (global and local nodal force vectors);
- the `h` matrix, that is, the elemental tensor used for Helmholtz filtering.

The respective matrices can be calculated by running the script with the
arguments `-k`, `-DB`, and `-Nf`. The result is printed and can be copy-pasted
into the C++ source file, only requiring the necessary variables to be properly
defined in the same scope.

Simpler elements (TRI3 and GT9) have `k` calculated using symbolic integration,
as they are much simpler. Others (like Q4) have much more complex and verbose
integrals, so the resulting matrix is not integrated, being integrated by
Gaussian quadrature during runtime.

As it involves symbolic computations, expect it to be a bit slow, especially
for `k` with high-order elements.

The element `TRI3` doesn't actually use the `k` and `DB` matrices calculated
from `TRI3.py`, as they're trivial to calculate, but it does use `Nf`. The two
matrices have been implemented in the script for reference, but are untested.
