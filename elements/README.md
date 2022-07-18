# Element creation scripts

Python scripts designed to aid the creation of element classes. They depend
on `numpy` and `sympy`.

The scripts are used mainly to calculate the element-dependent matrices:
- the stiffness matrix (`k`);
- the `DB` matrix, that is, the matrix resulting from the multiplication of the
  constitutive matrix (D) and the linear elasticity matrix (B);
- the `Nf` matrix, that is, the boundary integral of the interpolation matrix,
  used for calculating the `f` vector (global and local nodal force vectors).

The respective matrices can be calculated by running the script with the
arguments `-k`, `-DB`, and `-Nf`. The result is printed and can be copy-pasted
into the C++ source file, only requiring the necessary variables to be properly
defined in the same scope.

Unlike what's usually done, the matrices are calculated using symbolic methods.
The resulting formulas are therefore as exact as can be, but the results
printed tend to be quite verbose because of that.

It was the quickest method to obtain those matrices at the time I implemented
them, being based on the method used by [topy](https://github.com/williamhunter/topy).
It's unknown if it's faster than the usual method of Gaussian numeric
integration, but it's definitely more exact. It may or may not benefit from
the use of `-ffast-math`, but this is untested.

To extend it, just take one of the elements as a base and work on it from
there. You'll mostly need to change the interpolation functions and how they
work with derivation and integration.

As it involves symbolic computations, expect it to be a bit slow, especially
for `k` with high-order elements.

The element `TRI3` doesn't actually use the `k` and `DB` matrices calculated
from `TRI3.py`, as they're trivial to calculate, but it does use `Nf`. The two
matrices have been implemented in the script for reference, but are untested.
