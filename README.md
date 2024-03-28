# Hyperelliptic Curves

The current implementation of hyperelliptic curves in SageMath uses the
projective plane model. Although this works nicely enough for imaginary curves
with only one point at infinity, it is not descriptive enough for the real
models.

This repository is a total rewrite of the hyperelliptic curve classes to instead
use the smooth model for hyperelliptic curves which facilitates implementing
arithmetic of Jacobians of hyperelliptic curves for (almost) all cases.

The aim is ultimately to introduce the new class `HyperellipticCurveSmoothModel`
into SageMath as an alternative class to `HyperellipticCurve` with the potential
to deprecate and replace this model in the future.

## Usage

`HyperellipticCurveSmoothModel` has been designed as a drop in replacement for `HyperellipticCurve`, so from the
root of this repo run the following:

```py
from hyperelliptic_constructor import HyperellipticCurveSmoothModel
R.<x> = QQ[]
f = x^7 + 1
H = HyperellipticCurveSmoothModel(f)
```

and you should be able to work "as normal" with `H`, but without arithmetic bugs

## Progress

First we should copy everything from sage and make sure it
stills works:

### Hyperelliptic Curves

- [x] `Hyperelliptic_generic` ~~(need to implement `rational_points()`)~~
- [x] `Hyperelliptic_g2`
- [x] `Hyperelliptic_rational_field`
- [x] `Hyperelliptic_finite_field`
- [x] `Hyperelliptic_padic_field`

### Jacobians

- [ ] `jacobian_generic`
- [ ] `jacobian_g2`
- [ ] `jacobian_homset`
- [ ] `jacobian_morphism`

### Other

- [x] `mestre.py`
- [ ] `monsky_washnitzer.py` (failing examples as we don't support curves over rings)
- [ ] `kummer_surface.py` (to rewrite)
- [x] `hyperellfrob` library (no changes needed)

## Plan

To maintain functionality there are several files which need to be rewritten to
work with the smooth model.

### Hyperelliptic Curves

- Write a constructor function which picks an appropriate
  `HyperellipticSmothModel_*` class determined by the base ring of the defining
  polynomials
- Rewrite the `Hyperelliptic_generic` class as
  `HyperellipticSmoothModel_generic`
- Rewrite the `Hyperelliptic_rational_field` class as
  `HyperellipticSmoothModel_rational_field`
- Rewrite the `Hyperelliptic_finite_field` class as
  `HyperellipticSmoothModel_finite_field`
- Rewrite the `Hyperelliptic_padic_field` class as
  `HyperellipticSmoothModel_padic_field`

Additionally, we should maintain the fast computation of the matrix of frobenius
from `hypellrob.pyx` and also copy over the work from `mestre.py` which
computes:

- hyperelliptic curves from the Igusa-Clebsch invariants (over `\QQ` and finite
  fields)
- Mestre's conic from the Igusa-Clebsch invariants

The case of genus two should also be singled out as there are additional methods
here. For example, computing certain invariants of the derived Kummer surface.

as well as `monsky_washnitzer.py` which performs the computation of Frobenius
matrix on Monsky-Washnitzer cohomology

### Jacobians of Hyperelliptic Curves

- Rewrite `jacobian_generic.py` for the Jacobian of the curve. This needs only
  the most basic features and should be essentially unaffected by the change of
  model. I could however be wrong here though, maybe we will need child classes
  based on the number of points at infinty.
- Rewrite `jacobian_homset.py` which is the set of rational points on the
  Jacobian. This is where the creation of divisors and arithmetic is performed.
  This will need careful handling for the ramified, split and inert models. We
  need to be extra careful for the split/inert case as a change in the base ring
  may change the number of points at infinity.
- Rewrite `jacobian_morphism.py` which represents divisors on the jacobian. Here
  we need to again have a base class and then have children classes based on the
  number of points at infinity as the divisors themselves have different
  representations.
