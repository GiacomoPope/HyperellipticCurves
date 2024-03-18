"""
Constructor for Hyperelliptic Curves using the smooth model

Adapted from /hyperelliptic/constructor.py

AUTHORS:

- David Kohel (2006): initial version
- TODO
"""
# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2019 Anna Somoza <anna.somoza.henares@gmail.com>
#                2024 TODO
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from hyperelliptic_generic import HyperellipticCurveSmoothModel_generic
from hyperelliptic_g2 import (
    HyperellipticCurveSmoothModel_g2,
    HyperellipticCurveSmoothModel_g2_finite_field,
    HyperellipticCurveSmoothModel_g2_rational_field,
    HyperellipticCurveSmoothModel_g2_padic_field,
)
from hyperelliptic_finite_field import HyperellipticCurveSmoothModel_finite_field
from hyperelliptic_rational_field import HyperellipticCurveSmoothModel_rational_field
from hyperelliptic_padic_field import HyperellipticCurveSmoothModel_padic_field

from sage.categories.finite_fields import FiniteFields
from sage.rings.abc import pAdicField
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.rational_field import is_RationalField
from sage.schemes.toric.library import toric_varieties
from sage.rings.integer import Integer

"""
TODO:

- We currently cannot support the construction of curves over rings
- ...

"""


def HyperellipticCurveSmoothModel(f, h=0, check_squarefree=True):
    r"""
    TODO
    """

    # -----------------
    # Helper functions
    # -----------------

    def __genus(f, h):
        """
        TODO
        """
        # Some classes still have issues with degrees returning `int`
        df = Integer(f.degree())
        dh_2 = 2 * Integer(h.degree())
        if dh_2 < df:
            return (df - 1) // 2
        return (dh_2 - 1) // 2

    def __check_no_affine_singularities(f, h):
        """
        TODO
        """
        if f.base_ring().characteristic() == 2:
            if h.is_zero():
                return False
            elif h.is_constant():
                return True
            return h.gcd(f.derivative() ** 2 - f * h.derivative() ** 2).is_one()

        if h.is_zero():
            return f.gcd(f.derivative()).is_one()

        g1 = h**2 + 4 * f
        g2 = 2 * f.derivative() + h * h.derivative()
        return g1.gcd(g2).is_one()

    def __projective_model(f, h, genus):
        """
        Compute the weighted projective model (1 : g + 1 : 1)
        """
        T = toric_varieties.WP(
            [1, genus + 1, 1], base_ring=f.base_ring(), names="X, Y, Z"
        )
        (X, Y, Z) = T.gens()

        # Some classes still have issues with degrees returning `int`
        d = max(Integer(h.degree()), (Integer(f.degree()) / 2).ceil())
        F = sum(f[i] * X**i * Z ** (2 * d - i) for i in range(2 * d + 1))

        if h.is_zero():
            G = Y**2 - F
        else:
            H = sum(h[i] * X**i * Z ** (d - i) for i in range(d + 1))
            G = Y**2 + H * Y - F

        return T.subscheme(G)

    # -------------------------------------------
    # Typechecking and projective model creation
    # -------------------------------------------

    # Check the polynomials are of the right type
    F = h**2 + 4 * f
    if not isinstance(F, Polynomial):
        raise TypeError(f"arguments {f = } and {h = } must be polynomials")

    # Store the hyperelliptic polynomials as the correct type
    polynomial_ring = F.parent()
    base_ring = F.base_ring()
    f = polynomial_ring(f)
    h = polynomial_ring(h)

    # Ensure that there are no affine singular points
    if check_squarefree and not __check_no_affine_singularities(f, h):
        raise ValueError("singularity in the provided affine patch")

    # Compute the genus of the curve from f, h
    genus = __genus(f, h)

    # Compute the smooth model for the hyperelliptic curve
    # using a weighted projective space (via Toric Variety)
    projective_model = __projective_model(f, h, genus)

    # -----------------------
    # Class selection
    # -----------------------

    # Special class for finite fields
    if base_ring in FiniteFields():
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_finite_field
        else:
            cls = HyperellipticCurveSmoothModel_finite_field
    # Special class for padic fields
    elif isinstance(base_ring, pAdicField):
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_padic_field
        else:
            cls = HyperellipticCurveSmoothModel_padic_field
    # Special class for rational fields
    elif is_RationalField(base_ring):
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_rational_field
        else:
            cls = HyperellipticCurveSmoothModel_rational_field
    # Default class for all other fields
    else:
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2
        else:
            cls = HyperellipticCurveSmoothModel_generic

    return cls(projective_model, f, h, genus)
