"""
Constructor for Hyperelliptic Curves using the smooth model

Adapted from /hyperelliptic/constructor.py

AUTHORS:

- David Kohel (2006): initial version
- Hyperelliptic-Team (TODO)
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2019 Anna Somoza <anna.somoza.henares@gmail.com>
#                2024 Hyperelliptic-Team (TODO)
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

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.finite_fields import FiniteFields
from sage.rings.abc import pAdicField
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.rational_field import RationalField
from sage.rings.integer import Integer

"""
TODO:

- We currently cannot support the construction of curves over rings

"""


def HyperellipticCurveSmoothModel(f, h=0, check_squarefree=True):
    r"""
    Constructor function for creating a hyperelliptic curve with
    smooth model with polynomials f, h.

    EXAMPLES:

    TODO
    """

    # ---------------------------
    # Internal Helper functions
    # ---------------------------

    def __genus(f, h):
        """
        Helper function to compute the genus of a hyperelliptic curve
        defined by `y^2 + h(x)y = f(x)`.
        """
        # Some classes still have issues with degrees returning `int`
        # rather than Sage Integer types
        df = Integer(f.degree())
        dh_2 = 2 * Integer(h.degree())
        if dh_2 < df:
            return (df - 1) // 2
        return (dh_2 - 1) // 2

    def __check_no_affine_singularities(f, h):
        """
        Helper function which determines whether there are any
        affine singularities in the curve `y^2 + h(x)y = f(x)`.
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

    def __defining_polynomial(f, h):
        """
        Compute the weighted projective model (1 : g + 1 : 1)

        WARNING::

            Due to limitations of the toric_varieties.WP class, we
            cannot instantiate hyperelliptic curves over rings, but
            only fields.
        """
        X, Y, Z = PolynomialRing(f.base_ring(), names="X, Y, Z").gens()

        # Some classes still have issues with degrees returning `int`
        d = max(Integer(h.degree()), (Integer(f.degree()) + 1) // 2)
        F = sum(f[i] * X**i * Z ** (2 * d - i) for i in range(2 * d + 1))

        if h.is_zero():
            G = Y**2 - F
        else:
            H = sum(h[i] * X**i * Z ** (d - i) for i in range(d + 1))
            G = Y**2 + H * Y - F

        return G

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
    defining_polynomial = __defining_polynomial(f, h)

    # -----------------------
    # Class selection
    # -----------------------

    # Special class for finite fields
    if base_ring in FiniteFields():
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_finite_field
        else:
            cls = HyperellipticCurveSmoothModel_finite_field
    # Special class for pAdic fields
    elif isinstance(base_ring, pAdicField):
        if genus == 2:
            cls = HyperellipticCurveSmoothModel_g2_padic_field
        else:
            cls = HyperellipticCurveSmoothModel_padic_field
    # Special class for rational fields
    elif isinstance(base_ring, RationalField):
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

    return cls(defining_polynomial, f, h, genus)
