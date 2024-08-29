r"""
    Constructor for hyperelliptic curves using the smooth model

    In Sage, a hyperelliptic curve of genus `g` is always 
    specified by an (affine) equation in Weierstrass form

    .. MATH::

        y^2 + h(x) y = f(x),

    for some polynomials `h` and `f`. This defines a smooth
    model in weighted projective space `\PP(1 : g + 1 : 1)`

    .. MATH:: 
        Y^2 + H(X,Z) Y = F(X,Z),
    
    where `H` is the degree `g + 1` homogenization of `h`,
    and `F` is the degree `2 g + 2` homogenization of `f`. 

    There are either 0, 1 or 2 points at infinity (`Z=0`),
    in which case we say that the hyperelliptic curve is 
    inert, ramified or split, respectively. 


    EXAMPLES:
    
    We create a hyperelliptic curve over the rationals::

        sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
        sage: R.<x> = PolynomialRing(QQ)
        sage: H = HyperellipticCurveSmoothModel(x^8+1, x^4+1); H
        Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 1)*y = x^8 + 1

    This hyperelliptic curve has no points at infinity,  
    i.e. `H` is inert::
        sage: H.points_at_infinity()
        []
        sage: H.is_inert()
        True

    We can extend the base field to obtain a hyperelliptic curve
    with two points at infinity::
        sage: K.<alpha> = QQ.extension(x^2+x-1)
        sage: HK = H.change_ring(K)
        sage: HK.points_at_infinity()
        [[1 : alpha : 0], [1 : -alpha - 1 : 0]]
        sage: HK.is_split()
        True

    The construction of hyperelliptic curves is supported over different fields.
    The correct class is chosen automatically::
        sage: F = FiniteField(13)
        sage: S.<x> = PolynomialRing(F)
        sage: HF = HyperellipticCurve(x^5 + x^4 + x^3 + x^2 + x + 1, x^3 + x); HF
        Hyperelliptic Curve over Finite Field of size 13 defined by y^2 + (x^3 + x)*y = x^5 + x^4 + x^3 + x^2 + x + 1
        sage: type(HF)
        <class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g2_FiniteField_with_category'>
        sage: Q5 = Qp(5,10)
        sage: T.<x> = Q5[]
        sage: H5 = HyperellipticCurveSmoothModel(x^7 + 1); H5
        Hyperelliptic Curve over 5-adic Field with capped relative precision 10 defined by y^2 = x^7 + 1 + O(5^10)
        sage: type(H5)
        <class 'hyperelliptic_padic_field.HyperellipticCurveSmoothModel_padic_field_with_category'>
    
    The input polynomials need not be monic::
        sage: R.<x> = QQ[]
        sage: HyperellipticCurveSmoothModel(3*x^5+1)
        Hyperelliptic Curve over Rational Field defined by y^2 = 3*x^5 + 1

    The polynomials f and h need to define a smooth curve of genus at 
    least one. In particular polynomials defining elliptic curves are
    allowed as input.
        sage: E = HyperellipticCurveSmoothModel(x^3+1)
        sage: E.genus()
        1
        sage: HyperellipticCurveSmoothModel(x)
        Traceback (most recent call last):
        ...
        ValueError: arguments f = x and h = 0 must define a curve of genus at least one.

    The following polynomials define a singular curve and are 
    not allowed as input::
        sage: C = HyperellipticCurve(x^6 + 2*x - 1, 2*x - 2)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: singularity in the provided affine patch

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

from sage.categories.finite_fields import FiniteFields
from sage.rings.abc import pAdicField
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.rational_field import RationalField
from sage.schemes.toric.library import toric_varieties
from sage.rings.integer import Integer

"""
TODO:

- We currently cannot support the construction of curves over rings

"""


def HyperellipticCurveSmoothModel(f, h=0, check_squarefree=True):
    r"""
    Constructor function for creating a hyperelliptic curve with
    smooth model with polynomials f, h.

    In Sage, a hyperelliptic curve of genus `g` is always 
    specified by an (affine) equation in Weierstrass form

    .. MATH::

        y^2 + h(x) y = f(x),

    for some polynomials `h` and `f`. This defines a smooth
    model in weighted projective space `\PP(1 : g + 1 : 1)`

    .. MATH:: 
        Y^2 + H(X,Z) Y = F(X,Z),
    
    where `H` is the degree `g + 1` homogenization of `h`,
    and `F` is the degree `2 g + 2` homogenization of `f`. 

    INPUT:

    -  ``f`` -- polynomial

    -  ``h`` (default: ``0``) -- polynomial

    -  ``check_squarefree`` (default: ``True``) -- test if
       the input defines a hyperelliptic curve

    EXAMPLES:
    
    We create a hyperelliptic curve over the rationals::

        sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
        sage: R.<x> = PolynomialRing(QQ)
        sage: H = HyperellipticCurveSmoothModel(x^8+1, x^4+1); H
        Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 1)*y = x^8 + 1

    This hyperelliptic curve has no points at infinity,  
    i.e. `H` is inert::
        sage: H.points_at_infinity()
        []
        sage: H.is_inert()
        True

    We can extend the base field to obtain a hyperelliptic curve
    with two points at infinity::
        sage: K.<alpha> = QQ.extension(x^2+x-1)
        sage: HK = H.change_ring(K)
        sage: HK.points_at_infinity()
        [[1 : alpha : 0], [1 : -alpha - 1 : 0]]
        sage: HK.is_split()
        True

    The construction of hyperelliptic curves is supported over different fields.
    The correct class is chosen automatically::
        sage: F = FiniteField(13)
        sage: S.<x> = PolynomialRing(F)
        sage: HF = HyperellipticCurve(x^5 + x^4 + x^3 + x^2 + x + 1, x^3 + x); HF
        Hyperelliptic Curve over Finite Field of size 13 defined by y^2 + (x^3 + x)*y = x^5 + x^4 + x^3 + x^2 + x + 1
        sage: type(HF)
        <class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g2_FiniteField_with_category'>
        sage: Q5 = Qp(5,10)
        sage: T.<x> = Q5[]
        sage: H5 = HyperellipticCurveSmoothModel(x^7 + 1); H5
        Hyperelliptic Curve over 5-adic Field with capped relative precision 10 defined by y^2 = x^7 + 1 + O(5^10)
        sage: type(H5)
        <class 'hyperelliptic_padic_field.HyperellipticCurveSmoothModel_padic_field_with_category'>
    
    The input polynomials need not be monic::
        sage: R.<x> = QQ[]
        sage: HyperellipticCurveSmoothModel(3*x^5+1)
        Hyperelliptic Curve over Rational Field defined by y^2 = 3*x^5 + 1

    The polynomials f and h need to define a smooth curve of genus at 
    least one. In particular polynomials defining elliptic curves are
    allowed as input.
        sage: E = HyperellipticCurveSmoothModel(x^3+1)
        sage: E.genus()
        1
        sage: HyperellipticCurveSmoothModel(x)
        Traceback (most recent call last):
        ...
        ValueError: arguments f = x and h = 0 must define a curve of genus at least one.

    The following polynomials define a singular curve and are 
    not allowed as input::
        sage: C = HyperellipticCurve(x^6 + 2*x - 1, 2*x - 2)
        Traceback (most recent call last):
        ...
        ValueError: not a hyperelliptic curve: singularity in the provided affine patch
    
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
        d = max(Integer(h.degree()), (Integer(f.degree()) + 1) // 2)
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
    if genus == 0:
        raise ValueError(f"arguments {f = } and {h = } must define a curve of genus at least one.")

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

    return cls(projective_model, f, h, genus)
