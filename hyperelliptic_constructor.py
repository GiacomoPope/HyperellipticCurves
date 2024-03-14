"""
Constructor for Hyperelliptic Curves using the smooth model

Adapted from /hyperelliptic/constructor.py 

AUTHORS:

- David Kohel (2006): initial version
- Anna Somoza (2019-04): dynamic class creation
- TODO
"""
# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2019 Anna Somoza <anna.somoza.henares@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from hyperelliptic_generic import HyperellipticCurveSmoothModel_generic
from hyperelliptic_g2 import HyperellipticCurveSmoothModel_g2
from hyperelliptic_finite_field import HyperellipticCurveSmoothModel_finite_field
from hyperelliptic_rational_field import HyperellipticCurveSmoothModel_rational_field
from hyperelliptic_padic_field import HyperellipticCurveSmoothModel_padic_field

# TODO: rewrite this class for genus two specifically
# from .hyperelliptic_g2 import HyperellipticCurve_g2

# TODO: use categories to determine which field the base ring is?
from sage.categories.finite_fields import FiniteFields
from sage.rings.abc import pAdicField
from sage.rings.rational_field import is_RationalField

from sage.structure.dynamic_class import dynamic_class

from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.toric.library import toric_varieties


def HyperellipticCurveSmoothModel(f, h=0):
    r"""
    TODO
    """

    def __genus(f, h):
        """
        TODO
        """
        df = f.degree()
        dh_2 = 2 * h.degree()
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

        d = max(h.degree(), (f.degree() / 2).ceil())
        F = sum(f[i] * X**i * Z ** (2 * d - i) for i in range(2 * d + 1))

        if h.is_zero():
            G = Y**2 - F
        else:
            H = sum(h[i] * X**i * Z ** (d - i) for i in range(d + 1))
            G = Y**2 + H * Y - F

        return T.subscheme(G)

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
    if not __check_no_affine_singularities(f, h):
        raise ValueError("singularity in the provided affine patch")

    # Compute the genus of the curve from f, h
    genus = __genus(f, h)

    # Compute the smooth model for the hyperelliptic curve
    # using a weighted projective space (via Toric Variety)
    projective_model = __projective_model(f, h, genus)

    # ----------------------
    # Dynamic class creation
    # ----------------------

    superclass = []
    cls_name = ["HyperellipticCurveSmoothModel"]
    genus_classes = {2: HyperellipticCurveSmoothModel_g2}

    # Look to see if the curve genus puts us in a special case
    # currently only genus two has additional methods
    if genus in genus_classes:
        superclass.append(genus_classes[genus])
        cls_name.append(f"g{genus}")

    def is_FiniteField(x):
        return x in FiniteFields()

    def is_pAdicField(x):
        # TODO: is there a category way to do this?
        return isinstance(x, pAdicField)

    fields = [
        ("FiniteField", is_FiniteField, HyperellipticCurveSmoothModel_finite_field),
        (
            "RationalField",
            is_RationalField,
            HyperellipticCurveSmoothModel_rational_field,
        ),
        ("pAdicField", is_pAdicField, HyperellipticCurveSmoothModel_padic_field),
    ]

    # TODO:
    # Is this really the best way to do this construction???
    # Test the base field to check if it is in the singled out
    # fields with additional methods
    for name, test, cls in fields:
        if test(base_ring):
            superclass.append(cls)
            cls_name.append(name)
            break

    base_cls = HyperellipticCurveSmoothModel_generic
    if len(superclass) != 0:
        base_cls = None

    class_name = "_".join(cls_name)
    cls = dynamic_class(
        class_name,
        tuple(superclass),
        cls=base_cls,
        doccls=HyperellipticCurveSmoothModel,
    )
    return cls(projective_model, f, h, genus)
