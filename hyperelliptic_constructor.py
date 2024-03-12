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
from hyperelliptic_finite_field import HyperellipticCurveSmoothModel_finite_field
from hyperelliptic_rational_field import HyperellipticCurveSmoothModel_rational_field
from hyperelliptic_padic_field import HyperellipticCurveSmoothModel_padic_field

# TODO: rewrite this class for genus two specifically
# from .hyperelliptic_g2 import HyperellipticCurve_g2

# TODO: use categories to determine which field the base ring is
import sage.rings.abc
from sage.rings.finite_rings.finite_field_base import FiniteField
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
    
        return T.subscheme([G])

    # Check the polynomials are of the right type
    F = h**2 + 4 * f
    if not isinstance(F, Polynomial):
        raise TypeError(f"arguments {f = } and {h = } must be polynomials")

    # Ensure that there are no affine singular points
    if not __check_no_affine_singularities(f, h):
        raise ValueError("singularity in the provided affine patch")
    
    # Store the hyperelliptic polynomials as the correct type
    polynomial_ring = F.parent()
    f = polynomial_ring(f)
    h = polynomial_ring(h)

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

    # TODO:
    # This means we must compute the genus on construction
    # We could then take this and pass it as a parameter to the specific classes
    # TODO genus_classes = {2: HyperellipticCurve_g2}

    def is_FiniteField(x):
        return isinstance(x, FiniteField)

    def is_pAdicField(x):
        return isinstance(x, sage.rings.abc.pAdicField)

    fields = [
        ("FiniteField", is_FiniteField, HyperellipticCurveSmoothModel_finite_field),
        ("RationalField", is_RationalField, HyperellipticCurveSmoothModel_rational_field),
        ("pAdicField", is_pAdicField, HyperellipticCurveSmoothModel_padic_field)]

    # TODO
    # if g in genus_classes:
    #     superclass.append(genus_classes[g])
    #     cls_name.append("g%s" % g)

    # Computing F from f, h does the coercion for us
    F = h**2 + 4*f
    k = F.base_ring()

    # TODO:
    # Is this really the best way to do this construction???
    for name, test, cls in fields:
        if test(k):
            superclass.append(cls)
            cls_name.append(name)
            break

    class_name = "_".join(cls_name)
    cls = dynamic_class(class_name, tuple(superclass),
                        HyperellipticCurveSmoothModel_generic, doccls=HyperellipticCurveSmoothModel)
    return cls(projective_model, f, h, genus)

