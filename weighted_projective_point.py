<<<<<<< HEAD
r"""
Points on projective schemes

This module implements scheme morphism for points on weighted projective schemes.

AUTHORS:

- Gareth Ma (2024): initial version, based on unweighted version.
"""

# ****************************************************************************
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.integral_domains import IntegralDomains
from sage.categories.number_fields import NumberFields
from sage.rings.integer_ring import ZZ
from sage.rings.fraction_field import FractionField
from sage.rings.number_field.order import Order as NumberFieldOrder
from sage.rings.qqbar import number_field_elements_from_algebraics
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.rational_field import QQ
from sage.arith.misc import GCD as gcd
from sage.arith.functions import lcm
from sage.misc.misc_c import prod

from copy import copy
from sage.schemes.generic.morphism import (SchemeMorphism,
                                           is_SchemeMorphism,
                                           SchemeMorphism_point)
from sage.structure.element import AdditiveGroupElement
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp, op_EQ, op_NE


# --------------------
# Projective varieties
# --------------------

class SchemeMorphism_point_weighted_projective_ring(SchemeMorphism_point):
    """
    A rational point of projective space over a ring.

    INPUT:

    -  ``X`` -- a homset of a subscheme of an ambient projective space over a ring `K`.

    - ``v`` -- a list or tuple of coordinates in `K`.

    - ``check`` -- boolean (default:``True``). Whether to check the input for consistency.

    EXAMPLES::

        sage: P = ProjectiveSpace(2, ZZ)
        sage: P(2,3,4)
        (2 : 3 : 4)

    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        EXAMPLES::

            TODO
        """
        SchemeMorphism.__init__(self, X)

        if check:
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            from sage.rings.ring import CommutativeRing
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            else:
                try:
                    if isinstance(v.parent(), CommutativeRing):
                        v = [v]
                except AttributeError:
                    pass
            if not isinstance(v, (list, tuple)):
                raise TypeError("argument v (= %s) must be a scheme point, list, or tuple" % str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components" % (v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(R.one())

            if R in IntegralDomains():
                # Over integral domains, any tuple with at least one
                # non-zero coordinate is a valid projective point.
                if not any(v):
                    raise ValueError(f"{v} does not define a valid projective "
                                     "point since all entries are zero")
            else:
                # Over rings with zero divisors, a more careful check
                # is required: We test whether the coordinates of the
                # point generate the unit ideal. See #31576.
                if 1 not in R.ideal(v):
                    raise ValueError(f"{v} does not define a valid projective point "
                                     "since it is a multiple of a zero divisor")

            X.extended_codomain()._check_satisfies_equations(v)

        self._coords = tuple(v)
        self._normalized = False

    # def _repr_(self):
    #     return "({})".format(" : ".join(map(repr, self._coords)))

    def _richcmp_(self, right, op):
        """
        Test the projective equality of two points.

        INPUT:

        - ``right`` -- a point on projective space

        OUTPUT:

        Boolean

        EXAMPLES::

            TODO
        """
        if not isinstance(right, SchemeMorphism_point):
            try:
                right = self.codomain()(right)
            except TypeError:
                return NotImplemented
        if self.codomain() != right.codomain():
            return op == op_NE

        n = len(self._coords)
        if op in [op_EQ, op_NE]:
            b = all(self[i] * right[j] == self[j] * right[i]
                    for i in range(n) for j in range(i + 1, n))
            return b == (op == op_EQ)
        return richcmp(self._coords, right._coords, op)

    def __hash__(self):
        """
        Compute the hash value of this point.

        If the base ring has a fraction field, normalize the point in
        the fraction field and then hash so that equal points have
        equal hash values. If the base ring is not an integral domain,
        return the hash of the parent.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: hash(P([1, 1])) == hash(P.point([2, 2], False))
            True

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^2 + 3)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: hash(P([1 + w, 2])) == hash(P([2, 1 - w]))
            True

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(Zmod(10), 1)
            sage: Q.<x,y> = ProjectiveSpace(Zmod(10), 1)
            sage: hash(P([2, 5])) == hash(Q([2, 5]))
            True
            sage: hash(P([2, 5])) == hash(P([2, 5]))
            True
            sage: hash(P([3, 7])) == hash(P([2, 5]))
            True
        """
        R = self.codomain().base_ring()
        #if there is a fraction field normalize the point so that
        #equal points have equal hash values
        if R in IntegralDomains():
            P = self.change_ring(FractionField(R))
            P.normalize_coordinates()
            return hash(tuple(P))
        #if there is no good way to normalize return
        #a constant value
        return hash(self.codomain())

    def _matrix_times_point_(self, mat, dom):
        r"""
        Multiplies the point by a matrix ``mat`` on the left.

        INPUT:

        - ``mat`` -- a matrix

        - ``dom`` -- point set of the resulting point

        OUTPUT: a scheme point given by ``mat*self``

        EXAMPLES::

            sage: P = ProjectiveSpace(QQ, 1)
            sage: Q = P(1,1)
            sage: m = matrix(QQ, 2, 2, [1,1, 0,1])                                      # needs sage.modules
            sage: m*Q                                                                   # needs sage.modules
            (2 : 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x - y)
            sage: Q = X(1,1)
            sage: m = matrix(CC, 3, 3, [1,CC.0,0, CC.0,1,0, 1,1,1])                     # needs sage.modules
            sage: m*Q                                                                   # needs sage.modules
            (0.333333333333333 + 0.333333333333333*I : 0.333333333333333
            + 0.333333333333333*I : 1.00000000000000)

        ::

            sage: # needs sage.modules sage.rings.number_field sage.symbolic
            sage: P = ProjectiveSpace(QQbar, 1)
            sage: Q = P(QQbar(sqrt(2)),1)
            sage: m = matrix(ZZ, 2, 2, [1,-1, 0,1])
            sage: m*Q
            (0.4142135623730951? : 1)

        ::

            sage: P = ProjectiveSpace(QQ, 1)
            sage: Q = P(1,1)
            sage: m = matrix(QQ, 3, 2, [1,1, 0,1, 1,1])                                 # needs sage.modules
            sage: m*Q                                                                   # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square
        """
        from sage.modules.free_module_element import vector
        if not mat.is_square():
            raise ValueError("matrix must be square")
        if mat.ncols() != self.codomain().ngens():
            raise ValueError("matrix size is incompatible")
        X = mat * vector(list(self))
        return dom.codomain()(list(X))

    def scale_by(self, t):
        """
        Scale the coordinates of the point by ``t``.

        A :class:`TypeError` occurs if the point is not in the
        base_ring of the codomain after scaling.

        INPUT:

        - ``t`` -- a ring element.

        OUTPUT: None.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P = ProjectiveSpace(R, 2, 'x')
            sage: p = P([3/5*t^3, 6*t, t])
            sage: p.scale_by(1/t); p
            (3/5*t^2 : 6 : 1)

        ::

            sage: # needs sage.libs.pari
            sage: R.<t> = PolynomialRing(QQ)
            sage: S = R.quo(R.ideal(t^3))
            sage: P.<x,y,z> = ProjectiveSpace(S, 2)
            sage: Q = P(t, 1, 1)
            sage: Q.scale_by(t);Q
            (tbar^2 : tbar : tbar)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: Q = P(2, 2, 2)
            sage: Q.scale_by(1/2);Q
            (1 : 1 : 1)
        """
        if t.is_zero():  #what if R(t) == 0 ?
            raise ValueError("Cannot scale by 0")
        R = self.codomain().base_ring()
        if isinstance(R, QuotientRing_generic):
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                new_coords = [R(u.lift()*t) for u in self._coords]
        else:
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                new_coords = [R(u*t) for u in self._coords]
        self._coords = tuple(new_coords)
        self._normalized = False

    def normalize_coordinates(self):
        """
        Removes the gcd from the coordinates of this point (including `-1`).

        .. WARNING:: The gcd will depend on the base ring.

        OUTPUT: None.

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 2, 'x')
            sage: p = P([-5, -15, -20])
            sage: p.normalize_coordinates(); p
            (1 : 3 : 4)

        ::

            sage: # needs sage.rings.padics
            sage: P = ProjectiveSpace(Zp(7), 2, 'x')
            sage: p = P([-5, -15, -2])
            sage: p.normalize_coordinates(); p
            (5 + O(7^20) : 1 + 2*7 + O(7^20) : 2 + O(7^20))

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P = ProjectiveSpace(R, 2, 'x')
            sage: p = P([3/5*t^3, 6*t, t])
            sage: p.normalize_coordinates(); p
            (3/5*t^2 : 6 : 1)

        ::

            sage: P.<x,y> = ProjectiveSpace(Zmod(20), 1)
            sage: Q = P(3, 6)
            sage: Q.normalize_coordinates()
            sage: Q
            (1 : 2)

        Since the base ring is a polynomial ring over a field, only the
        gcd `c` is removed. ::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P = ProjectiveSpace(R, 1)
            sage: Q = P(2*c, 4*c)
            sage: Q.normalize_coordinates();Q
            (2 : 4)

        A polynomial ring over a ring gives the more intuitive result. ::

            sage: R.<c> = PolynomialRing(ZZ)
            sage: P = ProjectiveSpace(R, 1)
            sage: Q = P(2*c, 4*c)
            sage: Q.normalize_coordinates();Q
            (1 : 2)

        ::

            sage: # needs sage.libs.singular
            sage: R.<t> = QQ[]
            sage: S = R.quotient_ring(R.ideal(t^3))
            sage: P.<x,y> = ProjectiveSpace(S, 1)
            sage: Q = P(t + 1, t^2 + t)
            sage: Q.normalize_coordinates()
            sage: Q
            (1 : tbar)
        """
        if self._normalized:
            return
        R = self.codomain().base_ring()
        if isinstance(R, QuotientRing_generic):
            index = len(self._coords) - 1
            while not self._coords[index]:
                index -= 1
            last = self._coords[index].lift()
            mod, = R.defining_ideal().gens()
            unit = last
            while not (zdiv := mod.gcd(unit)).is_unit():
                unit //= zdiv
            self.scale_by(unit.inverse_mod(mod))
        else:
            GCD = R(gcd(self._coords[0], self._coords[1]))
            index = 2
            neg = self._coords[0] <= 0 and self._coords[1] <= 0
            while not GCD.is_one() and index < len(self._coords):
                neg = self._coords[index] <= 0
                GCD = R(gcd(GCD, self._coords[index]))
                index += 1
            if not GCD.is_one():
                self.scale_by(~GCD)
            if neg:
                self.scale_by(-ZZ.one())
        self._normalized = True

    def dehomogenize(self,n):
        r"""
        Dehomogenizes at the `n`-th coordinate.

        INPUT:

        - ``n`` -- non-negative integer.

        OUTPUT:

        - :class:`SchemeMorphism_point_affine`.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: Q = X(23, 23, 46)
            sage: Q.dehomogenize(2)                                                     # needs sage.libs.singular
            (1/2, 1/2)

        ::

            sage: # needs sage.libs.pari
            sage: R.<t> = PolynomialRing(QQ)
            sage: S = R.quo(R.ideal(t^3))
            sage: P.<x,y,z> = ProjectiveSpace(S, 2)
            sage: Q = P(t, 1, 1)
            sage: Q.dehomogenize(1)
            (tbar, 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: Q = P(1, 3, 1)
            sage: Q.dehomogenize(0)
            (3, 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: Q = P(1, 3, 0)
            sage: Q.dehomogenize(2)
            Traceback (most recent call last):
            ...
            ValueError: can...t dehomogenize at 0 coordinate
        """
        if self[n].is_zero():
            raise ValueError("can't dehomogenize at 0 coordinate")
        PS = self.codomain()
        A = PS.affine_patch(n)
        Q = []
        for i in range(PS.ambient_space().dimension_relative() + 1):
            if i != n:
                Q.append(self[i] / self[n])
        return A.point(Q)

class SchemeMorphism_point_weighted_projective_field(SchemeMorphism_point_weighted_projective_ring):
    """
    A rational point of projective space over a field.

    INPUT:

    -  ``X`` -- a homset of a subscheme of an ambient projective space
       over a field `K`.

    - ``v`` -- a list or tuple of coordinates in `K`.

    - ``check`` -- boolean (default:``True``). Whether to
      check the input for consistency.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr
        sage: P = ProjectiveSpace(3, RR)
        sage: P(2, 3, 4, 5)
        (0.400000000000000 : 0.600000000000000 : 0.800000000000000 : 1.00000000000000)
    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_point_projective_ring` for details.

        This function still normalizes points so that the rightmost non-zero coordinate is 1.
        This is to maintain functionality with current
        implementations of curves in projectives space (plane, conic, elliptic, etc).
        The :class:`SchemeMorphism_point_projective_ring` is for general use.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, QQ)
            sage: P(2, 3/5, 4)
            (1/2 : 3/20 : 1)

        ::

            sage: P = ProjectiveSpace(3, QQ)
            sage: P(0, 0, 0, 0)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0, 0, 0] does not define a valid projective point since all entries are zero

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X = P.subscheme([x^2 - y*z])
            sage: X([2, 2, 2])
            (1 : 1 : 1)

        ::

            sage: P = ProjectiveSpace(1, GF(7))
            sage: Q = P([2, 1])
            sage: Q[0].parent()
            Finite Field of size 7

        ::

            sage: P = ProjectiveSpace(QQ, 1)
            sage: P.point(Infinity)
            (1 : 0)
            sage: P(infinity)
            (1 : 0)

        ::

            sage: P = ProjectiveSpace(QQ, 2)
            sage: P(infinity)
            Traceback (most recent call last):
            ...
            ValueError: +Infinity not well defined in dimension > 1
            sage: P.point(infinity)
            Traceback (most recent call last):
            ...
            ValueError: +Infinity not well defined in dimension > 1
        """
        SchemeMorphism.__init__(self, X)

        self._normalized = False

        if check:
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            from sage.rings.ring import CommutativeRing
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            else:
                try:
                    if isinstance(v.parent(), CommutativeRing):
                        v = [v]
                except AttributeError:
                    pass
            if not isinstance(v, (list,tuple)):
                raise TypeError("argument v (= %s) must be a scheme point, list, or tuple" % str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components" % (v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(R.one())

            for last in reversed(range(len(v))):
                c = v[last]
                if c.is_one():
                    break
                if c:
                    for j in range(last):
                        v[j] /= c
                    v[last] = R.one()
                    break
            else:
                raise ValueError(f"{v} does not define a valid projective "
                                 "point since all entries are zero")
            self._normalized = True

            X.extended_codomain()._check_satisfies_equations(v)

        self._coords = tuple(v)

    def __hash__(self):
        """
        Computes the hash value of this point.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: hash(P([1/2, 1])) == hash(P.point([1, 2], False))
            True
        """
        P = copy(self)
        P.normalize_coordinates()
        return hash(tuple(P))

    def normalize_coordinates(self):
        r"""
        Normalizes the point so that the last non-zero coordinate is `1`.

        OUTPUT: None.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: Q = P.point([GF(5)(1), GF(5)(3), GF(5)(0)], False); Q
            (1 : 3 : 0)
            sage: Q.normalize_coordinates(); Q
            (2 : 1 : 0)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x^2 - y^2);
            sage: Q = X.point([23, 23, 46], False); Q
            (23 : 23 : 46)
            sage: Q.normalize_coordinates(); Q
            (1/2 : 1/2 : 1)
        """
        if self._normalized:
            return
        for index in reversed(range(len(self._coords))):
            c = self._coords[index]
            if c.is_one():
                break
            if c:
                inv = c.inverse()
                new_coords = [d * inv for d in self._coords[:index]]
                new_coords.append(self.base_ring().one())
                new_coords.extend(self._coords[index+1:])
                self._coords = tuple(new_coords)
                break
        else:
            assert False, 'bug: invalid projective point'
        self._normalized = True

    def _number_field_from_algebraics(self):
        r"""
        Given a projective point defined over ``QQbar``, return the same point, but defined
        over a number field.

        This is only implemented for points of projective space.

        OUTPUT: scheme point

        EXAMPLES::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<x> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: Q = P([-1/2*QQbar(sqrt(2)) + QQbar(I), 1])
            sage: S = Q._number_field_from_algebraics(); S
            (1/2*a^3 + a^2 - 1/2*a : 1)
            sage: S.codomain()
            Projective Space of dimension 1 over Number Field in a with defining
             polynomial y^4 + 1 with a = 0.7071067811865475? + 0.7071067811865475?*I

        The following was fixed in :issue:`23808`::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<x> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: Q = P([-1/2*QQbar(sqrt(2)) + QQbar(I), 1]);Q
            (-0.7071067811865475? + 1*I : 1)
            sage: S = Q._number_field_from_algebraics(); S
            (1/2*a^3 + a^2 - 1/2*a : 1)
            sage: T = S.change_ring(QQbar)  # Used to fail
            sage: T
            (-0.7071067811865475? + 1.000000000000000?*I : 1)
            sage: Q[0] == T[0]
            True
        """
        from sage.schemes.projective.projective_space import ProjectiveSpace_ring
        if not isinstance(self.codomain(), ProjectiveSpace_ring):
            raise NotImplementedError("not implemented for subschemes")

        # Issue #23808: Keep the embedding info associated with the number field K
        # used below, instead of in the separate embedding map phi which is
        # forgotten.
        K_pre,P,phi = number_field_elements_from_algebraics(list(self))
        if K_pre is QQ:
            K = QQ
        else:
            from sage.rings.number_field.number_field import NumberField
            K = NumberField(K_pre.polynomial(), embedding=phi(K_pre.gen()), name='a')
            psi = K_pre.hom([K.gen()], K) # Identification of K_pre with K
            P = [ psi(p) for p in P ] # The elements of P were elements of K_pre
        from sage.schemes.projective.projective_space import ProjectiveSpace
        PS = ProjectiveSpace(K,self.codomain().dimension_relative(),'z')
        return PS(P)

    def clear_denominators(self):
        r"""
        scales by the least common multiple of the denominators.

        OUTPUT: None.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: Q = P([t, 3/t^2, 1])
            sage: Q.clear_denominators(); Q
            (t^3 : 3 : t^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: Q = P([1/w, 3, 0])
            sage: Q.clear_denominators(); Q
            (w : 9 : 0)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: Q = X([1/2, 1/2, 1])
            sage: Q.clear_denominators(); Q
            (1 : 1 : 2)

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: Q = PS.point([1, 2/3], False); Q
            (1 : 2/3)
            sage: Q.clear_denominators(); Q
            (3 : 2)
        """
        self.scale_by(lcm([t.denominator() for t in self]))

    def intersection_multiplicity(self, X):
        r"""
        Return the intersection multiplicity of the codomain of this point and ``X`` at this point.

        This uses the intersection_multiplicity implementations for projective/affine subschemes. This
        point must be a point of a projective subscheme.

        INPUT:

        - ``X`` -- a subscheme in the same ambient space as that of the codomain of this point.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x*z - y^2])
            sage: Y = P.subscheme([x^3 - y*w^2 + z*w^2, x*y - z*w])
            sage: Q1 = X([1/2, 1/4, 1/8, 1])
            sage: Q1.intersection_multiplicity(Y)                                       # needs sage.libs.singular
            1
            sage: Q2 = X([0,0,0,1])
            sage: Q2.intersection_multiplicity(Y)                                       # needs sage.libs.singular
            5
            sage: Q3 = X([0,0,1,0])
            sage: Q3.intersection_multiplicity(Y)                                       # needs sage.libs.singular
            6

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x^2 - y^2])
            sage: Q = P([1,1,1,0])
            sage: Q.intersection_multiplicity(X)
            Traceback (most recent call last):
            ...
            TypeError: this point must be a point on a projective subscheme
        """
        from sage.schemes.projective.projective_space import ProjectiveSpace_ring
        if isinstance(self.codomain(), ProjectiveSpace_ring):
            raise TypeError("this point must be a point on a projective subscheme")
        return self.codomain().intersection_multiplicity(X, self)

    def multiplicity(self):
        r"""
        Return the multiplicity of this point on its codomain.

        Uses the subscheme multiplicity implementation. This point must be a point on
        a projective subscheme.

        OUTPUT: an integer.

        EXAMPLES::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(QQ, 4)
            sage: X = P.subscheme([y^6 - x^3*w^2*t + t^5*w, x^2 - t^2])
            sage: Q1 = X([1,0,2,1,1])
            sage: Q1.multiplicity()                                                     # needs sage.libs.singular
            1
            sage: Q2 = X([0,0,-2,1,0])
            sage: Q2.multiplicity()                                                     # needs sage.libs.singular
            8
        """
        from sage.schemes.projective.projective_space import ProjectiveSpace_ring
        if isinstance(self.codomain(), ProjectiveSpace_ring):
            raise TypeError("this point must be a point on a projective subscheme")
        return self.codomain().multiplicity(self)

    def as_subscheme(self):
        r"""
        Return the subscheme associated with this rational point.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: p1 = P2.point([0,0,1]).as_subscheme(); p1
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x, y
            sage: p2 = P2.point([1,1,1]).as_subscheme(); p2
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x - z, y - z
            sage: p1 + p2
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x - y, y^2 - y*z
        """
        P = self.codomain().ambient_space()
        g = P.gens()
        v = self._coords
        n = len(v)
        for i in range(n - 1, -1, -1):
            if v[i]:
                break
        a = v[i]
        x = g[i]
        return P.subscheme([a*g[j] - v[j]*x for j in range(n) if j != i])


class SchemeMorphism_point_projective_finite_field(SchemeMorphism_point_projective_field):

    def __hash__(self):
        r"""
        Returns the integer hash of this point.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: hash(P(2, 1, 2))
            41

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: hash(X(1, 1, 2))
            81

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13), 1)
            sage: hash(P(3, 4))
            17

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13^3,'t'), 1)                            # needs sage.rings.finite_rings
            sage: hash(P(3, 4))                                                         # needs sage.rings.finite_rings
            2201
        """
        p = self.codomain().base_ring().order()
        N = self.codomain().ambient_space().dimension_relative()
        return hash(sum(hash(self[i]) * p**i for i in range(N + 1)))


# -----------------
# Abelian varieties
# -----------------

class SchemeMorphism_point_abelian_variety_field(AdditiveGroupElement, SchemeMorphism_point_projective_field):
    """
    A rational point of an abelian variety over a field.

    EXAMPLES::

        sage: # needs sage.schemes
        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: origin = E(0)
        sage: origin.domain()
        Spectrum of Rational Field
        sage: origin.codomain()
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    """
    pass
