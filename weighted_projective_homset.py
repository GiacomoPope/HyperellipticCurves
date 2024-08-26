"""
Hom-sets of weighted projective schemes

AUTHORS:

- Gareth Ma (2024): initial version, based on unweighted version.
"""

# *****************************************************************************
#        Copyright (C) 2024 Gareth Ma <grhkm21@gmail.com>
#        Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#        Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from copy import copy

from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RR
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.schemes.generic.homset import SchemeHomset_points

lazy_import(
    "sage.schemes.projective.projective_rational_point",
    [
        "enum_projective_finite_field",
        "enum_projective_number_field",
        "enum_projective_rational_field",
        "sieve",
    ],
)

# TODO: remember to remove this
# this exists to make my linter happy
from sage.schemes.projective.projective_rational_point import (
    enum_projective_finite_field,
    enum_projective_number_field,
    enum_projective_rational_field,
    sieve,
)


class SchemeHomset_points_weighted_projective_ring(SchemeHomset_points):
    """
    Set of rational points of a weighted projective variety over a ring.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: # TODO
    """

    def points(self, **kwds):
        """
        Return some or all rational points of this weighted projective scheme.

        For dimension 0 subschemes points are determined through a groebner
        basis calculation. For schemes or subschemes with dimension greater than 1
        points are determined through enumeration up to the specified bound.

        INPUT:

        kwds:

        - ``bound`` -- real number (default: 0). The bound for the coordinates for
          subschemes with dimension at least 1.

        - ``precision`` -- integer (default: 53). The precision to use to
          compute the elements of bounded height for number fields.

        - ``point_tolerance`` -- positive real number (default: `10^{-10}`).
          For numerically inexact fields, two points are considered the same
          if their coordinates are within tolerance.

        - ``zero_tolerance`` -- positive real number (default: `10^{-10}`).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        - ``tolerance`` -- a rational number in (0,1] used in doyle-krumm algorithm-4
          for enumeration over number fields.

        OUTPUT:

        - a list of rational points of a projective scheme

        .. WARNING::

            For numerically inexact fields such as ComplexField or RealField the
            list of points returned is very likely to be incomplete. It may also
            contain repeated points due to tolerances.

        EXAMPLES::

            sage: P = ProjectiveSpace(QQ, 2, "x, y")
            sage: P(QQ).points(bound=4)
            [(-4 : 1), (-3 : 1), (-2 : 1), (-3/2 : 1), (-4/3 : 1), (-1 : 1),
             (-3/4 : 1), (-2/3 : 1), (-1/2 : 1), (-1/3 : 1), (-1/4 : 1), (0 : 1),
             (1/4 : 1), (1/3 : 1), (1/2 : 1), (2/3 : 1), (3/4 : 1), (1 : 0), (1 : 1),
             (4/3 : 1), (3/2 : 1), (2 : 1), (3 : 1), (4 : 1)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)                                          # needs sage.rings.number_field
            sage: WP.<x,y,z> = WeightedProjectiveSpace(K, [1, 3, 1])                    # needs sage.rings.number_field
            sage: len(P(K).points(bound=1.8))                                           # needs sage.rings.number_field
            309
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: len(P(K).points(bound=1.8))
            309

        ::

            sage: P1 = ProjectiveSpace(GF(2), 1)
            sage: F.<a> = GF(4, 'a')                                                    # needs sage.rings.finite_rings
            sage: P1(F).points()                                                        # needs sage.libs.singular sage.rings.finite_rings
            [(0 : 1), (1 : 0), (1 : 1), (a : 1), (a + 1 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([(y^3-y*z^2) - (x^3-x*z^2), (y^3-y*z^2) + (x^3-x*z^2)])
            sage: E(P.base_ring()).points()                                             # needs sage.libs.singular
            [(-1 : -1 : 1), (-1 : 0 : 1), (-1 : 1 : 1), (0 : -1 : 1), (0 : 0 : 1),
             (0 : 1 : 1), (1 : -1 : 1), (1 : 0 : 1), (1 : 1 : 1)]

        ::

            sage: # needs sage.rings.real_mpfr
            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: L = E(P.base_ring()).points(); sorted(L, key=str)                     # needs sage.libs.singular
            verbose 0 (...: projective_homset.py, points) Warning: computations in
            the numerical fields are inexact;points may be computed partially or incorrectly.
            [(-0.500000000000000 + 0.866025403784439*I : 1.00000000000000 : 0.000000000000000),
             (-0.500000000000000 - 0.866025403784439*I : 1.00000000000000 : 0.000000000000000),
             (-1.00000000000000*I : 0.000000000000000 : 1.00000000000000),
             (0.000000000000000 : 0.000000000000000 : 1.00000000000000),
             (1.00000000000000 : 1.00000000000000 : 0.000000000000000),
             (1.00000000000000*I : 0.000000000000000 : 1.00000000000000)]
            sage: L[0].codomain()                                                       # needs sage.libs.singular
            Projective Space of dimension 2 over Complex Field with 53 bits of precision

        ::

            sage: # needs sage.rings.complex_double
            sage: P.<x,y,z> = ProjectiveSpace(CDF, 2)
            sage: E = P.subscheme([y^2 + x^2 + z^2, x*y*z])
            sage: len(E(P.base_ring()).points())                                        # needs sage.libs.singular
            verbose 0 (...: projective_homset.py, points) Warning: computations in
            the numerical fields are inexact;points may be computed partially or incorrectly.
            6
        """
        from weighted_projective_space import WeightedProjectiveSpace_ring

        X = self.codomain()
        if not isinstance(X, WeightedProjectiveSpace_ring) and X.base_ring() in Fields():
            if hasattr(X.base_ring(), "precision"):
                numerical = True
                verbose(
                    "Warning: computations in the numerical fields are inexact;points may be"
                    " computed partially or incorrectly.",
                    level=0,
                )
                pt_tol = RR(kwds.pop("point_tolerance", 10 ** (-10)))
                zero_tol = RR(kwds.pop("zero_tolerance", 10 ** (-10)))
                if pt_tol <= 0 or zero_tol <= 0:
                    raise ValueError("tolerance must be positive")
            else:
                numerical = False

            # Then it must be a subscheme
            dim_ideal = X.defining_ideal().dimension()

            # if X has no points
            if dim_ideal < 1:
                return []

            # if X is zero-dimensional
            if dim_ideal == 1:
                rat_points = set()
                PS = X.ambient_space()
                N = PS.dimension_relative()
                BR = X.base_ring()
                # need a lexicographic ordering for elimination
                R = PolynomialRing(BR, N + 1, PS.variable_names(), order="lex")
                I = R.ideal(X.defining_polynomials())
                I0 = R.zero_ideal()

                # Determine the points through elimination
                # This is much faster than using the I.variety() function on each affine chart.
                for k, R_gen in enumerate(R.gens()):
                    # create the elimination ideal for the kth affine patch
                    G = I.substitute({R_gen: 1}).groebner_basis()
                    if G == [1]:
                        continue

                    # We have a nontrivial Groebner basis
                    P = {R_gen: 1}
                    points = [P]
                    # work backwards from solving each equation for the possible
                    # values of the next coordinate
                    for i in reversed(range(len(G))):
                        new_points = []
                        good = 0
                        for P in points:
                            # substitute in our dictionary entry that has the values
                            # of coordinates known so far. This results in a single
                            # variable polynomial (by elimination)
                            L = G[i].substitute(P)
                            if R(L).degree() <= 0:
                                new_points.append(P)
                                good = 1
                                continue

                            # R(L) is not constant
                            if numerical:
                                for pol in L.univariate_polynomial().roots(multiplicities=False):
                                    good = 1
                                    r = L.variables()[0]
                                    varindex = R.gens().index(r)
                                    P[R.gen(varindex)] = pol
                                    new_points.append(copy(P))
                            else:
                                L = L.factor()
                                # the linear factors give the possible rational values of
                                # this coordinate
                                for pol, _ in L:
                                    if pol.degree() == 1 and len(pol.variables()) == 1:
                                        good = 1
                                        r = pol.variables()[0]
                                        varindex = R.gens().index(r)
                                        # add this coordinates information to
                                        # each dictionary entry
                                        key = R.gen(varindex)
                                        val = (
                                            -pol.constant_coefficient()
                                            / pol.monomial_coefficient(r)
                                        )
                                        P[key] = val
                                        new_points.append(copy(P))
                        if good:
                            points = new_points

                    # the dictionary entries now have values for all coordinates
                    # they are the rational solutions to the equations
                    # make them into projective points
                    for i, pt in enumerate(points):
                        if numerical:
                            if len(pt) == N + 1:
                                S = PS([pt[R.gen(j)] for j in range(N + 1)])
                                S.normalize_coordinates()
                                if all(g(list(S)) < zero_tol for g in X.defining_polynomials()):
                                    rat_points.add(S)
                        else:
                            if len(pt) == N + 1 and I.subs(pt) == I0:
                                S = X([pt[R.gen(j)] for j in range(N + 1)])
                                S.normalize_coordinates()
                                rat_points.add(S)

                # remove duplicate element using tolerance
                if numerical:
                    dupl_points = list(rat_points)
                    rat_points = set()
                    for u in dupl_points:
                        for v in rat_points:
                            if all((u[k] - v[k]).abs() < pt_tol for k in range(len(u))):
                                break
                        else:
                            rat_points.add(u)

                return sorted(rat_points)

        R = self.value_ring()
        B = kwds.pop("bound", 0)
        tol = kwds.pop("tolerance", 1e-2)
        prec = kwds.pop("precision", 53)
        if isinstance(R, RationalField):
            if not B > 0:
                raise TypeError(f"a positive bound B (={B}) must be specified")

            if isinstance(X, AlgebraicScheme_subscheme):
                return sieve(X, B)

            return enum_projective_rational_field(self, B)

        if R in NumberFields():
            if not B > 0:
                raise TypeError(f"a positive bound B (={B}) must be specified")

            return enum_projective_number_field(self, bound=B, tolerance=tol, precision=prec)

        if isinstance(R, FiniteField):
            return enum_projective_finite_field(self.extended_codomain())

        raise TypeError("unable to enumerate points over %s" % R)
