import sage.all
from sage.schemes.toric.toric_subscheme import AlgebraicScheme_subscheme_toric
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.big_oh import O
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.real_mpfr import RR
from sage.functions.all import log
from sage.rings.integer import Integer

from weighted_projective_curve import WeightedProjectiveCurve

def is_HyperellipticCurveSmoothModel(C):
    """
    TODO
    """
    return isinstance(C, HyperellipticCurveSmoothModel_generic)

class HyperellipticCurveSmoothModel_generic(WeightedProjectiveCurve):
    def __init__(self, projective_model, f, h, genus):
        self._projective_model = projective_model
        self._genus = genus
        self._hyperelliptic_polynomials = (f, h)

        self._polynomial_ring = f.parent()
        self._base_ring = f.base_ring()

        # Some values which we will cache as a user asks for them
        self._alphas = None
        self._infinte_polynomials = None
        self._distinguished_point = None

        # TODO: is this simply genus + 1
        self._d = max(h.degree(), (f.degree() / Integer(2)).ceil())

        # Initalise the underlying curve
        A = self._projective_model.ambient_space()
        X = self._projective_model.defining_polynomials()
        WeightedProjectiveCurve.__init__(self, A, X)

    # TODO: _richcmp_ instead?
    def __richcmp__(self, other, op):
        return self._projective_model.__richcmp__(other._projective_model, op)

    def __ne__(self, other):
        return not self == other


    def _repr_(self):
        f, h = self._hyperelliptic_polynomials

        if h:
            if h.is_one():
                curve = f"y^2 + y = {f}"
            else:
                curve = f"y^2 + ({h})*y = {f}"
        else:
            curve = f"y^2 = {f}"

        # TODO:
        # The old class has these weird internal gens and then
        # printing polynomial rings to change output. This seems
        # dumb??
        # Will do something hacky here and we can talk about it.
        old_gen = str(self._polynomial_ring.gen())
        curve = curve.replace(old_gen, "x")
        return f"Hyperelliptic Curve over {self.base_ring()} defined by {curve}"

    def genus(self):
        """
        Compute the genus of the hyperelliptic curve
        """
        return self._genus

    def base_ring(self):
        """
        TODO
        """
        return self._base_ring

    def change_ring(self, R):
        """
        TODO
        """
        from hyperelliptic_constructor import HyperellipticCurveSmoothModel
        f, h = self._hyperelliptic_polynomials
        fR = f.change_ring(R)
        hR = h.change_ring(R)
        return HyperellipticCurveSmoothModel(fR, hR)

    base_extend = change_ring

    def polynomial_ring(self):
        """
        TODO
        """
        return self._polynomial_ring

    def point(self, coords, check=True):
        """
        TODO
        """
        if len(coords) == 2:
            X, Y = coords
            Z = self.base_ring().one()
        elif len(coords) == 3:
            X, Y, Z = coords
        else:
            raise ValueError("TODO")

        return self._projective_model.point([X, Y, Z], check=check)

    def hyperelliptic_polynomials(self):
        """
        Return the polynomials (f, h) such that
        C : y^2 + h*y = f
        """
        return self._hyperelliptic_polynomials

    def roots_at_infinity(self):
        """
        Compute the roots of: Y^2 + h[d]Y - d[2d] = 0
        When the curve is ramified, we expect one root, when
        the curve is split or inert we expect zero or two roots.
        """
        if self._alphas:
            return self._alphas

        f, h = self._hyperelliptic_polynomials
        x = f.parent().gen()
        d = self._d

        if h.is_zero():
            coeff = f[2 * d]
            # Handle the ramified case
            if coeff.is_zero():
                return [coeff]
            return f[2 * d].sqrt(all=True)

        self._alphas = (x**2 + x * h[d] - f[2 * d]).roots(multiplicities=False)
        return self._alphas

    def is_split(self):
        """
        Return True if the curve is split, i.e. there are two rational
        points at infinity.
        """
        return len(self.roots_at_infinity()) == 2

    def is_ramified(self):
        """
        Return True if the curve is ramified, i.e. there is one rational
        point at infinity.
        """
        return len(self.roots_at_infinity()) == 1

    def is_inert(self):
        """
        Return True if the curve is inert, i.e. there are no rational
        points at infinity.
        """
        return len(self.roots_at_infinity()) == 0

    def infinite_polynomials(self):
        """
        TODO: stupid name

        Computes G^±(x) for curves in the split degree model
        """
        if self._infinte_polynomials is not None:
            return self._infinte_polynomials

        alphas = self.roots_at_infinity()

        # This function only makes sense for the split model
        if not len(alphas) == 2:
            raise ValueError("hyperelliptic curve does not have the split model")

        f, h = self._hyperelliptic_polynomials
        alpha_plus, alpha_minus = alphas
        d = self._d

        # Construct G_plus from alpha_plus
        g = [None] * (d + 1)
        g[d] = alpha_plus
        for i in range(d - 1, -1, -1):
            # We need (g * (g + h))[x^(i + d)] to match f_{i + d}
            the_rest = g[d] * h[i] + sum(
                g[k] * (g[i + d - k] + h[i + d - k]) for k in range(i + 1, d)
            )
            g[i] = (f[i + d] - the_rest) / (2 * g[d] + h[d])

        G_plus = self._polynomial_ring(g)
        G_minus = -G_plus - h
        # Checks for the assumptions on G^±
        genus = self.genus()
        assert G_plus.degree() <= (genus + 1)
        assert (G_plus**2 + h * G_plus - f).degree() <= genus
        assert G_minus.leading_coefficient() == alpha_minus

        self._infinte_polynomials = G_plus, G_minus
        return self._infinte_polynomials

    def points_at_infinity(self):
        """
        Compute the points at infinity on the curve. Assumes we are using
        a weighted projective model for the curve
        """
        # TODO: check to False
        return [self.point([1, y, 0], check=True) for y in self.roots_at_infinity()]

    def is_x_coord(self, x):
        """
        Return True if ``x`` is the `x`-coordinate of a point on this curve.

        .. SEEALSO::

            See also :meth:`lift_x` to find the point(s) with a given
            `x`-coordinate.  This function may be useful in cases where
            testing an element of the base field for being a square is
            faster than finding its square root.

        INPUT:

        - ``x`` -- an element of the base ring of the curve

        OUTPUT:

        A bool stating whether or not `x` is a x-coordinate of a point on the curve

        EXAMPLES:

        When `x` is the `x`-coordinate of a rational point on the
        curve, we can request these::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(0)
            True

        There are no rational points with `x`-coordinate 3::

            sage: H.is_x_coord(3)
            False

        The function also handles the case when `h(x)` is not zero::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.is_x_coord(1)
            True

        We can perform these operations over finite fields too::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PolynomialRing(GF(163))
            sage: f = x^7 + x + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(13)
            True

        Including the case of characteristic two::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.is_x_coord(z4^3 + z4^2 + z4)
            True

        AUTHORS:

        - Giacomo Pope (2024): adapted from :meth:`lift_x`

        TESTS:

        The `x`-coordinate must be defined over the base field of the curve::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: p = 11
            sage: F = GF(11)
            sage: F_ext = GF(11^2)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(F_ext.gen())
            Traceback (most recent call last):
            ...
            TypeError: x must be coercible into the base ring of the curve
        """
        f, h = self.hyperelliptic_polynomials()
        K = self.base_ring()
        try:
            x = K(x)
        except (ValueError, TypeError):
            raise TypeError('x must be coercible into the base ring of the curve')

        # When h is zero then x is a valid coordinate if y2 is square
        if not h:
            y2 = f(x)
            return y2.is_square()
        # Generic case for h != 0
        a = f(x)
        b = h(x)
        # Special case for char 2
        if K.characteristic() == 2:
            R = f.parent()  # Polynomial ring K[x]
            F = R([-a, b, 1])
            return bool(F.roots())
        # Otherwise x is a point on the curve if the discriminant is a square
        D = b*b + 4*a
        return D.is_square()

    def lift_x(self, x, all=False):
        """
        Return one or all points with given `x`-coordinate.

        This method is deterministic: It returns the same data each
        time when called again with the same `x`.

        INPUT:

        - ``x`` -- an element of the base ring of the curve

        - ``all`` (bool, default ``False``) -- if ``True``, return a
          (possibly empty) list of all points; if ``False``, return
          just one point, or raise a :class:`ValueError` if there are none.

        OUTPUT:

        A point or list of up to two points on this curve.

        .. SEEALSO::

            :meth:`is_x_coord`

        AUTHORS:

        - Giacomo Pope (2024): Allowed for the case of characteristic two

        EXAMPLES:

        When `x` is the `x`-coordinate of a rational point on the
        curve, we can request these::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(0)
            [0 : -1 : 1]
            sage: H.lift_x(4, all=True)
            [[4 : -33 : 1], [4 : 33 : 1]]

        # TODO: these break in 10.2 because QQ(3, all=True) != [] but instead errors
        # Fixed in 10.3
        # There are no rational points with `x`-coordinate 3::

        #     sage: H.lift_x(3)
        #     Traceback (most recent call last):
        #     ...
        #     ValueError: No point with x-coordinate 3 on Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x^3 + 1

        # An empty list is returned when there are no points and ``all=True``::

        #     sage: H.lift_x(3, all=True)
        #     []

        The function also handles the case when `h(x)` is not zero::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(1)
            [1 : -3 : 1]

        We can perform these operations over finite fields too::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PolynomialRing(GF(163))
            sage: f = x^7 + x + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(13)
            [13 : 41 : 1]

        Including the case of characteristic two::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(z4^3 + z4^2 + z4, all=True)
            [[z4^3 + z4^2 + z4 : z4^2 + z4 + 1 : 1], [z4^3 + z4^2 + z4 : z4^3 : 1]]

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F1 = GF(11)
            sage: F2 = GF(13)
            sage: R.<x> = PolynomialRing(F1)
            sage: f = x^7 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(F2.random_element())
            Traceback (most recent call last):
            ...
            ValueError: x must have a common parent with the base ring

        Ensure that :issue:`37097` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(z4^3 + z4^2 + z4, all=True)
            [[z4^3 + z4^2 + z4 : z4^2 + z4 + 1 : 1], [z4^3 + z4^2 + z4 : z4^3 : 1]]
        """
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()

        f, h = self.hyperelliptic_polynomials()
        K = self.base_ring()

        # Compute the common parent between the base ring of the curve and
        # the parent of the input x-coordinate.
        try:
            L = cm.common_parent(x.parent(), K)
            x = L(x)
        except (TypeError, ValueError):
            raise ValueError('x must have a common parent with the base ring')

        # First we compute the y-coordinates the given x-coordinate
        ys = []
        one = L.one()

        # When h is zero we find all y-coordinates with a single sqrt
        if not h:
            y2 = f(x)
            # When y2 is not a square, ys will be an empty list
            ys = y2.sqrt(all=True, extend=False)
        # Otherwise we need roots of the discriminant
        else:
            a = f(x)
            b = h(x)
            # Special case for char 2
            if K.characteristic() == 2:
                R = f.parent()
                F = R([-a, b, 1])
                ys = F.roots(L, multiplicities=False)
            else:
                D = b*b + 4*a
                # When D is not a square, ys will be an empty list
                ys = [(-b+d)/2 for d in D.sqrt(all=True, extend=False)]

        if ys:
            ys.sort()  # Make lifting deterministic
            if all:
                return [self.point([x, y, one], check=False) for y in ys]
            else:
                return self.point([x, ys[0], one], check=False)

        if all:
            return []
        else:
            raise ValueError(f"No point with x-coordinate {x} on {self}")

    def distinguished_point(self):
        """
        Return the distinguished point of the hyperelliptic curve.
        By default, this is one of the points at infinity if possible.
        """
        if self._distinguished_point is None:
            if not self.is_inert():
                # For the the split and ramified case, a point at infinity is chosen,
                self._distinguished_point = self.points_at_infinity()[0]
            else:
                assert (
                    self.base_ring().characteristic() > 0
                ), "in characteristic 0, a distinguished_point needs to be specified"
                # in the inert case we choose a point with minimal x-coordinate
                for x0 in self.base_ring():
                    try:
                        self._distinguished_point = self.lift_x(x0)
                        break
                    except ValueError:
                        pass

        return self._distinguished_point

    def set_distinguished_point(self, P0):
        """
        Change the distinguished point of the hyperelliptic curve to P0.
        """
        assert isinstance(
            P0, AlgebraicScheme_subscheme_toric
        ), "the input has to be a point on the curve"
        self._distinguished_point = P0
        return None

    def __call__(self, *args):
        return self.point(args)

    def jacobian(self):
        from jacobian_generic import HyperellipticJacobian_generic
        return HyperellipticJacobian_generic(self)

    # -------------------------------------------
    # Hacky things
    # -------------------------------------------

    def is_singular(self):
        r"""
        Returns False, because hyperelliptic curves are smooth projective
        curves, as checked on construction.

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + 1)
            sage: H.is_singular()
            False
        """
        return False

    def is_smooth(self):
        r"""
        Returns True, because hyperelliptic curves are smooth projective
        curves, as checked on construction.

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 1)
            sage: H.is_smooth()
            True
        """
        return True

    # -------------------------------------------
    # Odd degree model functions
    # -------------------------------------------

    def odd_degree_model(self):
        r"""
        Return an odd degree model of self, or raise ValueError if one does not exist over the field of definition.

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: x = QQ['x'].gen()
            sage: H = HyperellipticCurveSmoothModel((x^2 + 2)*(x^2 + 3)*(x^2 + 5)); H
            Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + 10*x^4 + 31*x^2 + 30
            sage: H.odd_degree_model()
            Traceback (most recent call last):
            ...
            ValueError: No odd degree model exists over field of definition

            sage: K2 = QuadraticField(-2, 'a')                                          # needs sage.rings.number_field
            sage: Hp2 = H.change_ring(K2).odd_degree_model(); Hp2                       # needs sage.rings.number_field
            Hyperelliptic Curve over Number Field in a
             with defining polynomial x^2 + 2 with a = 1.414213562373095?*I
             defined by y^2 = 6*a*x^5 - 29*x^4 - 20*x^2 + 6*a*x + 1

            sage: K3 = QuadraticField(-3, 'b')                                          # needs sage.rings.number_field
            sage: Hp3 = H.change_ring(QuadraticField(-3, 'b')).odd_degree_model(); Hp3  # needs sage.rings.number_field
            Hyperelliptic Curve over Number Field in b
             with defining polynomial x^2 + 3 with b = 1.732050807568878?*I
             defined by y^2 = -4*b*x^5 - 14*x^4 - 20*b*x^3 - 35*x^2 + 6*b*x + 1

            Of course, ``Hp2`` and ``Hp3`` are isomorphic over the composite
            extension.  One consequence of this is that odd degree models
            reduced over "different" fields should have the same number of
            points on their reductions.  43 and 67 split completely in the
            compositum, so when we reduce we find:

            sage: # needs sage.rings.number_field
            sage: P2 = K2.factor(43)[0][0]
            sage: P3 = K3.factor(43)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849

            sage: H.change_ring(GF(43)).odd_degree_model().frobenius_polynomial()       # needs sage.rings.finite_rings
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849

            sage: # needs sage.rings.number_field
            sage: P2 = K2.factor(67)[0][0]
            sage: P3 = K3.factor(67)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489

            sage: H.change_ring(GF(67)).odd_degree_model().frobenius_polynomial()       # needs sage.rings.finite_rings
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489

        TESTS::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: HyperellipticCurveSmoothModel(x^5 + 1, 1).odd_degree_model()
            Traceback (most recent call last):
            ...
            NotImplementedError: odd_degree_model only implemented for curves in Weierstrass form

            sage: # TODO this used to have names="U, V"
            sage: HyperellipticCurveSmoothModel(x^5 + 1).odd_degree_model()
            Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + 1
        """
        f, h = self._hyperelliptic_polynomials
        if h:
            raise NotImplementedError("odd_degree_model only implemented for curves in Weierstrass form")
        if f.degree() % 2:
            # already odd, so just yield self
            return self

        rts = f.roots(multiplicities=False)
        if not rts:
            raise ValueError("No odd degree model exists over field of definition")
        rt = rts[0]
        x = f.parent().gen()
        fnew = f((x * rt + 1) / x).numerator()  # move rt to "infinity"

        from hyperelliptic_constructor import HyperellipticCurveSmoothModel
        return HyperellipticCurveSmoothModel(fnew, 0)

    def has_odd_degree_model(self):
        r"""
        Return True if an odd degree model of self exists over the field of definition; False otherwise.

        Use ``odd_degree_model`` to calculate an odd degree model.

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: x = QQ['x'].0
            sage: HyperellipticCurveSmoothModel(x^5 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurveSmoothModel(x^6 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurveSmoothModel(x^6 + x + 1).has_odd_degree_model()
            False
        """
        try:
            return bool(self.odd_degree_model())
        except ValueError:
            return False

    # -------------------------------------------
    # Magma (BOO!!) TODO
    # -------------------------------------------

    def _magma_init_(self, magma):
        """
        Internal function. Returns a string to initialize this elliptic
        curve in the Magma subsystem.

        EXAMPLES::

            sage: # optional - magma
            sage: R.<x> = QQ[]; C = HyperellipticCurveSmoothModel(x^3 + x - 1, x); C
            Hyperelliptic Curve over Rational Field
            defined by y^2 + x*y = x^3 + x - 1
            sage: magma(C)
            Hyperelliptic Curve defined by y^2 + x*y = x^3 + x - 1 over Rational Field
            sage: R.<x> = GF(9,'a')[]; C = HyperellipticCurveSmoothModel(x^3 + x - 1, x^10); C     # needs sage.rings.finite_rings
            Hyperelliptic Curve over Finite Field in a of size 3^2
            defined by y^2 + x^10*y = x^3 + x + 2
            sage: D = magma(C); D                                                       # needs sage.rings.finite_rings
            Hyperelliptic Curve defined by y^2 + (x^10)*y = x^3 + x + 2 over GF(3^2)
            sage: D.sage()                                                              # needs sage.rings.finite_rings
            Hyperelliptic Curve over Finite Field in a of size 3^2
            defined by y^2 + x^10*y = x^3 + x + 2
        """
        f, h = self._hyperelliptic_polynomials
        return 'HyperellipticCurve(%s, %s)' % (f._magma_init_(magma), h._magma_init_(magma))

    # -------------------------------------------
    # monsky washnitzer things... TODO
    # -------------------------------------------

    def monsky_washnitzer_gens(self):
        import monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        return S.gens()

    def invariant_differential(self):
        """
        Returns `dx/2y`, as an element of the Monsky-Washnitzer cohomology
        of self

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: C.invariant_differential()
            1 dx/2y

        """
        import monsky_washnitzer as m_w # TODO global import
        S = m_w.SpecialHyperellipticQuotientRing(self)
        MW = m_w.MonskyWashnitzerDifferentialRing(S)
        return MW.invariant_differential()


    def local_coordinates_at_nonweierstrass(self, P, prec=20, name='t'):
        """
        For a non-Weierstrass point `P = (a,b)` on the hyperelliptic
        curve `y^2 = f(x)`, return `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`,
        where `t = x - a` is the local parameter.

        INPUT:

        - ``P = (a, b)`` -- a non-Weierstrass point on self
        - ``prec`` --  desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = x - a`
        is the local parameter at `P`

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: P = H(1, 6)
            sage: x, y = H.local_coordinates_at_nonweierstrass(P, prec=5)
            sage: x
            1 + t + O(t^5)
            sage: y
            6 + t - 7/2*t^2 - 1/2*t^3 - 25/48*t^4 + O(t^5)
            sage: Q = H(-2, 12)
            sage: x, y = H.local_coordinates_at_nonweierstrass(Q, prec=5)
            sage: x
            -2 + t + O(t^5)
            sage: y
            12 - 19/2*t - 19/32*t^2 + 61/256*t^3 - 5965/24576*t^4 + O(t^5)

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        d = P[1]
        if d == 0:
            raise TypeError("P = %s is a Weierstrass point. Use local_coordinates_at_weierstrass instead!" % P)
        pol = self.hyperelliptic_polynomials()[0]
        L = PowerSeriesRing(self.base_ring(), name, default_prec=prec)
        t = L.gen()
        K = PowerSeriesRing(L, 'x')
        pol = K(pol)
        b = P[0]
        f = pol(t+b)
        for i in range((RR(log(prec)/log(2))).ceil()):
            d = (d + f/d)/2
        return t+b+O(t**(prec)), d + O(t**(prec))

    def local_coordinates_at_weierstrass(self, P, prec=20, name='t'):
        """
        For a finite Weierstrass point on the hyperelliptic
        curve `y^2 = f(x)`, returns `(x(t), y(t))` such that
        `(y(t))^2 = f(x(t))`, where `t = y` is the local parameter.

        INPUT:

        - ``P`` -- a finite Weierstrass point on self
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: `t`)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = y`
        is the local parameter at `P`

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: A = H(4, 0)
            sage: x, y = H.local_coordinates_at_weierstrass(A, prec=7)
            sage: x
            4 + 1/360*t^2 - 191/23328000*t^4 + 7579/188956800000*t^6 + O(t^7)
            sage: y
            t + O(t^7)
            sage: B = H(-5, 0)
            sage: x, y = H.local_coordinates_at_weierstrass(B, prec=5)
            sage: x
            -5 + 1/1260*t^2 + 887/2000376000*t^4 + O(t^5)
            sage: y
            t + O(t^5)

        AUTHOR:

            - Jennifer Balakrishnan (2007-12)
            - Francis Clarke (2012-08-26)
        """
        if P[1] != 0:
            raise TypeError("P = %s is not a finite Weierstrass point. Use local_coordinates_at_nonweierstrass instead!" % P)
        L = PowerSeriesRing(self.base_ring(), name)
        t = L.gen()
        pol = self.hyperelliptic_polynomials()[0]
        pol_prime = pol.derivative()
        b = P[0]
        t2 = t**2
        c = b + t2/pol_prime(b)
        c = c.add_bigoh(prec)
        for _ in range(int(1 + log(prec, 2))):
            c -= (pol(c) - t2)/pol_prime(c)
        return (c, t.add_bigoh(prec))

    def local_coordinates_at_infinity(self, prec=20, name='t'):
        """
        For the genus `g` hyperelliptic curve `y^2 = f(x)`, return
        `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`, where `t = y/x^{g+1}` is
        the local parameter at infinity.

        TODO/NOTE: In the previous implementation `t = x^g/y` was used.
        This is not a valid parameter on the smooth model, and the output is necessarily different.

        INPUT:

        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = y/x^{g+1}`
        is the local parameter at infinity

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 5*x^2 + 1)
            sage: x, y = H.local_coordinates_at_infinity(10)
            sage: x
            t^-2 - 5*t^4 + t^8 - 75*t^10 + O(t^12)
            sage: y
            t^-5 - 15*t + 3*t^5 - 150*t^7 + 90*t^11 + O(t^12)

        ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^3 - x + 1)
            sage: x, y = H.local_coordinates_at_infinity(10)
            sage: x
            t^-2 - t^2 + t^4 - 2*t^6 + 5*t^8 - 10*t^10 + O(t^12)
            sage: y
            t^-3 - 2*t + 2*t^3 - 3*t^5 + 8*t^7 - 15*t^9 + 42*t^11 + O(t^12)

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """

        if not self.is_ramified():
            raise NotImplementedError("only implemented for hyperelliptic curves with exactly one point at infinity")

        g = self.genus()
        pol, h = self.hyperelliptic_polynomials()
        if h:
            raise NotImplementedError("h = %s is nonzero" % h)
        K = LaurentSeriesRing(self.base_ring(), name, default_prec=prec+2)
        t = K.gen()
        L = PolynomialRing(K,'x')
        x = L.gen()
        i = 0
        w = (x**(g+1)*t)**2-pol
        wprime = w.derivative(x)
        x = t**-2
        for i in range((RR(log(prec+2)/log(2))).ceil()):
            x = x - w(x)/wprime(x)
        y = x**(g+1)*t
        return x+O(t**(prec+2)) , y+O(t**(prec+2))


    def local_coord(self, P, prec=20, name='t'):
        """
        Calls the appropriate local_coordinates function

        INPUT:

        - ``P`` -- a point on self
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))`, where `t`
        is the local parameter at `P`

        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: H.local_coord(H(1 ,6), prec=5)
            (1 + t + O(t^5), 6 + t - 7/2*t^2 - 1/2*t^3 - 25/48*t^4 + O(t^5))
            sage: H.local_coord(H(4, 0), prec=7)
            (4 + 1/360*t^2 - 191/23328000*t^4 + 7579/188956800000*t^6 + O(t^7), t + O(t^7))
            sage: H.local_coord(H(1, 0, 0), prec=5)
            (t^-2 - 23*t^2 + 18*t^4 - 1018*t^6 + O(t^7),
            t^-5 - 69*t^-1 + 54*t - 1467*t^3 + 3726*t^5 + O(t^6))

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """

        if P[2] == 0:
            return self.local_coordinates_at_infinity(prec, name)
        elif P[1] == 0:
            return self.local_coordinates_at_weierstrass(P, prec, name)
        else:
            return self.local_coordinates_at_nonweierstrass(P, prec, name)

    # TODO this is bad now we're in the smooth model

    # def rational_points(self, **kwds):
    #     r"""
    #     Find rational points on the hyperelliptic curve, all arguments are passed
    #     on to :meth:`sage.schemes.generic.algebraic_scheme.rational_points`.

    #     EXAMPLES:

    #     For the LMFDB genus 2 curve `932.a.3728.1 <https://www.lmfdb.org/Genus2Curve/Q/932/a/3728/1>`_::

    #         sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel # TODO Remove this after global import
    #         sage: R.<x> = PolynomialRing(QQ)
    #         sage: C = HyperellipticCurveSmoothModel(R([0, -1, 1, 0, 1, -2, 1]), R([1]))
    #         sage: C.rational_points(bound=8)
    #         [(-1 : -3 : 1),
    #         (-1 : 2 : 1),
    #         (0 : -1 : 1),
    #         (0 : 0 : 1),
    #         (0 : 1 : 0),
    #         (1/2 : -5/8 : 1),
    #         (1/2 : -3/8 : 1),
    #         (1 : -1 : 1),
    #         (1 : 0 : 1)]

    #     Check that :issue:`29509` is fixed for the LMFDB genus 2 curve
    #     `169.a.169.1 <https://www.lmfdb.org/Genus2Curve/Q/169/a/169/1>`_::

    #         sage: C = HyperellipticCurveSmoothModel(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]))
    #         sage: C.rational_points(bound=10)
    #         [(-1 : 0 : 1),
    #         (-1 : 1 : 1),
    #         (0 : -1 : 1),
    #         (0 : 0 : 1),
    #         (0 : 1 : 0)]

    #     An example over a number field::

    #         sage: R.<x> = PolynomialRing(QuadraticField(2))                             # needs sage.rings.number_field
    #         sage: C = HyperellipticCurveSmoothModel(R([1, 0, 0, 0, 0, 1]))              # needs sage.rings.number_field
    #         sage: C.rational_points(bound=2)                                            # needs sage.rings.number_field
    #         [(-1 : 0 : 1),
    #          (0 : -1 : 1),
    #          (0 : 1 : 0),
    #          (0 : 1 : 1),
    #          (1 : -a : 1),
    #          (1 : a : 1)]
    #     """
    #     from sage.schemes.curves.constructor import Curve
    #     # we change C to be a plane curve to allow the generic rational
    #     # points code to reduce mod any prime, whereas a HyperellipticCurve
    #     # can only be base changed to good primes.
    #     C = self
    #     if 'F' in kwds:
    #         C = C.change_ring(kwds['F'])

    #     return [C(pt) for pt in Curve(self).rational_points(**kwds)]
