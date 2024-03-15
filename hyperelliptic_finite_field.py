from sage.misc.prandom import choice
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

import hyperelliptic_generic

from sage.rings.real_mpfr import RR
from sage.arith.misc import binomial
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
from sage.libs.pari.all import pari

class HyperellipticCurveSmoothModel_finite_field(
    hyperelliptic_generic.HyperellipticCurveSmoothModel_generic
):
    def __init__(self, projective_model, f, h, genus):
        super().__init__(projective_model, f, h, genus)

    def random_point(self):
        """
        Return a random point on this hyperelliptic curve, uniformly chosen
        among all rational points.
        """
        k = self.base_ring()
        n = 2 * k.order() + 1

        while True:
            # Choose the point at infinity with probability 1/(2q + 1)
            i = ZZ.random_element(n)
            if not i and not self.is_inert():
                # TODO: deal with that there is more than one point at infinity
                return choice(self.points_at_infinity())
            v = self.lift_x(k.random_element(), all=True)
            try:
                return v[i % 2]
            except IndexError:
                pass

    @cached_method
    def points(self):
        """
        TODO: couldn't be more stupid
        """
        # TODO: this is very silly but works
        points = self.points_at_infinity()
        for x in self.base_ring():
            points.extend(self.lift_x(x, all=True))
        return points

    @cached_method
    def cardinality(self):
        """
        TODO: couldn't be more stupid
        """
        return len(self.points())

    order = cardinality

    # -------------------------------------------
    # Frobenius Computations
    # -------------------------------------------

    def _frobenius_coefficient_bound_charpoly(self):
        r"""
        Computes bound on number of `p`-adic digits needed to recover
        frobenius polynomial computing the characteristic polynomial
        of the frobenius matrix, i.e. returns `B` so that knowledge of
        `a_1`, ..., `a_g` modulo `p^B` determine frobenius polynomial
        uniquely.

        The bound used here stems from the expression of the coefficients
        of the characteristic polynomial of the Frobenius as sums
        of products of its eigenvalues:

        .. MATH::

            \| a_i \| \leq \binom{2g}{i} \sqrt{q}^i

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_charpoly()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_charpoly()
            3

            sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_charpoly()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^9 + t + 1)._frobenius_coefficient_bound_charpoly()
            3
            sage: HyperellipticCurve(t^11 + t + 1)._frobenius_coefficient_bound_charpoly()
            3
            sage: HyperellipticCurve(t^13 + t + 1)._frobenius_coefficient_bound_charpoly()
            4
        """
        assert self.base_ring().is_finite()
        p = self.base_ring().characteristic()
        q = self.base_ring().order()
        sqrtq = RR(q).sqrt()
        g = self.genus()

        # note: this bound is from Kedlaya's paper, but he tells me it's not
        # the best possible
        M = 2 * binomial(2*g, g) * sqrtq**g
        B = ZZ(M.ceil()).exact_log(p)
        if p**B < M:
            B += 1
        return B


    def _frobenius_coefficient_bound_traces(self, n=1):
        r"""
        Computes bound on number of `p`-adic digits needed to recover
        the number of rational points on `n` extensions computing
        traces of the frobenius matrix powers, i.e. returns `B` so that
        knowledge of `\tr(M^1)`, ..., `\tr(M^n)` modulo `p^B` determine
        `N_1`, ..., `N_n` uniquely.

        The formula stems from the expression of the trace of the Frobenius
        as a sum of `i`-th powers of its eigenvalues.

        .. MATH::

            \| \tr(M^i) \| \leq 2 g \sqrt{q}^i

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_traces()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_traces()
            2

            sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^9 + t + 1)._frobenius_coefficient_bound_traces(n=3)
            2
            sage: HyperellipticCurve(t^11 + t + 1)._frobenius_coefficient_bound_traces(n=3)
            2
            sage: HyperellipticCurve(t^13 + t + 1)._frobenius_coefficient_bound_traces(n=5)
            3

            sage: R.<t> = PolynomialRing(GF(11))
            sage: H = HyperellipticCurve(t^5 - t + 1)
            sage: H._frobenius_coefficient_bound_traces()
            2
        """
        p = self.base_ring().characteristic()
        q = self.base_ring().order()
        sqrtq = RR(q).sqrt()
        g = self.genus()

        M = 4 * g * sqrtq**n
        B = ZZ(M.ceil()).exact_log(p)
        if p**B < M:
            B += 1
        return B

    def frobenius_matrix_hypellfrob(self, N=None):
        r"""
        Compute `p`-adic frobenius matrix to precision `p^N`.
        If `N` not supplied, a default value is selected, which is the
        minimum needed to recover the charpoly unambiguously.

        .. note::

            Implemented using ``hypellfrob``, which means it only works
            over the prime field `GF(p)`, and requires `p > (2g+1)(2N-1)`.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_matrix_hypellfrob()
            [1258 + O(37^2)  925 + O(37^2)  132 + O(37^2)  587 + O(37^2)]
            [1147 + O(37^2)  814 + O(37^2)  241 + O(37^2) 1011 + O(37^2)]
            [1258 + O(37^2) 1184 + O(37^2) 1105 + O(37^2)  482 + O(37^2)]
            [1073 + O(37^2)  999 + O(37^2)  772 + O(37^2)  929 + O(37^2)]

        The ``hypellfrob`` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_matrix_hypellfrob()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented
            for hyperelliptic curves defined over prime fields.

        nor too small characteristic::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_matrix_hypellfrob()
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
        """
        p = self.base_ring().characteristic()
        e = self.base_ring().degree()
        if e != 1:
            raise NotImplementedError("Computation of Frobenius matrix only implemented for hyperelliptic curves defined over prime fields.")

        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            # need y^2 = f(x)
            raise NotImplementedError("only implemented for curves y^2 = f(x)")

        sign = 1
        if not f.is_monic():
            # at this time we need a monic f
            c = f.leading_coefficient()
            f = f / c
            if c.is_square():
                # solutions of $y^2 = c * f(x)$ correspond naturally to
                # solutions of $(sqrt(c) y)^2 = f(x)$
                pass
            else:
                # we'll count points on a twist and then correct by untwisting...
                sign = -1
        assert f.is_monic()

        # By default, use precision enough to be able to compute the
        # frobenius minimal polynomial
        if N is None:
            N = self._frobenius_coefficient_bound_charpoly()

        matrix_of_frobenius = hypellfrob(p, N, f)
        matrix_of_frobenius = sign * matrix_of_frobenius
        return matrix_of_frobenius


    def frobenius_matrix(self, N=None, algorithm='hypellfrob'):
        r"""
        Compute `p`-adic frobenius matrix to precision `p^N`.
        If `N` not supplied, a default value is selected, which is the
        minimum needed to recover the charpoly unambiguously.

        .. note::

            Currently only implemented using ``hypellfrob``,
            which means it only works over the prime field `GF(p)`,
            and requires `p > (2g+1)(2N-1)`.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_matrix()
            [1258 + O(37^2)  925 + O(37^2)  132 + O(37^2)  587 + O(37^2)]
            [1147 + O(37^2)  814 + O(37^2)  241 + O(37^2) 1011 + O(37^2)]
            [1258 + O(37^2) 1184 + O(37^2) 1105 + O(37^2)  482 + O(37^2)]
            [1073 + O(37^2)  999 + O(37^2)  772 + O(37^2)  929 + O(37^2)]

        The ``hypellfrob`` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented
            for hyperelliptic curves defined over prime fields.

        nor too small characteristic::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
        """
        if algorithm != 'hypellfrob':
            raise ValueError("Unknown algorithm")

        # By default, use precision enough to be able to compute the
        # frobenius minimal polynomial
        if N is None:
            N = self._frobenius_coefficient_bound_charpoly()

        return self.frobenius_matrix_hypellfrob(N=N)

    def frobenius_polynomial_cardinalities(self, a=None):
        r"""
        Compute the charpoly of frobenius, as an element of `\ZZ[x]`,
        by computing the number of points on the curve over `g` extensions
        of the base field where `g` is the genus of the curve.

        .. WARNING::

            This is highly inefficient when the base field or the genus of the
            curve are large.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial_cardinalities()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial_cardinalities()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Curve over a non-prime field::

            sage: K.<z> = GF(7**2)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z^2)
            sage: H.frobenius_polynomial_cardinalities()
            x^4 + 8*x^3 + 70*x^2 + 392*x + 2401

        This method may actually be useful when ``hypellfrob`` does not work::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_polynomial_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
            sage: H.frobenius_polynomial_cardinalities()
            x^8 - 5*x^7 + 7*x^6 + 36*x^5 - 180*x^4 + 252*x^3 + 343*x^2 - 1715*x + 2401
        """
        g = self.genus()
        q = self.base_ring().cardinality()

        if a is None:
            # this may actually call frobenius_polynomial()
            a = self.count_points(g)
            # maybe calling count_points_exhaustive() would make more sense
            # but the method is currently only called with a precomputed list
            # of number of points so it does not really matter

        # computation of the reciprocal polynomial
        s = [ai - q**(i+1) - 1 for i, ai in enumerate(a)]
        coeffs = [1]
        for i in range(1, g + 1):
            c = 0
            for j in range(i):
                c += s[i-1-j]*coeffs[j]
            coeffs.append(c/i)
        coeffs = coeffs + [coeffs[g-i] * q**(i) for i in range(1, g + 1)]

        return ZZ['x'](coeffs).reverse()

    def frobenius_polynomial_matrix(self, M=None, algorithm='hypellfrob'):
        r"""
        Compute the charpoly of frobenius, as an element of `\ZZ[x]`,
        by computing the charpoly of the frobenius matrix.

        This is currently only supported when the base field is prime
        and large enough using the ``hypellfrob`` library.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial_matrix()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial_matrix()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Curves defined over larger prime fields::

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^5 + 1)
            sage: H.frobenius_polynomial_matrix()
            x^8 + 281*x^7 + 55939*x^6 + 14144175*x^5 + 3156455369*x^4 + 707194605825*x^3
            + 139841906155939*x^2 + 35122892542149719*x + 6249500014999800001
            sage: H = HyperellipticCurve(t^15 + t^5 + 1)
            sage: H.frobenius_polynomial_matrix()  # long time, 8s on a Corei7
            x^14 - 76*x^13 + 220846*x^12 - 12984372*x^11 + 24374326657*x^10 - 1203243210304*x^9
            + 1770558798515792*x^8 - 74401511415210496*x^7 + 88526169366991084208*x^6
            - 3007987702642212810304*x^5 + 3046608028331197124223343*x^4
            - 81145833008762983138584372*x^3 + 69007473838551978905211279154*x^2
            - 1187357507124810002849977200076*x + 781140631562281254374947500349999

        This ``hypellfrob`` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_polynomial_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented
            for hyperelliptic curves defined over prime fields.
        """
        K = self.base_ring()
        p = K.characteristic()
        q = K.cardinality()
        g = self.genus()
        N = self._frobenius_coefficient_bound_charpoly()
        # compute charpoly over ZZ and then reduce back
        # (because charpoly of p-adic matrices sometimes loses precision)
        M = self.frobenius_matrix(N=N, algorithm=algorithm).change_ring(ZZ)

        # get a_g, ..., a_0 in ZZ (i.e. with correct signs)
        f = M.charpoly().list()[g:2*g+1]
        ppow = p**N
        f = [x % ppow for x in f]
        f = [x if 2*x < ppow else x - ppow for x in f]

        # get a_{2g}, ..., a_{g+1}
        f = [f[g-i] * q**(g-i) for i in range(g)] + f

        return ZZ['x'](f)

    def frobenius_polynomial_pari(self):
        r"""
        Compute the charpoly of frobenius, as an element of `\ZZ[x]`,
        by calling the PARI function ``hyperellcharpoly``.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial_pari()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial_pari()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Slightly larger example::

            sage: K = GF(2003)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 487*t^5 + 9*t + 1)
            sage: H.frobenius_polynomial_pari()
            x^6 - 14*x^5 + 1512*x^4 - 66290*x^3 + 3028536*x^2 - 56168126*x + 8036054027

        Curves defined over a non-prime field are supported as well::

            sage: K.<a> = GF(7^2)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + a*t + 1)
            sage: H.frobenius_polynomial_pari()
            x^4 + 4*x^3 + 84*x^2 + 196*x + 2401

            sage: K.<z> = GF(23**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^3 + z*t + 4)
            sage: H.frobenius_polynomial_pari()
            x^2 - 15*x + 12167

        Over prime fields of odd characteristic, `h` may be non-zero::

            sage: K = GF(101)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + 27*t + 3, t)
            sage: H.frobenius_polynomial_pari()
            x^4 + 2*x^3 - 58*x^2 + 202*x + 10201

        TESTS:

        Check that :issue:`28789` is fixed::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: u = x^10 + x^9 + x^8 + x
            sage: C = HyperellipticCurve(u)
            sage: C.frobenius_polynomial_pari()
            x^8 + 2*x^7 + 6*x^6 + 9*x^5 + 18*x^4 + 27*x^3 + 54*x^2 + 54*x + 81
        """
        f, h = self.hyperelliptic_polynomials()
        return ZZ['x'](pari([f, h]).hyperellcharpoly())

    @cached_method
    def frobenius_polynomial(self):
        r"""
        Compute the charpoly of frobenius, as an element of `\ZZ[x]`.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Slightly larger example::

            sage: K = GF(2003)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 487*t^5 + 9*t + 1)
            sage: H.frobenius_polynomial()
            x^6 - 14*x^5 + 1512*x^4 - 66290*x^3 + 3028536*x^2 - 56168126*x + 8036054027

        Curves defined over a non-prime field of odd characteristic,
        or an odd prime field which is too small compared to the genus,
        are supported via PARI::

            sage: K.<z> = GF(23**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^3 + z*t + 4)
            sage: H.frobenius_polynomial()
            x^2 - 15*x + 12167

            sage: K.<z> = GF(3**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z**3)
            sage: H.frobenius_polynomial()
            x^4 - 3*x^3 + 10*x^2 - 81*x + 729

        Over prime fields of odd characteristic, `h` may be non-zero::

            sage: K = GF(101)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + 27*t + 3, t)
            sage: H.frobenius_polynomial()
            x^4 + 2*x^3 - 58*x^2 + 202*x + 10201

        Over prime fields of odd characteristic, `f` may have even degree::

            sage: H = HyperellipticCurve(t^6 + 27*t + 3)
            sage: H.frobenius_polynomial()
            x^4 + 25*x^3 + 322*x^2 + 2525*x + 10201

        In even characteristic, the naive algorithm could cover all cases
        because we can easily check for squareness in quotient rings of
        polynomial rings over finite fields but these rings unfortunately
        do not support iteration::

            sage: K.<z> = GF(2**5)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z**3, t)
            sage: H.frobenius_polynomial()
            x^4 - x^3 + 16*x^2 - 32*x + 1024

        TESTS:

        Check that :issue:`28789` is fixed::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: u = x^10 + x^9 + x^8 + x
            sage: C = HyperellipticCurve(u)
            sage: C.frobenius_polynomial()
            x^8 + 2*x^7 + 6*x^6 + 9*x^5 + 18*x^4 + 27*x^3 + 54*x^2 + 54*x + 81
        """
        K = self.base_ring()
        e = K.degree()
        q = K.cardinality()

        g = self.genus()
        f, h = self.hyperelliptic_polynomials()

        if (e == 1 and
            q >= (2*g+1)*(2*self._frobenius_coefficient_bound_charpoly()-1) and
            h == 0 and f.degree() % 2):
            return self.frobenius_polynomial_matrix()
        elif q % 2 == 1:
            return self.frobenius_polynomial_pari()
        else:
            return self.frobenius_polynomial_cardinalities()
