from sage.rings.polynomial.polynomial_element import Polynomial
from hyperelliptic_split import HyperellipticCurveSplit, HyperellipticPoint
from sage.rings.integer import Integer
from sage.misc.prandom import randint
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.misc.prandom import choice
from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method

# Needed until https://github.com/sagemath/sage/pull/37118 is merged.
from uniform_random_sampling import uniform_random_polynomial

class JacobianSplit:
    def __init__(self, H):
        if not isinstance(H, HyperellipticCurveSplit):
            raise ValueError("TODO")
        self._curve = H
        self._element = MumfordDivisorSplit

    def __repr__(self) -> str:
        return f"Jacobian of {self._curve}"

    def curve(self):
        return self._curve

    def base_ring(self):
        return self._curve.base_ring()

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self._curve.genus()
        R = self._curve.polynomial_ring()
        n = (g/2).ceil()
        return self._element(self, R.one(), R.zero(), n)

    def __call__(self, *args):
        if isinstance(args[0], HyperellipticPoint):
            R, x = self._curve.polynomial_ring().objgen()
            g = self._curve.genus()
            P = args[0]
            [X,Y,Z] = P.coords()
            # TODO:
            # we use the embedding P \mapsto P - \infty_+
            # note: P - \inft_+ = P + n*\infty_+ + m*\infty_- - D_\infty,
            # where n = ((g-1)/2).floor() 
            # do we want a different embedding by default?
            n = ((g-1)/2).floor()
            if Z == 0:
                alpha = Y/X
                if alpha == self._curve._alphas[0]:
                    n = n + 1
                u = R(1)
                v = R(0)
            else:
                u = x - X
                v = R(Y)
        elif len(args) == 1:
            return self.zero()
        elif len(args) == 2:
            u, v = args
            n = 0
        elif len(args) == 3:
            u, v, n = args
        else:
            raise NotImplementedError
        return self._element(self, u, v, n=n)

    @cached_method
    def cardinality(self):
        """
        TODO: currently using lazy methods by calling sage
        """
        return sum(self.curve().frobenius_polynomial())

    order = cardinality

    def _random_element_cover(self, degree=None):
        r"""
        Returns a random element from the Jacobian. Distribution is not
        uniformly random, but returns the entire group.
        """
        H = self.curve()
        K = H.base_ring()
        R = H.polynomial_ring()
        g = H.genus()

        f, h = H.hyperelliptic_polynomials()

        if degree is None:
            # until the random sampling is uniform (open PR for sage)
            # this is bad.
            degree = (-1, g)
            # degree = g # TODO: for odd genus this may need to be g + 1

        while True:
            u = uniform_random_polynomial(R, degree=degree)
            if u.is_zero():
                n = randint(0, g)
                return self(R(1),R(0),n)
            u = u.monic()
            try:
                D = self.zero()
                for x, e in u.factor():
                    # Solve y^2 + hy - f = 0 mod x
                    from sage.rings.polynomial.polynomial_ring import polygen
                    K_ext = K.extension(modulus=x, names="a")
                    y_ext = polygen(K_ext, "y_ext")
                    h_ = K_ext(h % x)
                    f_ = K_ext(f % x)
                    y = choice((y_ext**2 + h_ * y_ext - f_).roots(multiplicities=False))
                    try:
                        # Quotient ring elements
                        y = y.lift()
                    except AttributeError:
                        pass

                    try:
                        g = self.curve().genus()
                        for _ in range(e):
                            n = randint(0, g - x.degree())
                            D += self(x, R(y), n)
                    except (ValueError, TypeError):
                        raise IndexError

                return D

            except IndexError:
                pass

    def _random_element_rational(self):
        r"""
        Returns a random element from the Jacobian. It does not necessarily
        return the entire group.
        """
        H = self.curve()
        g = H.genus()

        # We randomly sample 2g + 1 points on the hyperelliptic curve
        points = [H.random_point() for _ in range(2*g + 1)]

        # We create 2g + 1 divisors of the form (P) - infty
        divisors = [self(P) for P in points if P[2] != 0]

        # If we happened to only sample the point at infinity, we return this
        # Otherwise we compute the sum of all divisors.
        if not divisors:
            return self.zero()
        return sum(divisors, start=self.zero())

    def random_element(self, fast=True, *args, **kwargs):
        r"""
        Returns a random element from the Jacobian. Distribution is not
        uniformly random.

        INPUT:

        - ``fast`` -- (boolean, default ``True``) If set to ``True``, a fast
          algorithm is used, but the output is **NOT** guaranteed to cover the
          entire Jacobian. See examples below. If set to ``False``, a slower
          algorithm is used, but covers the entire Jacobian.
        """
        if not isinstance(self.base_ring(), FiniteField_generic):
            raise NotImplementedError("random element of Jacobian is only implemented over Finite Fields")

        if fast:
            return self._random_element_rational(*args, **kwargs)
        return self._random_element_cover(*args, **kwargs)


    def cantor_composition(self, u1, v1, u2, v2):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()

        # Ensure D1 and D2 are semi-reduced divisors
        assert (
            v1.degree() < u1.degree() and v2.degree() < u2.degree()
        ), "The degree of bi must be smaller than ai"
        assert (
            u1.degree() <= f.degree() and u2.degree() <= f.degree()
        ), f"The degree of ai must be smaller than f, {u1.degree()}, {u2.degree()}"
        assert (v1**2 + v1 * h - f) % u1 == 0, "D1 is not a valid divisor"
        assert (v2**2 + v2 * h - f) % u2 == 0, "D2 is not a valid divisor"

        # TODO: when gcd(u0, u1) == 1 we can speed this up
        # TODO: when h = 0 we can speed this up.

        # Step One
        s0, e1, e2 = u1.xgcd(u2)
        assert s0 == e1 * u1 + e2 * u2

        # Step Two
        s, c1, c2 = s0.xgcd(v1 + v2 + h)
        assert s == c1 * s0 + c2 * (v1 + v2 + h)

        # Step Three
        f1 = c1 * e1
        f2 = c1 * e2
        f3 = c2
        assert s == f1 * u1 + f2 * u2 + f3 * (v1 + v2 + h)

        # Step Four
        u3 = (u1 * u2) // (s**2)
        v3 = (f1 * u1 * v2 + f2 * u2 * v1 + f3 * (v1 * v2 + f))

        # computation of v/s (mod u)
        if s.divides(v3):
            v3 //= s
        else:
            s_inverse = s.inverse_mod(u3)
            v3 *= s_inverse

        v3 = v3 % u3

        return (u3, v3), (s.degree(), s.degree())

    def cantor_reduction(self, u0, v0):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Ensure D is a semi-reduced divisor
        assert u0.degree() >= g + 2, "Divisor has incorrect degree"
        assert (v0**2 + v0 * h - f) % u0 == 0, "D is not a valid divisor"

        # Compute u' and v'
        u1 = (v0**2 + h * v0 - f) // u0
        u1 = u1.monic()
        v1 = (-h - v0) % u1
        
        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        a_plus, a_minus = H.roots_at_infinity()

        leading_coefficient = v0[g+1] # check coefficient of x^(g+1)
        if leading_coefficient == a_plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        elif leading_coefficient == a_minus:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)
        else:
            omega = (d0 - d1) // 2
            omega_plus, omega_minus = (omega, omega)
        return (u1, v1), (omega_plus, omega_minus)


    def cantor_compose_at_infinity(self, u0, v0, plus=True):
        """
        Follows algorithm 3.6 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.infinite_polynomials()
        if plus:
            G = G_plus
        else:
            G = G_minus

        v1_prime = G + ((v0 - G) % u0)
        u1 = (v1_prime**2 + h * v1_prime - f) // u0
        u1 = u1.monic()
        v1 = (-h - v1_prime) % u1

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        if plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        else:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)

        return (u1, v1), (omega_plus, omega_minus)




class MumfordDivisorSplit():
    def __init__(self, parent, u, v, n=0, check=True):
        if not isinstance(parent, JacobianSplit):
            raise TypeError("parent must be of type ")
        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")
        #TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))
        # 2. define J(0) = (1,0:(g/2).ceil())

        self._parent = parent
        g = parent.curve().genus()

        # Ensure the divisor is valid
        if check:
            f, h = self._parent.curve().hyperelliptic_polynomials()
            assert (v**2 + v * h - f) % u == 0, f"{u = }, {v = } do not define a divisor on the Jacobian"

        self._u = u
        self._v = v

        assert 0 <= n <= (g - u.degree())
        self._n = n
        self._m = g - u.degree() - n

    def parent(self):
        return self._parent

    def __repr__(self) -> str:
        return f"({self._u}, {self._v} : {self._n})"

    def uv(self):
        return (self._u, self._v)

    def nm(self):
        return (self._n, self._m)

    def is_zero(self):
        g = self._parent.curve().genus()
        if self._n != (g/2).ceil():
            return False
        return self._u.is_one() and self._v.is_zero()

    def __eq__(self, other):
        if not isinstance(other, MumfordDivisorSplit):
            return False

        n1, _ = self.nm()
        n2, _ = other.nm()

        if n1 != n2:
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __hash__(self):
        data = (self._u, self._v, self._n)
        return hash(data)

    @cached_method
    def order(self):
        n = self.parent().order()
        return order_from_multiple(self, n)

    def degree(self):
        """
        Returns the degree of the affine part of the divisor.
        """
        return self._u.degree()

    def __add__(self, other):
        r"""
        Follows algorithm 3.7 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Ensure we are adding two divisors
        if not isinstance(other, MumfordDivisorSplit):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Extract out integers for weights
        n1, m1 = self.nm()
        n2, m2 = other.nm()

        # Compute weights for the sum
        omega_plus = n1 + n2
        omega_minus = m1 + m2

        # Step one: cantor composition of the two divisors
        (u3, v3), (a, b) = self._parent.cantor_composition(u1, v1, u2, v2)
        omega_plus += a
        omega_minus += b

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > (g + 1):
            (u3, v3), (a, b) = self._parent.cantor_reduction(u3, v3)
            omega_plus += a
            omega_minus += b

        # Step three: compose and then reduce at infinity to ensure
        # unique representation of D
        while (omega_plus < g/2) or (omega_minus < (g-1)/2):
            # TODO the paper has omega_plus > omega_minus
            # is this a typo? what happens for equal weight?
            # For example 3.2 after composition and reduction
            # the weights are (1, 1) and they reduce using + infty
            # which makes me thing this should be >= instead of >
            if omega_plus > omega_minus:
                (u3, v3), (a, b) = self._parent.cantor_compose_at_infinity(u3, v3)
            else:
                (u3, v3), (a, b) = self._parent.cantor_compose_at_infinity(u3, v3, plus=False)
            omega_plus += a
            omega_minus += b

        # Computing n3:
        # Algorithm states
        # E := D + omega^+ \infty^+ + omega^- \infty^- - D_\infty
        # Write E = D + n3 \infty^+ + m3 \infty^-
        # D_\infty = (g/2)(\infty^+ + \infty^-) when g is even
        #            (g+1)/2 \infty^+ + (g-1)/2 \infty_- when g is odd
        # So:
        # n3 = (omega^+ - (g/2).ceil())
        # m3 = (omega^- - (g/2).floor())

        n3 = omega_plus - (g/2).ceil()
        assert g - u3.degree() - n3 == omega_minus - (g/2).floor()

        return self._parent(u3, v3, n3)

    def __neg__(self):
        r"""
        Follows algorithm 3.8 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        _, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u1, v1 = self.uv()
        n1, m1 = self.nm()

        # Case for even genus
        if not (g % 2):
            v1 = (-h - v1) % u1
            n1 = m1
            return self.parent()(u1, v1, n1)
        # Odd genus, positive n1
        elif n1 > 0:
            v1 = (-h - v1) % u1
            #Note: I think this is a typo in the paper
            #n1 = g - m1 - u1.degree() + 1
            n1 = m1 + 1
            return self.parent()(u1, v1, n1)
        else:
            # Composition at infinity always with plus=True.
            # want to "substract" \infty_+ - \infty_-
            (u, v), (a, b) = self._parent.cantor_compose_at_infinity(u1, -h - v1, plus=True)
            return self.parent()(u, v, m1 + 1 + a)
        
    def __sub__(self, other):
        return self.__add__(-other)
    
    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, n):
        """
        """
        if not isinstance (n, (int, Integer)):
            raise ValueError

        if not n:
            return self._parent().zero()

        P = self

        # Handle negative scalars
        if n < 0:
            n = - n
            P = - P

        # Double and Add
        Q = P
        R = self.parent().zero()
        while n > 0:
            if n % 2 == 1:
                R = R + Q
            Q = Q + Q
            n = n // 2
        return R

    __rmul__ = __mul__
