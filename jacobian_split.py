from sage.misc.prandom import randint
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.misc.prandom import choice

# Needed until https://github.com/sagemath/sage/pull/37118 is merged.
from uniform_random_sampling import uniform_random_polynomial

# Base classes
from hyperelliptic import HyperellipticPoint
from jacobian import HyperellipticJacobian, MumfordDivisor

class HyperellipticJacobianSplit(HyperellipticJacobian):
    def __init__(self, H):
        super(HyperellipticJacobianSplit, self).__init__(H)
        self._element = MumfordDivisorSplit

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self._curve.genus()
        R = self._curve.polynomial_ring()
        n = (g / 2).ceil()
        return self._element(self, R.one(), R.zero(), n)

    def __call__(self, *args, check=True):
        if isinstance(args[0], HyperellipticPoint):
            R, x = self._curve.polynomial_ring().objgen()
            g = self._curve.genus()
            P = args[0]
            [X, Y, Z] = P.coords()
            # TODO:
            # we use the embedding P \mapsto P - \infty_+
            # note: P - \inft_+ = P + n*\infty_+ + m*\infty_- - D_\infty,
            # where n = ((g-1)/2).floor()
            # do we want a different embedding by default?
            n = ((g - 1) / 2).floor()
            if Z == 0:
                alpha = Y / X
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
        return self._element(self, u, v, n=n, check=check)

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
                return self(R(1), R(0), n)
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
                    except (ValueError, AttributeError):
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
        points = [H.random_point() for _ in range(2 * g + 1)]

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
            raise NotImplementedError(
                "random element of Jacobian is only implemented over Finite Fields"
            )

        if fast:
            return self._random_element_rational(*args, **kwargs)
        return self._random_element_cover(*args, **kwargs)

    def cantor_composition(self, u1, v1, n1, u2, v2, n2):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        TODO: when h = 0 we can speed this up.
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Cantor composition
        u3, v3, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)
        
        # Compute new weight
        n3 = n1 + n2 + s_deg - (g / 2).ceil()

        return u3, v3, n3

    def cantor_reduction(self, u0, v0, n0):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Perform regular cantor reduction
        u1, v1 = self._cantor_reduction_generic(u0, v0)

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        a_plus, a_minus = H.roots_at_infinity()

        if v0.degree() == g + 1:
            leading_coefficient = v0[g + 1]  # check coefficient of x^(g+1)
            if leading_coefficient == a_plus:
                n1 = n0 + d0 - g - 1
            elif leading_coefficient == a_minus:
                n1 = n0 + g + 1 - d1
            else:
                n1 = n0 + (d0 - d1) // 2
        else:
            n1 = n0 + (d0 - d1) // 2
        return u1, v1, n1

    def cantor_compose_at_infinity(self, u0, v0, n0, plus=True):
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
        if plus:
            n1 = n0 + u0.degree() - g - 1
        else:
            n1 = n0 + g + 1 - u1.degree()

        return u1, v1, n1


class MumfordDivisorSplit(MumfordDivisor):
    def __init__(self, parent, u, v, n=0, check=True):
        if not isinstance(parent, HyperellipticJacobianSplit):
            raise TypeError(f"parent must be of type {HyperellipticJacobianSplit}")

        # Ensure the weight is set correctly
        g = parent.curve().genus()
        assert 0 <= n <= (g - u.degree())
        self._n = n
        self._m = g - u.degree() - n

        super(MumfordDivisorSplit, self).__init__(parent, u, v, check=check)

    def __repr__(self):
        return f"({self._u}, {self._v} : {self._n})"

    def is_zero(self):
        g = self._parent.curve().genus()
        if self._n != (g / 2).ceil():
            return False
        return self._u.is_one() and self._v.is_zero()

    def __eq__(self, other):
        if not isinstance(other, MumfordDivisorSplit):
            return False

        n1, n2 = self._n, other._n

        if n1 != n2:
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __hash__(self):
        data = (self._u, self._v, self._n)
        return hash(data)

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
        n1, n2 = self._n, other._n

        # Step one: cantor composition of the two divisors
        u3, v3, n3 = self._parent.cantor_composition(u1, v1, n1, u2, v2, n2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > (g + 1):
            u3, v3, n3 = self._parent.cantor_reduction(u3, v3, n3)
        u3 = u3.monic()

        # Step three: compose and then reduce at infinity to ensure
        # unique representation of D
        while n3 < 0 or n3 > g - u3.degree():
            if n3 < 0:
                u3, v3, n3 = self._parent.cantor_compose_at_infinity(
                    u3, v3, n3, plus=False
                )
            else:
                u3, v3, n3 = self._parent.cantor_compose_at_infinity(
                    u3, v3, n3, plus=True
                )

        return self._parent(u3, v3, n3, check=False)

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

        u0, v0 = self.uv()
        n0, m0 = self._n, self._m

        # Case for even genus
        if not (g % 2):
            v1 = (-h - v0) % u0
            u1 = u0
            n1 = m0
        # Odd genus, positive n0
        elif n0 > 0:
            v1 = (-h - v0) % u0
            # Note: I think this is a typo in the paper
            # n1 = g - m0 - u1.degree() + 1
            u1 = u0
            n1 = m0 + 1
        else:
            # Composition at infinity always with plus=True.
            # want to "substract" \infty_+ - \infty_-
            (u1, v1, n1) = self._parent.cantor_compose_at_infinity(
                u0, -h - v0, n0, plus=True
            )
            n1 = n1 - n0 + m0 + 1
        
        return self.parent()(u1, v1, n1, check=False)
