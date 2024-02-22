from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from hyperelliptic_split import HyperellipticCurveSplit, HyperellipticPoint

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

    def __call__(self, *args):
        if isinstance(args[0], HyperellipticPoint):
            P = args[0]
            R, x = self._curve.polynomial_ring().objgen()
            X, Y = P.xy()
            u = x - X
            v = R(Y)
        elif len(args) == 2:
            u, v = args
        else:
            raise NotImplementedError
        return self._element(self, u, v)


class MumfordDivisorSplit():
    def __init__(self, parent, u, v):
        if not isinstance(parent, JacobianSplit):
            raise TypeError(f"parent must be of type ")
        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")

        self._parent = parent
        self._u = u
        self._v = v

    def parent(self):
        return self._parent

    def __repr__(self) -> str:
        return f"({self._u}, {self._v})"

    def uv(self):
        return (self._u, self._v)

    def cantor_composition(self, other):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Ensure we are composing two divisors
        if not isinstance(other, MumfordDivisorSplit):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Ensure D1 and D2 are semi-reduced divisors
        assert (
            v1.degree() < u1.degree() and v2.degree() < u2.degree()
        ), "The degree of bi must be smaller than ai"
        assert (
            u1.degree() <= f.degree() and u2.degree() <= f.degree()
        ), f"The degree of ai must be smaller than f, {u1.degree()}, {u2.degree()}"
        assert (v1**2 + v1 * h - f) % u1 == 0, "D1 is not a valid divisor"
        assert (v2**2 + v2 * h - f) % u2 == 0, "D2 is not a valid divisor"

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
        u = (u1 * u2) // (s**2)
        v = (f1 * u1 * v2 + f2 * u2 * v1 + f3 * (v1 * v2 + f))

        # TODO:
        # Sometimes I cannot invert d mod a, but it seems like in these cases
        # I can always just divide it out (at least experimentally). Potentially
        # edge cases...
        if s.divides(v):
            v = (v // s) % u
        else:
            d_inverse = s.inverse_mod(u)
            v = (v * d_inverse) % u
        D3 = self.parent()(u, v)
        return D3, (s.degree(), s.degree())

    def cantor_reduction(self):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Mumford polynomials
        u0, v0 = self.uv()

        # Ensure D is a semi-reduced divisor
        assert u0.degree() >= g + 2, "Divisor has incorrect degree"
        assert (v0**2 + v0 * h - f) % u0 == 0, "D is not a valid divisor"

        # Compute u' and v'
        u1 = (f - v0 * h - v0**2) // u0
        u1 = u1.monic()
        v1 = (-h - v0) % u1
        D1 = self.parent()(u1, v1)

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        a_plus, a_minus = self.parent().curve().roots_at_infinity()

        if v0.leading_coefficient() == a_plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        elif v0.leading_coefficient() == a_minus:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)
        else:
            omega = (d0 - d1) // 2
            omega_plus, omega_minus = (omega, omega)
        return D1, (omega_plus, omega_minus)


    def cantor_reduction_at_infinity(self, plus=True):
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.infinite_polynomials()
        if plus:
            G = G_plus
        else:
            G = G_minus

        u0, v0 = self.uv()
        v1_prime = G + (v0 - G) % u0

        u1 = (v1_prime**2 + h * v1_prime - f) // u0
        u1 = u1.monic()
        v1 = (-h - v1_prime) % u1
        D1 = self.parent()(u1, v1)

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        if plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        else:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)

        return D1, (omega_plus, omega_minus)
