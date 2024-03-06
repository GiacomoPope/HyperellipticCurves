from sage.rings.polynomial.polynomial_element import Polynomial
from hyperelliptic import HyperellipticCurveNew
from sage.rings.integer import Integer

from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method


class HyperellipticJacobian:
    def __init__(self, H):
        if not isinstance(H, HyperellipticCurveNew):
            raise ValueError("TODO")
        self._curve = H
        self._element = MumfordDivisor

    def __repr__(self) -> str:
        return f"Jacobian of {self._curve}"

    def curve(self):
        return self._curve

    def base_ring(self):
        return self._curve.base_ring()

    @cached_method
    def cardinality(self):
        """
        TODO: currently using lazy methods by calling sage
        """
        return sum(self.curve().frobenius_polynomial())

    order = cardinality

    def __cantor_double_generic(self, u1, v1):
        """
        Efficient cantor composition for doubling an affine divisor

        Returns the Cantor composition of (u1, v1) with (u1, v1)
        together with the degree of the polynomial ``s`` which is
        needed for computing weights for the split and inert models 
        """
        f, h = self._curve.hyperelliptic_polynomials()

        # New mumford coordinates
        if h.is_zero():
            s, _, e2 = u1.xgcd(v1 + v1)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1**2) // s
        else:
            s, _, e2 = u1.xgcd(v1 + v1 + h)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1 * h - v1**2) // s
        v3 = v3 % u3

        return u3, v3, s.degree()

    def _cantor_composition_generic(self, u1, v1, u2, v2):
        """
        Cantor composition

        Returns the Cantor composition of (u1, v1) with (u2, v2)
        together with the degree of the polynomial ``s`` which is
        needed for computing weights for the split and inert models 
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

        # Special case: duplication law
        if u1 == u2 and v1 == v2:
            return self.__cantor_double_generic(u1, v1)

        # Step One
        s0, _, e2 = u1.xgcd(u2)

        # Special case: when gcd(u0, u1) == 1 we can
        # avoid many expensive steps as we have s = 1
        if s0.is_one():
            u3 = u1 * u2
            v3 = v2 + e2 * u2 * (v1 - v2)
            v3 = v3 % u3
            return u3, v3, 0

        # Step Two
        w0 = v1 + v2 + h

        # Another special case, when w0 is zero we skip
        # a xgcd and can return early
        if w0.is_zero():
            u3 = (u1 * u2) // (s0**2)
            v3 = v2 + e2 * (v1 - v2) * (u2 // s0)
            v3 = v3 % u3
            return u3, v3, s0.degree()

        # Step Three
        s, c1, c2 = s0.xgcd(w0)
        u3 = (u1 * u2) // (s**2)
        v3 = v2 + (c1 * e2 * (v1 - v2) * u2 + c2 * (f - h * v2 - v2**2)) // s
        v3 = v3 % u3
        return u3, v3, s.degree()

    def _cantor_reduction_generic(self, u0, v0):
        """
        TODO
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

        return u1, v1

class MumfordDivisor:
    def __init__(self, parent, u, v, check=True):
        if not isinstance(parent, HyperellipticJacobian):
            raise TypeError(f"parent must be of type {HyperellipticJacobian}")
        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")
 
        # TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))

        # Ensure the divisor is valid
        if check:
            f, h = parent.curve().hyperelliptic_polynomials()
            assert (
                v**2 + v * h - f
            ) % u == 0, f"{u = }, {v = } do not define a divisor on the Jacobian"

        self._parent = parent
        self._u = u
        self._v = v

    def parent(self):
        return self._parent

    def __repr__(self) -> str:
        return f"({self._u}, {self._v})"

    def uv(self):
        return (self._u, self._v)

    def __eq__(self, other):
        if not isinstance(other, MumfordDivisor):
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __hash__(self):
        data = (self._u, self._v)
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

    def __add__(self, *args):
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, n):
        """ """
        # TODO: is there a better handlings for this?
        if not isinstance(n, (int, Integer)):
            raise ValueError

        if not n:
            return self._parent().zero()

        P = self

        # Handle negative scalars
        if n < 0:
            n = -n
            P = -P

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
