

# Base classes
from hyperelliptic import HyperellipticPoint
from jacobian import HyperellipticJacobian, MumfordDivisor

class HyperellipticJacobianRamified(HyperellipticJacobian):
    def __init__(self, H):
        if not H.is_ramified():
            raise ValueError(f"the hyperelliptic curve {H} must have only one point at infinity")
        super(HyperellipticJacobianRamified, self).__init__(H)
        self._element = MumfordDivisorRamified

    def __call__(self, *args, check=True):
        if isinstance(args[0], HyperellipticPoint):
            R, x = self._curve.polynomial_ring().objgen()
            P = args[0]
            [X, Y, Z] = P.coords()
            if Z == 0:
                return self.zero()
            else:
                u = x - X
                v = R(Y)
        # TODO handle this better!!
        elif len(args) == 1:
            return self.zero()
        elif len(args) == 2:
            u, v = args
            n = 0
        elif len(args) == 3:
            u, v, n = args
        else:
            raise NotImplementedError
        return self._element(self, u, v, check=check)

    def cantor_composition(self, u1, v1, u2, v2):
        u3, v3, _ = self._cantor_composition_generic(u1, v1, u2, v2)
        return u3, v3

    def cantor_reduction(self, u0, v0):
        return self._cantor_reduction_generic(u0, v0)

class MumfordDivisorRamified(MumfordDivisor):
    def __init__(self, parent, u, v, check=True):
        if not isinstance(parent, HyperellipticJacobianRamified):
            raise TypeError(f"parent must be of type {HyperellipticJacobianRamified}")

        super(MumfordDivisorRamified, self).__init__(parent, u, v, check=check)

    def is_zero(self):
        return self._u.is_one() and self._v.is_zero()

    def __add__(self, other):
        # Ensure we are adding two divisors
        if not isinstance(other, MumfordDivisorRamified):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Step one: cantor composition of the two divisors
        u3, v3 = self._parent.cantor_composition(u1, v1, u2, v2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > g:
            u3, v3 = self._parent.cantor_reduction(u3, v3)
        u3 = u3.monic()
        v3 = v3 % u3

        return MumfordDivisorRamified(self._parent, u3, v3, check=False)

    def __neg__(self):
        _, h = self.parent().curve().hyperelliptic_polynomials()
        u0, v0 = self.uv()

        if h.is_zero():
            return MumfordDivisorRamified(self._parent, u0, -v0, check=False)

        v1 = (-h - v0) % u0
        return MumfordDivisorRamified(self._parent, u0, v1, check=False)
