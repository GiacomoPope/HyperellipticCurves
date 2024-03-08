from hyperelliptic import HyperellipticPoint
from jacobian import HyperellipticJacobian, MumfordDivisor

class HyperellipticJacobianInert(HyperellipticJacobian):
    def __init__(self, H):
        if not H.is_inert():
            raise ValueError(f"the hyperelliptic curve {H} must have only no rational points at infinity")
        super(HyperellipticJacobianInert, self).__init__(H)
        self._element = MumfordDivisorInert

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self._curve.genus()
        if (g % 2):
            raise ValueError("unable to perform arithmetic for inert models of odd genus")
        R = self._curve.polynomial_ring()
        return self._element(self, R.one(), R.zero(), g // 2)
    
    def point_to_mumford_coordinates(self, P):
        """
        TODO
        """
        R, x = self._curve.polynomial_ring().objgen()
        # we use the embedding P \mapsto P - P0
        # where P0 is the distinguished point of the curve
        # note: P - \inft_+ = P + n*\infty_+ + m*\infty_- - D_\infty,
        # where n = ((g-1)/2).floor()
        P0 = self._curve.distinguished_point()
        if P == P0:
            return R.one(), R.zero()
        [X, Y, Z] = P.coords()
        [X0, Y0, Z0] = P0.coords()
        assert Z != 0 and Z0 != 0, "there should be no points at infinity"
        f, h = self._curve.hyperelliptic_polynomials()
        u, v, _ = self._cantor_composition_generic(x-X, R(Y), x-X0, R(-Y0-h(X0)))
        return u, v

class MumfordDivisorInert(MumfordDivisor):
    def __init__(self, parent, u, v, check=True):
        if not isinstance(parent, HyperellipticJacobianInert):
            raise TypeError(f"parent must be of type {HyperellipticJacobianInert}")

        if (u.degree() % 2):
            raise ValueError(f"mumford coordinate {u} must have even degree")
        
        g = parent.curve().genus()
        self._n = (g - u.degree()) // 2
        super(MumfordDivisorInert, self).__init__(parent, u, v, check=check)

    def __repr__(self):
        return f"({self._u}, {self._v} : {self._n})"

    def __hash__(self):
        data = (self._u, self._v, self._n)
        return hash(data)
