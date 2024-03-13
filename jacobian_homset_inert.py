from jacobian_homset_generic import HyperellipticJacobianHomset
from jacobian_morphism import MumfordDivisorClassFieldInert

class HyperellipticJacobianHomsetInert(HyperellipticJacobianHomset):
    def __init__(self, Y, X):
        super(HyperellipticJacobianHomsetInert, self).__init__(Y, X)
        self._morphism_element = MumfordDivisorClassFieldInert

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self._curve.genus()
        if g % 2:
            raise ValueError(
                "unable to perform arithmetic for inert models of odd genus"
            )
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
        [X, Y, Z] = P._coords
        [X0, Y0, Z0] = P0._coords
        assert Z != 0 and Z0 != 0, "there should be no points at infinity"
        _, h = self._curve.hyperelliptic_polynomials()
        u, v, _ = self._cantor_composition_generic(x - X, R(Y), x - X0, R(-Y0 - h(X0)))
        return u, v
