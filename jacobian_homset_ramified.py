from jacobian_homset_generic import HyperellipticJacobianHomset
from jacobian_morphism import MumfordDivisorClassFieldRamified

class HyperellipticJacobianHomsetRamified(HyperellipticJacobianHomset):
    def __init__(self, Y, X):
        super(HyperellipticJacobianHomsetRamified, self).__init__(Y, X)
        self._element = MumfordDivisorClassFieldRamified
