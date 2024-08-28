from jacobian_homset_generic import HyperellipticJacobianHomset
from jacobian_morphism import MumfordDivisorClassFieldRamified


class HyperellipticJacobianHomsetRamified(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldRamified
