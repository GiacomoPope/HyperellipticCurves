from jacobian import HyperellipticJacobian, MumfordDivisor

class HyperellipticJacobianRamified(HyperellipticJacobian):
    def __init__(self, H):
        if not H.is_ramified():
            raise ValueError(f"the hyperelliptic curve {H} must have only one point at infinity")
        super(HyperellipticJacobianRamified, self).__init__(H)
        self._element = MumfordDivisorRamified

class MumfordDivisorRamified(MumfordDivisor):
    def __init__(self, parent, u, v, check=True):
        if not isinstance(parent, HyperellipticJacobianRamified):
            raise TypeError(f"parent must be of type {HyperellipticJacobianRamified}")

        super(MumfordDivisorRamified, self).__init__(parent, u, v, check=check)

