from jacobian_generic import HyperellipticJacobian_generic

import jacobian_g2_homset_ramified
import jacobian_g2_homset_split
import jacobian_g2_homset_inert


class HyperellipticJacobian_g2_generic(HyperellipticJacobian_generic):
    """
    Special class to handle optimisations for jacobian computations
    in genus two
    """

    def _point_homset(self, *args, **kwds):
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_g2_homset_ramified.HyperellipticJacobianHomsetRamified_g2(
                *args, **kwds
            )
        elif H.is_split():
            return jacobian_g2_homset_split.HyperellipticJacobianHomsetSplit_g2(
                *args, **kwds
            )
        return jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2(
            *args, **kwds
        )
