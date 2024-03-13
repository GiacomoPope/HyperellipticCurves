"""
Jacobian of a general hyperelliptic curve
"""

# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer import Integer
from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic

import jacobian_homset_ramified
import jacobian_homset_split
import jacobian_homset_inert

class HyperellipticJacobian_generic(Jacobian_generic):
    """
    TODO
    """
    def dimension(self):
        """
        Return the dimension of this Jacobian.
        """
        return Integer(self.curve().genus())

    # TODO why is check passed here and not used
    def point(self, mumford, check=True):
        try:
            return self(self.base_ring())(mumford)
        except AttributeError:
            raise ValueError("Arguments must determine a valid Mumford divisor.")

    def _point_homset(self, *args, **kwds):
        # TODO: make a constructor for this??
        H = self.curve()
        if H.is_ramified():
            return jacobian_homset_ramified.JacobianHomset_divisor_classes(*args, **kwds)
        elif  H.is_split():
            return jacobian_homset_split.JacobianHomset_divisor_classes(*args, **kwds)
        return jacobian_homset_inert.JacobianHomset_divisor_classes(*args, **kwds)

    def _point(self, *args, **kwds):
        return jacobian_morphism.JacobianMorphism_divisor_class_field(*args, **kwds)
