from sage.misc.prandom import choice
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

import hyperelliptic_generic

class HyperellipticCurveSmoothModel_finite_field(hyperelliptic_generic.HyperellipticCurveSmoothModel_generic):
    def __init__(self, projective_model, f, h, genus):
        super().__init__(projective_model, f, h, genus)

    def random_point(self):
        """
        Return a random point on this hyperelliptic curve, uniformly chosen
        among all rational points.
        """
        k = self.base_ring()
        n = 2 * k.order() + 1

        while True:
            # Choose the point at infinity with probability 1/(2q + 1)
            i = ZZ.random_element(n)
            if not i and not self.is_inert():
                # TODO: deal with that there is more than one point at infinity
                return choice(self.points_at_infinity())
            v = self.lift_x(k.random_element(), all=True)
            try:
                return v[i % 2]
            except IndexError:
                pass

    def frobenius_polynomial(self):
        """
        TODO: very lazy but just for now
        """
        f, h = self._hyperelliptic_polynomials
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve

        H_tmp = HyperellipticCurve(f, h)
        return H_tmp.frobenius_polynomial()

    @cached_method
    def points(self):
        """
        TODO: couldn't be more stupid
        """
        # TODO: this is very silly but works
        points = self.points_at_infinity()
        for x in self.base_ring():
            points.extend(self.lift_x(x, all=True))
        return points

    @cached_method
    def cardinality(self):
        """
        TODO: couldn't be more stupid
        """
        return len(self.points())

    order = cardinality
