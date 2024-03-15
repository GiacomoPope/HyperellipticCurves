import hyperelliptic_generic


class HyperellipticCurveSmoothModel_rational_field(
    hyperelliptic_generic.HyperellipticCurveSmoothModel_generic
):
    def __init__(self, projective_model, f, h, genus):
        super().__init__(projective_model, f, h, genus)

    # The below code is what is in sagemath right now, but it's a total mess??
    #
    # I dont think we should copy paste this but simply fix it?

    # def matrix_of_frobenius(self, p, prec=20):

    #     # BUG: should get this method from HyperellipticCurve_generic
    #     def my_chage_ring(self, R):
    #         from .constructor import HyperellipticCurve
    #         f, h = self._hyperelliptic_polynomials
    #         y = self._printing_ring.gen()
    #         x = self._printing_ring.base_ring().gen()
    #         return HyperellipticCurve(f.change_ring(R), h, "%s,%s" % (x,y))

    #     import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
    #     if isinstance(p, (sage.rings.abc.pAdicField, sage.rings.abc.pAdicRing)):
    #         K = p
    #     else:
    #         K = pAdicField(p, prec)
    #     frob_p, forms = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(my_chage_ring(self, K))
    #     return frob_p

    # def lseries(self, prec=53):
    #     """
    #     Return the L-series of this hyperelliptic curve of genus 2.

    #     EXAMPLES::

    #         sage: x = polygen(QQ, 'x')
    #         sage: C = HyperellipticCurve(x^2+x, x^3+x^2+1)
    #         sage: C.lseries()
    #         PARI L-function associated to Hyperelliptic Curve
    #         over Rational Field defined by y^2 + (x^3 + x^2 + 1)*y = x^2 + x
    #     """
    #     from sage.lfunctions.pari import LFunction, lfun_genus2
    #     L = LFunction(lfun_genus2(self), prec=prec)
    #     L.rename('PARI L-function associated to %s' % self)
    #     return L
