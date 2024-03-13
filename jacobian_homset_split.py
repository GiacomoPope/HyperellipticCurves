from jacobian_homset_generic import HyperellipticJacobianHomset
from jacobian_morphism import MumfordDivisorClassFieldSplit

# TODO should we make a hyperelliptic point class?
# at the moment, this is the type we get from calling a point from the projective model
from sage.schemes.toric.morphism import SchemeMorphism_point_toric_field


class HyperellipticJacobianHomsetSplit(HyperellipticJacobianHomset):
    def __init__(self, Y, X):
        super(HyperellipticJacobianHomsetSplit, self).__init__(Y, X)
        self._element = MumfordDivisorClassFieldSplit

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self._curve.genus()
        R = self._curve.polynomial_ring()
        n = (g / 2).ceil()
        return self._element(self, R.one(), R.zero(), n)

    def point_to_mumford_coordinates(self, P):
        """
        TODO
        """
        R, x = self._curve.polynomial_ring().objgen()
        g = self._curve.genus()

        [X, Y, Z] = P._coords
        # we use the embedding P \mapsto P - P0
        # where P0 is the distinguished point of the curve
        # TODO: at the moment, we assume P0 = infty_+
        # note: P - \inft_+ = P + n*\infty_+ + m*\infty_- - D_\infty,
        # where n = ((g-1)/2).floor()
        n = ((g - 1) / 2).floor()
        if Z == 0:
            alpha = Y / X
            if alpha == self._curve._alphas[0]:
                n = n + 1
            u = R.one()
            v = R.zero()
        else:
            u = x - X
            v = R(Y)

        return u, v, n

    def __call__(self, *args, check=True):
        """
        TODO: fix code reuse?
        """
        if isinstance(args[0], SchemeMorphism_point_toric_field) and len(args) == 1:
            u, v, n = self.point_to_mumford_coordinates(args[0])
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
        return self._element(self, u, v, n=n, check=check)

    def cantor_composition(self, u1, v1, n1, u2, v2, n2):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        TODO: when h = 0 we can speed this up.
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Cantor composition
        u3, v3, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)

        # Compute new weight
        n3 = n1 + n2 + s_deg - (g / 2).ceil()

        return u3, v3, n3

    def cantor_reduction(self, u0, v0, n0):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Perform regular cantor reduction
        u1, v1 = self._cantor_reduction_generic(u0, v0)

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        a_plus, a_minus = H.roots_at_infinity()

        if v0.degree() <= g + 1:
            leading_coefficient = v0[g + 1]  # check coefficient of x^(g+1)
            if leading_coefficient == a_plus:
                n1 = n0 + d0 - g - 1
            elif leading_coefficient == a_minus:
                n1 = n0 + g + 1 - d1
            else:
                n1 = n0 + (d0 - d1) // 2
        else:
            n1 = n0 + (d0 - d1) // 2
        return u1, v1, n1

    def cantor_compose_at_infinity(self, u0, v0, n0, plus=True):
        """
        Follows algorithm 3.6 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.infinite_polynomials()
        if plus:
            G = G_plus
        else:
            G = G_minus

        v1_prime = G + ((v0 - G) % u0)
        u1 = (v1_prime**2 + h * v1_prime - f) // u0
        u1 = u1.monic()
        v1 = (-h - v1_prime) % u1

        # Compute the counter weights
        if plus:
            n1 = n0 + u0.degree() - g - 1
        else:
            n1 = n0 + g + 1 - u1.degree()

        return u1, v1, n1
