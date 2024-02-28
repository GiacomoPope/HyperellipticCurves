from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from hyperelliptic_split import HyperellipticCurveSplit, HyperellipticPoint

class JacobianSplit:
    def __init__(self, H):
        if not isinstance(H, HyperellipticCurveSplit):
            raise ValueError("TODO")
        self._curve = H
        self._element = MumfordDivisorSplit

    def __repr__(self) -> str:
        return f"Jacobian of {self._curve}"

    def curve(self):
        return self._curve

    def __call__(self, *args):
        if isinstance(args[0], HyperellipticPoint):
            P = args[0]
            R, x = self._curve.polynomial_ring().objgen()
            X, Y = P.xy()
            u = x - X
            v = R(Y)
            n = 0 # FIXME
        elif len(args) == 2:
            u, v = args
            n = 0 # FIXME
        elif len(args) == 3:
            u, v, n = args
        else:
            raise NotImplementedError
        return self._element(self, u, v, n=n)


class MumfordDivisorSplit():
    def __init__(self, parent, u, v, n=0):
        if not isinstance(parent, JacobianSplit):
            raise TypeError(f"parent must be of type ")
        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")
        #TODO: 
        # 1. allow elements of the base field as input 
        #   (in particular something like (u,v) = (x-alpha, 0))
        # 2. define J(0) = (1,0:(g/2).ceil())

        self._parent = parent
        g = parent.curve().genus()

        self._u = u
        self._v = v

        # assert 0 <= n <= (g - u.degree())
        self._n = n

        # FIXME: is having n=0 by default OK? I have no idea!
        self._m = g - u.degree() - n

    def parent(self):
        return self._parent

    def __repr__(self) -> str:
        return f"({self._u}, {self._v} : {self._n})"

    def uv(self):
        return (self._u, self._v)

    def nm(self):
        return (self._n, self._m)

    def degree(self):
        """
        TODO: is this correct?
        """
        return self._u.degree()

    def cantor_composition(self, other):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Ensure we are composing two divisors
        if not isinstance(other, MumfordDivisorSplit):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Ensure D1 and D2 are semi-reduced divisors
        assert (
            v1.degree() < u1.degree() and v2.degree() < u2.degree()
        ), "The degree of bi must be smaller than ai"
        assert (
            u1.degree() <= f.degree() and u2.degree() <= f.degree()
        ), f"The degree of ai must be smaller than f, {u1.degree()}, {u2.degree()}"
        assert (v1**2 + v1 * h - f) % u1 == 0, "D1 is not a valid divisor"
        assert (v2**2 + v2 * h - f) % u2 == 0, "D2 is not a valid divisor"

        # Step One
        s0, e1, e2 = u1.xgcd(u2)
        assert s0 == e1 * u1 + e2 * u2

        # Step Two
        s, c1, c2 = s0.xgcd(v1 + v2 + h)
        assert s == c1 * s0 + c2 * (v1 + v2 + h)

        # Step Three
        f1 = c1 * e1
        f2 = c1 * e2
        f3 = c2
        assert s == f1 * u1 + f2 * u2 + f3 * (v1 + v2 + h)

        # Step Four
        u = (u1 * u2) // (s**2)
        v = (f1 * u1 * v2 + f2 * u2 * v1 + f3 * (v1 * v2 + f))

        # computation of v/s (mod u)
        c, d_inverse, _ = s.xgcd(u)
        assert c.divides(v)
        v = (v // c * d_inverse) % u
        D3 = self.parent()(u, v)

        return D3, (s.degree(), s.degree())

    def cantor_reduction(self):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Mumford polynomials
        u0, v0 = self.uv()

        # Ensure D is a semi-reduced divisor
        # assert u0.degree() >= g + 2, "Divisor has incorrect degree"
        assert (v0**2 + v0 * h - f) % u0 == 0, "D is not a valid divisor"

        # Compute u' and v'
        u1 = (f - v0 * h - v0**2) // u0
        u1 = u1.monic()
        v1 = (-h - v0) % u1
        D1 = self.parent()(u1, v1)

        # Compute the counter weights
        d0 = self.degree()
        d1 = D1.degree()
        a_plus, a_minus = self.parent().curve().roots_at_infinity()

        if v0.leading_coefficient() == a_plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        elif v0.leading_coefficient() == a_minus:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)
        else:
            omega = (d0 - d1) // 2
            omega_plus, omega_minus = (omega, omega)
        return D1, (omega_plus, omega_minus)


    def cantor_compose_at_infinity(self, plus=True):
        """
        Follows algorithm 3.6 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.infinite_polynomials()
        if plus:
            G = G_plus
        else:
            G = G_minus

        u0, v0 = self.uv()
        v1_prime = G + (v0 - G) % u0

        u1 = (v1_prime**2 + h * v1_prime - f) // u0
        u1 = u1.monic()
        v1 = (-h - v1_prime) % u1
        D1 = self.parent()(u1, v1)

        # Compute the counter weights
        d0 = self.degree()
        d1 = D1.degree()
        if plus:
            omega_plus, omega_minus = (d0 - g - 1, g + 1 - d1)
        else:
            omega_plus, omega_minus = (g + 1 - d1, d0 - g - 1)

        return D1, (omega_plus, omega_minus)

    def __add__(self, other):
        r"""
        TODO: this is not going to work until I really understand
        the relationship between (u,v) with n and m in the representation
        div(u,v) + n*\infty^+ + m*\infty^- - F_\infty

        Follows algorithm 3.7 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Ensure we are adding two divisors
        if not isinstance(other, MumfordDivisorSplit):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        _, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Extract out Mumford polynomials
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Extract out integers for weights
        n1, m1 = self.nm()
        n2, m2 = other.nm()

        omega_plus = n1 + n2
        omega_minus = m1 + m2

        # Step one: cantor composition of the two divisors
        D, (a, b) = self.cantor_composition(other)
        omega_plus += a
        omega_minus += b

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while D.degree() > (g + 1):
            D, (a, b) = D.cantor_reduction()
            omega_plus += a
            omega_minus += b

        # Step three: compose and then reduce at infinity to ensure
        # unique representation of D
        while (omega_plus < g/2) or (omega_minus < (g-1)/2):
            # TODO the paper has omega_plus > omega_minus
            # is this a typo? what happens for equal weight?
            # For example 3.2 after composition and reduction
            # the weights are (1, 1) and they reduce using + infty
            # which makes me thing this should be >= instead of >
            if omega_plus > omega_minus:
                D, (a, b) = D.cantor_compose_at_infinity()
            else:
                D, (a, b) = D.cantor_compose_at_infinity(plus=False)
            omega_plus += a
            omega_minus += b

        # TODO:
        # How to do compute n3 from this?
        # Algorithm states
        # E := D + omega^+ \infty^+ + omega^- \infty^- - D_\infty
        # Write E = D + n3 \infty^+ + m3 \infty^-
        # I believe D_\infty = (g/2)(\infty^+ + \infty^-) when g is even
        #                      ((g+1)/2)(\infty^+ + \infty^-) when g is odd
        # So I think we want:
        # n3 = (omega^+ - (g/2).ceil())
        # m3 = (omega^- - (g/2).ceil())

        # EDIT:
        # In the odd genus case, D_\infty = (g+1)/2 \infty_+ + (g-1)/2 \infty_-
        # so m3 = (omega^- - (g/2).floor())
        u3, v3 = D.uv()
        n3 = omega_plus - (g/2).ceil() # FIX ME?
        m3 = omega_minus - (g/2).floor()
        return self._parent(u3, v3, n3)

    def __neg__(self):
        r"""
        TODO: this is not going to work until I really understand
        the relationship between (u,v) with n and m in the representation
        div(u,v) + n*\infty^+ + m*\infty^- - D_\infty

        Follows algorithm 3.8 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        _, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u1, v1 = self.uv()
        n1, m1 = self.nm()

        # Case for even genus
        if not (g % 2):
            v1 = (-h - v1) % u1
            n1 = m1
            return self.parent()(u1, v1, n1)
        # Odd genus, positive n1
        elif n1 > 0:
            v1 = (-h - v1) % u1
            #Note: I think this is a typo in the paper 
            #n1 = g - m1 - u1.degree() + 1 
            n1 = m1 + 1
            return self.parent()(u1, v1, n1)
        else:
            # Composition at infinity always with plus=True. 
            # want to "substract" \infty_+ - \infty_-
            D = self.parent()(u1, -h-v1)
            D, (a, b) = D.cantor_compose_at_infinity(plus=True)
            u, v = D.uv()
            D = self.parent()(u, v, m1 + 1 + a)
            return D