from sage.misc.prandom import choice
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

from sage.schemes.toric.library import toric_varieties
from sage.schemes.toric.toric_subscheme import AlgebraicScheme_subscheme_toric

# Named HyperellipticCurveSmoothModel to differentiate from the global 
# HyperellipticCurve in Sage
class HyperellipticCurveSmoothModel(AlgebraicScheme_subscheme_toric):
    def __init__(self, f, h=0):
        # Some values which we will cache as a user asks for them
        self._alphas = None
        self._genus = None
        self._infinte_polynomials = None
        self._distinguished_point = None

        # Check the curve is well formed
        disc = h**2 + 4*f
        if not isinstance(disc, Polynomial):
            raise TypeError(f"arguments {f = } and {h = } must be polynomials")
        self._polynomial_ring = disc.parent()
        self._base_ring = disc.base_ring()

        # TODO: check singularity from disc

        # Store the hyperelliptic polynomials as the correct type
        f = self._polynomial_ring(f)
        h = self._polynomial_ring(h)
        self._hyperelliptic_polynomials = (f, h)
        self._d = max(h.degree(), (f.degree() / 2).ceil())

        # Compute the smooth model for the hyperelliptic curve
        # using a weighted projective space (via Toric Variety)
        self._projective_model = self.projective_model()

    def projective_model(self):
        """
        Compute the weighted projective model (1 : g + 1 : 1)
        """
        if hasattr(self, "_projective_model"):
            return self._projective_model

        f, h = self._hyperelliptic_polynomials
        g = self.genus()

        T = toric_varieties.WP([1, g + 1, 1], base_ring=self._base_ring, names="X, Y, Z")
        (X, Y, Z) = T.gens()

        d = max(h.degree(), (f.degree() / 2).ceil())
        F = sum(f[i] * X**i * Z**(2*d - i) for i in range(2*d + 1))

        if h.is_zero():
            G = Y**2 - F
        else:
            H = sum(h[i] * X**i * Z**(d - i) for i in range(d + 1))
            G = Y**2 + H*Y - F
        
        self._projective_model = T.subscheme([G])
        return self._projective_model

    def __repr__(self):
        f, h = self._hyperelliptic_polynomials
        if h:
            if h.is_one():
                curve = f"y^2 + y = {f}"
            else:
                curve = f"y^2 + ({h})*y = {f}"
        else:
            curve = f"y^2 = {f}"
        return f"Hyperelliptic Curve over {self.base_ring()} defined by {curve}"

    def genus(self):
        """
        Compute the genus of the hyperelliptic curve
        """
        if self._genus:
            return self._genus

        f, h = self._hyperelliptic_polynomials
        df = f.degree()
        dh_2 = 2*h.degree()
        if dh_2 < df:
            self._genus = (df-1)//2
        else:
            self._genus = (dh_2-1)//2

        return self._genus

    def base_ring(self):
        return self._base_ring

    def polynomial_ring(self):
        return self._polynomial_ring

    def point(self, coords, check=True):
        if len(coords) == 2:
            X, Y = coords
            Z = self.base_ring().one()
        elif len(coords) == 3:
            X, Y, Z = coords
        else:
            raise ValueError("TODO")

        return self._projective_model.point([X, Y, Z], check=check)

    def hyperelliptic_polynomials(self):
        """
        Return the polynomials (f, h) such that
        C : y^2 + h*y = f
        """
        return self._hyperelliptic_polynomials

    def roots_at_infinity(self):
        """
        Compute the roots of: Y^2 + h[d]Y - d[2d] = 0
        When the curve is ramified, we expect one root, when
        the curve is split or inert we expect zero or two roots.
        """
        if self._alphas:
            return self._alphas

        f, h = self._hyperelliptic_polynomials
        x = f.parent().gen()
        d = self._d
        
        if not h:
            coeff = f[2*d]
            # Handle the ramified case
            if coeff.is_zero():
                return [coeff]
            return f[2*d].sqrt(all=True)

        self._alphas = (x**2 + x * h[d] - f[2*d]).roots(multiplicities=False)
        return self._alphas
    
    def is_split(self):
        """
         Return True if the curve is split, i.e. there are two rational
         points at infinity.
         """
        return len(self.roots_at_infinity()) == 2
    
    def is_ramified(self):
        """
        Return True if the curve is ramified, i.e. there is one rational
        point at infinity.
        """
        return len(self.roots_at_infinity()) == 1

    def is_inert(self):
        """
        Return True if the curve is inert, i.e. there are no rational
        points at infinity.
        """
        return len(self.roots_at_infinity()) == 0

    def infinite_polynomials(self):
        """
        TODO: stupid name

        Computes G^±(x) for curves in the split degree model
        """
        if self._infinte_polynomials is not None:
            return self._infinte_polynomials

        alphas = self.roots_at_infinity()

        # This function only makes sense for the split model
        if not len(alphas) == 2:
            raise ValueError("hyperelliptic curve does not have the split model")

        f, h = self._hyperelliptic_polynomials
        alpha_plus, alpha_minus = alphas
        d = self._d

        # Construct G_plus from alpha_plus
        g = [None] * (d + 1)
        g[d] = alpha_plus
        for i in range(d - 1, -1, -1):
            # You just want (g * (g + h))[x^(i + d)] to match f_{i + d}
            the_rest = g[d] * h[i] + sum(g[k] * (g[i + d - k] + h[i + d - k]) for k in range(i + 1, d))
            g[i] = (f[i + d] - the_rest) / (2 * g[d] + h[d])

        G_plus = self._polynomial_ring(g)
        G_minus = - G_plus - h
        # Checks for the assumptions on G^±
        genus = self.genus()
        assert G_plus.degree() <= (genus + 1)
        assert (G_plus**2 + h * G_plus - f).degree() <= genus
        assert G_minus.leading_coefficient() == alpha_minus

        self._infinte_polynomials = G_plus, G_minus
        return self._infinte_polynomials

    def points_at_infinity(self):
        """
        Compute the points at infinity on the curve. Assumes we are using
        a weighted projective model for the curve
        """
        # TODO: check to False
        return [self.point([1, y, 0], check=True) for y in self.roots_at_infinity()]

    def is_x_coord(self, x):
        """
        Return True if ``x`` is the `x`-coordinate of a point on this curve.
        """
        f, h = self.hyperelliptic_polynomials()
        K = self.base_ring()
        try:
            x = K(x)
        except (ValueError, TypeError):
            raise TypeError('x must be coercible into the base ring of the curve')

        # When h is zero then x is a valid coordinate if y2 is square
        if not h:
            y2 = f(x)
            return y2.is_square()
        # Generic case for h != 0
        a = f(x)
        b = h(x)
        # Special case for char 2
        if K.characteristic() == 2:
            R = f.parent()  # Polynomial ring K[x]
            F = R([-a, b, 1])
            return bool(F.roots())
        # Otherwise x is a point on the curve if the discriminant is a square
        D = b*b + 4*a
        return D.is_square()

    def lift_x(self, x, all=False):
        """
        Return one or all points with given `x`-coordinate.

        This method is deterministic: It returns the same data each
        time when called again with the same `x`.
        """
        f, h = self._hyperelliptic_polynomials
        K = self.base_ring()

        # Compute the common parent between the base ring of the curve and
        # the parent of the input x-coordinate.
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        try:
            L = cm.common_parent(x.parent(), K)
            x = L(x)
        except (TypeError, ValueError):
            raise ValueError('x must have a common parent with the base ring')

        # First we compute the y-coordinates the given x-coordinate
        ys = []
        one = L.one()

        # When h is zero we find all y-coordinates with a single sqrt
        if not h:
            y2 = f(x)
            # When y2 is not a square, ys will be an empty list
            ys = y2.sqrt(all=True, extend=False)
        # Otherwise we need roots of the discriminant
        else:
            a = f(x)
            b = h(x)
            # Special case for char 2
            if K.characteristic() == 2:
                R = f.parent()
                F = R([-a, b, 1])
                ys = F.roots(L, multiplicities=False)
            else:
                D = b*b + 4*a
                # When D is not a square, ys will be an empty list
                ys = [(-b+d)/2 for d in D.sqrt(all=True, extend=False)]

        if ys:
            ys.sort()  # Make lifting deterministic
            if all:
                return [self.point([x, y, one], check=False) for y in ys]
            else:
                return self.point([x, ys[0], one], check=False)

        if all:
            return []
        else:
            raise ValueError(f"No point with x-coordinate {x} on {self}")

    def distinguished_point(self):
        """
        Return the distinguished point of the hyperelliptic curve.
        By default, this is one of the points at infinity if possible.
        """
        if self._distinguished_point is None:
            if not self.is_inert():
                # For the the split and ramified case, a point at infinity is chosen,
                self._distinguished_point = self.points_at_infinity()[0]
            else:
                assert self.base_ring().characteristic() > 0, "in characteristic 0, a distinguished_point needs to be specified"
                #in the inert case we choose a point with minimal x-coordinate
                for x0 in self.base_ring():
                    try:
                        self._distinguished_point = self.lift_x(x0)
                        break
                    except ValueError:
                        pass

        return self._distinguished_point

    def set_distinguished_point(self, P0):
        """
        Change the distinguished point of the hyperelliptic curve to P0.
        """
        assert isinstance(P0, HyperellipticPoint), "the input has to be a point on the curve"
        self._distinguished_point = P0
        return None


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

    def __call__(self, *args):
        return self.point(*args)

    def jacobian(self):
        if self.is_ramified():
            from jacobian_ramified import HyperellipticJacobianRamified
            return HyperellipticJacobianRamified(self)
        if self.is_split():
            from jacobian_split import HyperellipticJacobianSplit
            return HyperellipticJacobianSplit(self)

        assert self.is_inert()
        from jacobian_inert import HyperellipticJacobianInert
        return HyperellipticJacobianInert(self)

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
