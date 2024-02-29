from sage.misc.prandom import choice
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.functions.other import binomial
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer

class HyperellipticCurveSplit:
    def __init__(self, f, h=0):
        # Some values which we will cache as a user asks for them
        self._alphas = None
        self._genus = None
        self._projective_model = None
        self._infinte_polynomials = None

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

        # Class which represents points
        self._point = HyperellipticPoint

    def __repr__(self):
        f, h = self._hyperelliptic_polynomials
        if h:
            if h.is_one():
                curve = f"y^2 + y = {f}"
            else:
                curve = f"y^2 + y*{h} = {f}"
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

    def projective_model(self):
        """
        Compute the weighted projective model (1 : 3 : 1)
        """
        if self._projective_model:
            return self._projective_model

        f, h = self._hyperelliptic_polynomials
        R, (X, Y, Z) = PolynomialRing(self.base_ring(), names="X, Y, Z").objgens()

        if h is None:
            d = f.degree()
            F = sum(f[i] * X**i * Z**(d - i) for i in range(d + 1))
            return Y**2 - F

        d = max(h.degree(), (f.degree() / 2).ceil())
        H = sum(h[i] * X**i * Z**(d - i) for i in range(d + 1))
        F = sum(f[i] * X**i * Z**(2*d - i) for i in range(2*d + 1))
        self._projective_model = Y**2 + H*Y - F
        return self._projective_model

    def point(self, coords, check=True):
        if len(coords) == 2:
            X, Y = coords
            Z = self.base_ring().one()
        elif len(coords) == 3:
            X, Y, Z = coords
        else:
            raise ValueError("TODO")

        if check:
            F = self.projective_model()
            if not F(X, Y, Z).is_zero():
                raise ValueError("TODO")

        return self._point(X, Y, Z)

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
        if not h:
            return f.leading_coefficient().sqrt(all=True)

        d = self._d
        self._alphas = (x**2 + x * h[d] - f[2*d]).roots(multiplicities=False)
        return self._alphas

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
            if not i:
                # TODO: deal with that there is more than one point at infinity
                return choice(self.points_at_infinity())
            v = self.lift_x(k.random_element(), all=True)
            try:
                return v[i % 2]
            except IndexError:
                pass

    def __call__(self, coords):
        return self.point(coords)

    def jacobian(self):
        from jacobian_split import JacobianSplit
        return JacobianSplit(self)

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

class HyperellipticPoint:
    def __init__(self, X, Y, Z=1):
        self._X = X
        self._Y = Y
        self._Z = Z

    def coords(self):
        return (self._X, self._Y, self._Z)

    def xy(self):
        if not self._Z:
            raise ValueError("TODO")

        x = self._X / self._Z
        y = self._Y / self._Z
        return x, y

    def __repr__(self):
        return f"({self._X} : {self._Y} : {self._Z})"

    def __eq__(self, other):
        if not isinstance(other, HyperellipticPoint):
            return False
        # TODO: this is dumb but works fine as Z is always either 0 or 1 currently
        return self.coords() == other.coords()

    def __hash__(self):
        return hash(self.coords())

    def __getitem__(self, n):
        try:
            return self.coords()[n]
        except IndexError:
            raise IndexError("todo")
