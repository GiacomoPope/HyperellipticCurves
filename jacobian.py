from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.rings.integer import Integer
from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import choice, randint

from hyperelliptic import HyperellipticCurveNew, HyperellipticPoint

# Needed until https://github.com/sagemath/sage/pull/37118 is merged.
from uniform_random_sampling import uniform_random_polynomial

class HyperellipticJacobian:
    def __init__(self, H):
        if not isinstance(H, HyperellipticCurveNew):
            raise ValueError("TODO")
        self._curve = H
        self._element = MumfordDivisor

    def __repr__(self) -> str:
        return f"Jacobian of {self._curve}"

    def curve(self):
        return self._curve

    def base_ring(self):
        return self._curve.base_ring()
    
    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        R = self._curve.polynomial_ring()
        return self._element(self, R.one(), R.zero())

    @cached_method
    def cardinality(self):
        """
        TODO: currently using lazy methods by calling sage
        """
        return sum(self.curve().frobenius_polynomial())

    order = cardinality

    def point_to_mumford_coordinates(self, P):
        """
        TODO
        """
        R, x = self._curve.polynomial_ring().objgen()
        [X, Y, Z] = P.coords()        
        if Z == 0:
            return R.one(), R.zero()
        u = x - X
        v = R(Y)
        return u, v

    def __call__(self, *args, check=True):
        if isinstance(args[0], HyperellipticPoint) and len(args) == 1:
            u, v = self.point_to_mumford_coordinates(args[0])
        # TODO handle this better!!
        elif len(args) == 1:
            return self.zero()
        elif len(args) == 2:
            u, v = args
        else:
            raise NotImplementedError
        return self._element(self, u, v, check=check)

    def __cantor_double_generic(self, u1, v1):
        """
        Efficient cantor composition for doubling an affine divisor

        Returns the Cantor composition of (u1, v1) with (u1, v1)
        together with the degree of the polynomial ``s`` which is
        needed for computing weights for the split and inert models 
        """
        f, h = self._curve.hyperelliptic_polynomials()

        # New mumford coordinates
        if h.is_zero():
            s, _, e2 = u1.xgcd(v1 + v1)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1**2) // s
        else:
            s, _, e2 = u1.xgcd(v1 + v1 + h)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1 * h - v1**2) // s
        v3 = v3 % u3

        return u3, v3, s.degree()

    def _cantor_composition_generic(self, u1, v1, u2, v2):
        """
        Cantor composition

        Returns the Cantor composition of (u1, v1) with (u2, v2)
        together with the degree of the polynomial ``s`` which is
        needed for computing weights for the split and inert models 
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Ensure D1 and D2 are semi-reduced divisors
        assert (
            v1.degree() < u1.degree() and v2.degree() < u2.degree()
        ), "The degree of bi must be smaller than ai"
        assert (
            u1.degree() <= 2*g+2 and u2.degree() <= 2*g+2
        ), f"The degree of ai must be smaller than 2g+2, {u1.degree()}, {u2.degree()}"

        # Special case: duplication law
        if u1 == u2 and v1 == v2:
            return self.__cantor_double_generic(u1, v1)

        # Step One
        s0, _, e2 = u1.xgcd(u2)
        v1_m_v2 = (v1 - v2)

        # Special case: when gcd(u0, u1) == 1 we can
        # avoid many expensive steps as we have s = 1
        if s0.is_one():
            u3 = u1 * u2
            v3 = v2 + e2 * u2 * v1_m_v2
            v3 = v3 % u3
            return u3, v3, 0

        # Step Two
        w0 = v1 + v2 + h

        # Another special case, when w0 is zero we skip
        # a xgcd and can return early
        if w0.is_zero():
            u3 = (u1 * u2) // (s0**2)
            v3 = v2 + e2 * v1_m_v2 * (u2 // s0)
            v3 = v3 % u3
            return u3, v3, s0.degree()

        # Step Three
        s, c1, c2 = s0.xgcd(w0)
        u3 = (u1 * u2) // (s**2)
        v3 = v2 + (c1 * e2 * v1_m_v2 * u2 + c2 * (f - h * v2 - v2**2)) // s
        v3 = v3 % u3
        return u3, v3, s.degree()

    def _cantor_reduction_generic(self, u0, v0):
        """
        TODO

        NOTE: we do not make u1 monic, but leave this for *after* reduction
        to save on calls to this method.
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()

        # Compute u' and v'
        u1 = (v0**2 + h * v0 - f) // u0
        v1 = (-h - v0) % u1

        return u1, v1
    
    def cantor_composition(self, u1, v1, u2, v2):
        u3, v3, _ = self._cantor_composition_generic(u1, v1, u2, v2)
        return u3, v3

    def cantor_reduction(self, u0, v0):
        return self._cantor_reduction_generic(u0, v0)

    def _random_element_cover(self, degree=None):
        r"""
        Returns a random element from the Jacobian. Distribution is not
        uniformly random, but returns the entire group.
        """
        H = self.curve()
        K = H.base_ring()
        R = H.polynomial_ring()
        g = H.genus()

        f, h = H.hyperelliptic_polynomials()
        
        # For the inert case, the genus must be even
        if H.is_inert():
            assert not (g % 2)

        if degree is None:
            degree = (-1, g)

        while True:
            u = uniform_random_polynomial(R, degree=degree)
            if u.is_zero():
                if H.is_split():
                    n = randint(0, g)
                    return self._element(self, R.one(), R.zero(), n)
                return self.zero()
            
            u = u.monic()

            # TODO: i think we can skip this and simply ensure u
            #       is even degree with composition with the distinguished
            #       point?
            # if H.is_inert() and (u.degree() % 2) == 1:
            #     #TODO: better method to sample even degree polynomials
            #     continue
            
            try:
                u1, v1 = R.one(), R.zero()
                for x, e in u.factor():
                    # Solve y^2 + hy - f = 0 mod x
                    from sage.rings.polynomial.polynomial_ring import polygen

                    # TODO: is this the most efficient method? Maybe we should write
                    # a helper function which computes y^2 + hy - f = 0 mod x which
                    # properly handles trivial cases like when x is linear? 
                    K_ext = K.extension(modulus=x, names="a")
                    y_ext = polygen(K_ext, "y_ext")
                    h_ = K_ext(h % x)
                    f_ = K_ext(f % x)
                    y = choice((y_ext**2 + h_ * y_ext - f_).roots(multiplicities=False))
                    try:
                        # Attempt to coerce quotient ring element to the
                        # polynomial ring
                        v = R(y)

                        # Sum for the multiplicity of the root x of u
                        for _ in range(e):
                            u1, v1, _ = self._cantor_composition_generic(u1, v1, x, v)

                    # v is not rational, so we skip it
                    except (ValueError, AttributeError):
                        pass

                if H.is_split():
                    g = self._curve.genus()
                    n = randint(0, g - u1.degree())
                    return self._element(self, u1, v1, n, check=False)
                
                # We need to ensure the degree of u is even
                if H.is_inert():
                    if (u1.degree() % 2):
                        # TODO: make composition with distinguished_point
                        #       its own function?
                        P0 = self._curve.distinguished_point()
                        X0, Y0 = P0.xy()
                        X = R.gen() # TODO use better variable names in this function
                        _, h = self._curve.hyperelliptic_polynomials()
                        u0 = X - X0
                        v0 = R(-Y0 - h(X0))
                        u1, v1, _ = self._cantor_composition_generic(u1, v1, u0, v0)
                    assert not (u1.degree() % 2), f"{u1} must have even degree"
                return self._element(self, u1, v1, check=False)
            
            # TODO: better handling rather than looping with try / except?
            except IndexError:
                pass

    def _random_element_rational(self):
        r"""
        Returns a random element from the Jacobian. It does not necessarily
        return the entire group.
        """
        H = self.curve()
        g = H.genus()

        # We randomly sample 2g + 1 points on the hyperelliptic curve
        points = [H.random_point() for _ in range(2 * g + 1)]

        # We create 2g + 1 divisors of the form (P) - infty
        divisors = [self(P) for P in points if P[2] != 0]

        # If we happened to only sample the point at infinity, we return this
        # Otherwise we compute the sum of all divisors.
        if not divisors:
            return self.zero()
        return sum(divisors, start=self.zero())

    def random_element(self, fast=True, *args, **kwargs):
        r"""
        Returns a random element from the Jacobian. Distribution is not
        uniformly random.

        INPUT:

        - ``fast`` -- (boolean, default ``True``) If set to ``True``, a fast
          algorithm is used, but the output is **NOT** guaranteed to cover the
          entire Jacobian. See examples below. If set to ``False``, a slower
          algorithm is used, but covers the entire Jacobian.
        """
        if not isinstance(self.base_ring(), FiniteField_generic):
            raise NotImplementedError(
                "random element of Jacobian is only implemented over Finite Fields"
            )

        if fast:
            return self._random_element_rational(*args, **kwargs)
        return self._random_element_cover(*args, **kwargs)

class MumfordDivisor:
    def __init__(self, parent, u, v, check=True):
        if not isinstance(parent, HyperellipticJacobian):
            raise TypeError(f"parent must be of type {HyperellipticJacobian}")
        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")
 
        # TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))

        # Ensure the divisor is valid
        if check:
            f, h = parent.curve().hyperelliptic_polynomials()
            assert (
                v**2 + v * h - f
            ) % u == 0, f"{u = }, {v = } do not define a divisor on the Jacobian"

        self._parent = parent
        self._u = u
        self._v = v

    def parent(self):
        return self._parent

    def __repr__(self) -> str:
        return f"({self._u}, {self._v})"
    
    def is_zero(self):
        return self._u.is_one() and self._v.is_zero()

    def uv(self):
        return (self._u, self._v)

    def __eq__(self, other):
        if not isinstance(other, MumfordDivisor):
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __hash__(self):
        data = (self._u, self._v)
        return hash(data)

    @cached_method
    def order(self):
        n = self.parent().order()
        return order_from_multiple(self, n)

    def degree(self):
        """
        Returns the degree of the affine part of the divisor.
        """
        return self._u.degree()

    def __add__(self, other):
        # Ensure we are adding two divisors
        if not isinstance(other, type(self)):
            raise ValueError("TODO")

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Step one: cantor composition of the two divisors
        u3, v3 = self._parent.cantor_composition(u1, v1, u2, v2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > g:
            u3, v3 = self._parent.cantor_reduction(u3, v3)
        u3 = u3.monic()
        v3 = v3 % u3

        return self._parent._element(self._parent, u3, v3, check=False)

    def __neg__(self):
        _, h = self.parent().curve().hyperelliptic_polynomials()
        u0, v0 = self.uv()

        if h.is_zero():
            return self._parent._element(self._parent, u0, -v0, check=False)

        v1 = (-h - v0) % u0
        return self._parent._element(self._parent, u0, v1, check=False)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, n):
        """ """
        # TODO: is there a better handlings for this?
        if not isinstance(n, (int, Integer)):
            raise ValueError

        if not n:
            return self._parent().zero()

        P = self

        # Handle negative scalars
        if n < 0:
            n = -n
            P = -P

        # Double and Add
        Q = P
        R = self.parent().zero()
        while n > 0:
            if n % 2 == 1:
                R = R + Q
            Q = Q + Q
            n = n // 2
        return R

    __rmul__ = __mul__
