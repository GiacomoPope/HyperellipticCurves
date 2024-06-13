from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.projective.projective_space import ProjectiveSpace, _CommRings
from sage.structure.all import UniqueRepresentation
from sage.structure.category_object import normalize_names

try:
    # TODO: Remove this
    from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
except ImportError:
    print("Please merge #38207 to your local Sage installation.")
    from sage.rings.polynomial.polynomial_ring import (
        PolynomialRing_general as PolynomialRing_generic,
    )

    # exit(1)


def WeightedProjectiveSpace(weights, R=None, names=None):
    r"""
    Return a weighted projective space with the given ``weights`` over the ring ``R``.

    EXAMPLES::

    TODO
    """
    if isinstance(weights, (MPolynomialRing_base, PolynomialRing_generic)) and R is None:
        if names is not None:
            # Check for the case that the user provided a variable name
            # That does not match what we wanted to use from R
            names = normalize_names(weights.ngens(), names)
            if weights.variable_names() != names:
                # The provided name doesn't match the name of R's variables
                raise NameError(
                    "variable names passed to ProjectiveSpace conflict with names in ring"
                )
        A = WeightedProjectiveSpace(
            weights.ngens() - 1, weights.base_ring(), names=weights.variable_names()
        )
        A._coordinate_ring = weights
        return A

    if isinstance(R, (int, Integer, list, tuple)):
        weights, R = R, weights
    elif R is None:
        R = ZZ

    # WeightedProjectiveSpace(5) -> just return unweighted version
    if isinstance(weights, (int, Integer)):
        return ProjectiveSpace(weights, R=R, names=names)

    if names is None:
        names = "x"

    # TODO: Specialise implementation to projective spaces over non-rings
    # But since we don't really implement extra functionalities, I don't think
    # we care.
    if R in _CommRings:
        return WeightedProjectiveSpace_ring(weights, R=R, names=names)

    raise TypeError(f"R (={R}) must be a commutative ring")


class WeightedProjectiveSpace_ring(UniqueRepresentation, AmbientSpace):
    @staticmethod
    def __classcall__(cls, weights, R=ZZ, names=None):
        # __classcall_ is the "preprocessing" step for UniqueRepresentation
        # see docs of CachedRepresentation
        if not isinstance(weights, list):
            raise TypeError(
                f"weights(={weights}) is not a list. Please use the `WeightedProjectiveSpace`"
                " constructor"
            )
        normalized_names = normalize_names(len(weights) + 1, names)
        return super().__classcall__(cls, weights, R, normalized_names)

    def __init__(self, weights, R=ZZ, names=None):
        """
        Initialization function.

        EXAMPLES::

            sage: WeightedProjectiveSpace([1, 3, 1], Zp(5), 'y')                        # needs sage.rings.padics
            Weighted Projective Space of dimension 3 with weights (1, 3, 1) over 5-adic Ring with
            capped relative precision 20
            sage: WeightedProjectiveSpace(5, QQ, 'y')
            Projective Space of dimension 5 over Rational Field
        """
        AmbientSpace.__init__(self, len(weights), R)
        self.weights = weights
        self._assign_names(names)

    def ngens(self):
        """
        Return the number of generators of this weighted projective space.

        This is the number of variables in the coordinate ring of ``self``.

        EXAMPLES::

            sage: WeightedProjectiveSpace([1, 3, 1], QQ).ngens()
            4
            sage: WeightedProjectiveSpace(5, ZZ).ngens()
            6
        """
        return self.dimension_relative() + 1

