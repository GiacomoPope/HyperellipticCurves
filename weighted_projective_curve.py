from sage.schemes.curves.curve import Curve_generic
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme

class WeightedProjectiveCurve(Curve_generic, AlgebraicScheme_subscheme):
    def __init__(self, A, X, category=None):
        # TODO: ProjectiveCurve allows a category, which must
        # be allows for them, but Curve_generic crashes if we
        # try and send a category kwarg
        # TODO ensure that A is the right type?
        # Something like a `is_WeightProjectiveSpace` which means making a
        # WeightProjectiveSpace class  lmao
        Curve_generic.__init__(self, A, X)
