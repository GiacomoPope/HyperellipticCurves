from sage.schemes.generic.morphism import SchemeMorphism, SchemeMorphism_point

# This file imitates `projective/projective_point.py`

class SchemeMorphism_point_weighted_projective_ring(SchemeMorphism_point):
    def __init__(self, X, v, check=True):
        SchemeMorphism.__init__(self, X)
        self._coords = tuple(v)
        self._normalized = False

    def _repr_(self):
        return "({})".format(" : ".join(map(repr, self._coords)))

    def points(self):
        pass
