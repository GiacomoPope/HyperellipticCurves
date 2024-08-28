from sage.schemes.generic.homset import SchemeHomset_points

# This file imitates `projective/projective_homset.py`

class SchemeHomset_points_weighted_projective_field(SchemeHomset_points):
    def points(self):
        raise NotImplementedError
