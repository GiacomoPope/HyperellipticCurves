import hyperelliptic_generic


class HyperellipticCurveSmoothModel_rational_field(
    hyperelliptic_generic.HyperellipticCurveSmoothModel_generic
):
    def __init__(self, projective_model, f, h, genus):
        super().__init__(projective_model, f, h, genus)
