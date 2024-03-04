from hyperelliptic_split import HyperellipticCurveSplit

F = GF(211)
R.<x> = PolynomialRing(F)

f = x^8 + 53*x^5 + 158*x^4 + 12*x^3 + x + 187
H = HyperellipticCurveSplit(f)
J = H.jacobian()

u1 = x^3 + 40*x^2 + 28*x + 134
v1 = x^4 + 91*x^2 + 143*x + 92 # v1 has a weird degree?
v1 = v1 % u1

u2 = x^3 + 110*x^2 + 104*x + 197
v2 = x^4 + 93*x^2 + 52*x + 50
v2 = v2 % u2

(u3, v3), (omega_plus, omega_minus) = J.cantor_composition(u1, v1, u2, v2)

assert u3 == x^6 + 150*x^5 + 101*x^4 + 186*x^3 + x^2 + 40*x + 23
assert v3 == 47*x^5 + 169*x^4 + 155*x^3 + 209*x^2 + 161*x + 166
assert omega_plus == omega_minus == 0

(u4, v4), (omega_plus, omega_minus) = J.cantor_reduction(u3, v3)

assert u4 == x^4 + 149*x^3 + 129*x^2 + 152*x + 198
assert v4 == 20*x^3 + 155*x^2 + 57*x + 56
assert omega_plus == omega_minus == 1

(u5, v5), (omega_plus, omega_minus) = J.cantor_compose_at_infinity(u4, v4, plus=True)
assert u5 == x^3 + 195*x^2 + 181*x + 5
assert v5 == 102*x^2 + 154*x + 38
# assert omega_plus == 1
# assert omega_minus == 0
print(omega_plus, omega_minus)
