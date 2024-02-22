from hyperelliptic_split import HyperellipticCurveSplit

F = GF(97)
R.<x> = PolynomialRing(F)

# Example curve and divisors D1, D2 from Example 3.3.1
f = x^6 + 13*x^2 + 92*x + 7
H = HyperellipticCurveSplit(f)
J = H.jacobian()

u1 = x^2 + 75*x + 57
v1 = x + 13

u2 = x^2 + 38*x + 41
v2 = x + 25

D1 = J(u1, v1)
D2 = J(u2, v2)

D3, (omega_plus, omega_minus) = D1.cantor_composition(D2)
u3, v3 = D3.uv()

assert u3 == x^4 + 16*x^3 + 38*x^2 + 3*x + 9
assert v3 == 20*x^3 + 2*x^2 + 50*x + 84
assert omega_plus == omega_minus == 0

D4, (omega_plus, omega_minus) = D3.cantor_reduction()
u4, v4 = D4.uv()

assert u4 == x^2 + 53*x + 81
assert v4 == 10*x + 63
assert omega_plus == omega_minus == 1