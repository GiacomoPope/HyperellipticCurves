from hyperelliptic_split import HyperellipticCurveSplit

F = GF(127)
R.<x> = PolynomialRing(F)

# Example curve and divisors D1, D2 from Example 3.2.1
f = x^8 + 2*x^5 + x^4 + 4*x^2 + 88*x + 45
H = HyperellipticCurveSplit(f)
J = H.jacobian()

u1 = x^3 + 35*x^2 + 47*x + 51
v1 = 68*x^2 + x + 41
D1 = J(u1, v1)

u2 = x^2 + 121*x + 100
v2 = 37*x + 113
D2 = J(u2, v2)

D3, (omega_plus, omega_minus) = D1.cantor_composition(D2)
u3, v3 = D3.uv()

# Ensure cantor composition matches D3 from Example 3.2.1
assert u3 == x^5 + 29*x^4 + 64*x^3 + 94*x^2 + 76*x + 20
assert v3 == 6*x^4 + 40*x^3 + 115*x^2 + 64*x + 7
assert omega_plus == omega_minus == 0

# Example 3.2.2
D4, (omega_plus, omega_minus) = D3.cantor_reduction()
u4, v4 = D4.uv()

# Ensure cantor reduction matches D4 from Example 3.2.2
assert u4 == x^3 + 21*x^2 + 29*x + 45
assert v4 == 31*x^2 + 125*x + 60
assert omega_plus == omega_minus == 1

# Ensure reduction at infinity at G^+ matches example 3.2.3
#
D5, (omega_plus, omega_minus) = D4.cantor_reduction_at_infinity()
u5, v5 = D5.uv()

assert u5 == x^2 + 8*x + 57
assert v5 == 39*x + 14
assert omega_plus == -1
assert omega_minus == 2
