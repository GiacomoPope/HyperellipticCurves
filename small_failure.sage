from hyperelliptic_split import HyperellipticCurveSplit
R.<x> = PolynomialRing(GF(3))
f = x^8 + 2*x^6 + 2*x^5 + x^4 + 2*x^3 + x^2 + 1
h = x^2 + 2*x + 2
H = HyperellipticCurveSplit(f, h)
J = H.jacobian()
D = J(x^3 + x^2 + x + 1, x + 2, 0)
3*D + D
4*D
2*D + 2*D
