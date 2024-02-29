from hyperelliptic_split import HyperellipticCurveSplit

R.<x> = PolynomialRing(GF(11))
f = x^6 + 6*x^5 + 4*x^4 + 3*x^3 + 8*x^2 + 4*x + 1

H = HyperellipticCurveSplit(f)

J = H.jacobian()

D1 = J(x^2,2*x+1)
D2 = J(x,10+x*0)

# TODO:
# Sometimes J.random_element(fast=True) errors as the divisor is not valid
# We also have over counting, so something is definitely wrong
hmm = []
for _ in range(5000):
    D1 = J.random_element()
    hmm.append(D1)
print(f"{len(set(hmm)) = }")
print(f"{J.order() = }")

# The slow method seems OK though
hmm = []
for _ in range(5000):
    try:
        D1 = J.random_element(fast=False)
    except:
        pass
    hmm.append(D1)
print(f"{len(set(hmm)) = }")
print(f"{J.order() = }")
