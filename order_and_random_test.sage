from hyperelliptic_split import HyperellipticCurveSplit
R.<x> = PolynomialRing(GF(3))

def random_sample(J, n=2000, fast=True):
    p = [J.random_element(fast=fast) for _ in range(n)]
    return len(set(p))

def random_curve(use_h=True, genus=2):
    d = 2*genus + 2
    while True:
        while True:
            f = R.random_element(degree=d)
            if f.degree() == d:
                f = f.monic()
                break
        if use_h:
            h = R.random_element(degree=2)
        else:
            h = 0
        try:
            HyperellipticCurve(f, h)
            return f, h, HyperellipticCurveSplit(f, h)
        except:
            continue

# Test that randomly sampling gets all elements in the group
for _ in range(3):
    f, h, H = random_curve(genus=2)
    J = H.jacobian()
    o = J.order()

    fast = random_sample(J)
    slow = random_sample(J, fast=False)

    print(f"{f = }")
    print(f"{h = }")
    print(f"Order: {o}")
    print(f"Fast method: {fast}")
    print(f"Slow method: {slow}")
    print(f"")

# Test all points have order dividing the Jacobian order
# for _ in range(1):
#     f, h, H = random_curve(genus=3)
#     J = H.jacobian()
#     o = J.order()
#     assert all([(o * J.random_element()).is_zero() for _ in range(100)])

# Test inversion works
for _ in range(10):
    f, h, H = random_curve(genus=4)
    J = H.jacobian()
    o = J.order()

    D = J.random_element()
    assert all([(D - D).is_zero() for _ in range(100)])


# Test all points have order dividing the Jacobian order
R.<x> = PolynomialRing(GF(3))
g = 3
for _ in range(3):
    f, h, H = random_curve(use_h=False, genus=g)
    print(f"{f.factor() = }")
    J = H.jacobian()
    o = J.order()

    print(f"Testing with: {J}")
    bad = []
    for _ in range(1000):
        D = J.random_element()
        if not (o * D).is_zero():
            bad.append(D)
    print(f"Bad elements of 1000: {len(bad)}")
    print(f"Number of unique bad elements: {len(set(bad))}")
    print(f"Order of the Jacobian: {o.factor()}")

    print()
