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
for genus in [2, 3, 4, 5]:
    for _ in range(10):
        f, h, H = random_curve(genus=3)
        J = H.jacobian()
        o = J.order()

        # Test order
        assert all([(o * J.random_element()).is_zero() for _ in range(100)])

        # Test inversion
        for _ in range(100):
            D = J.random_element()
            assert (D - D).is_zero()

