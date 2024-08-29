from hyperelliptic_constructor import HyperellipticCurveSmoothModel

R.<x> = PolynomialRing(GF(5))

def random_sample(J, n=500, fast=True):
    p = [J.random_element(fast=fast) for _ in range(n)]
    return len(set(p))

def random_curve(use_h=True, genus=2):
    d = 2*genus + 2
    while True:
        # Find a polynomial for f
        while True:
            f = R.random_element(degree=d)
            if f.degree() == d:
                f = f.monic()
                break
        # Find a polynomial for h
        if use_h:
            if R.base_ring().characteristic() == 2:
                h = R.random_element(degree=(genus + 1))
            else:
                h = R.random_element(degree=genus)
        else:
            h = 0

        # Ensure that there are two points at infinity and the curve is non-singular
        try:
            H = HyperellipticCurveSmoothModel(f, h)
            if not H.is_split():
                continue
            return f, h, HyperellipticCurveSmoothModel(f, h)
        except Exception as e:
            # print(f"{e = }")
            continue

# Test that randomly sampling gets all elements in the group
for _ in range(1):
    f, h, H = random_curve(genus=1)
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


for _ in range(1):
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
for g in [1, 2, 3, 4, 5]:
    print(f"Testing arithmetic for genus: {g}")
    for _ in range(5):
        f, h, H = random_curve(genus=g)
        J = H.jacobian()
        o = J.order()

        # Test order
        assert all([(o * J.random_element()).is_zero() for _ in range(100)])
        assert all([(o * J.random_element(fast=False)).is_zero() for _ in range(100)])

        # Test order on divisor
        for _ in range(100):
            D = J.random_element()
            order = D.order()
            assert order * D == J.zero()

        # Test inversion
        for _ in range(100):
            D = J.random_element()
            assert (D - D).is_zero()

        # Test inversion on non-rational
        for _ in range(100):
            D = J.random_element(fast=False)
            assert (D - D).is_zero()

    print(f"Passed for genus: {g}")

