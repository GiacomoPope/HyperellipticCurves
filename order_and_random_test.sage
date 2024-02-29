from hyperelliptic_split import HyperellipticCurveSplit
R.<x> = PolynomialRing(GF(7))

def random_sample(J, n=2000, fast=True):
    p = []
    for _ in range(n):
        p.append(J.random_element(fast=fast))
    return len(set(p))

def random_curve(use_h=True):
    while True:
        while True:
            f = R.random_element(degree=6)
            if f.degree() == 6:
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
for _ in range(1):
    f, h, H = random_curve()
    J = H.jacobian()
    o = J.order()

    fast = random_sample(J)
    slow = random_sample(J, n=5000, fast=False)

    print(f"{f = }")
    print(f"{h = }")
    print(f"Order: {o}")
    print(f"Fast method: {fast}")
    print(f"Slow method: {slow}")
    print(f"")

# Test all points have order dividing the Jacobian order
for _ in range(1):
    f, h, H = random_curve()
    J = H.jacobian()
    o = J.order()
    assert all([(o * J.random_element()).is_zero() for _ in range(100)])
