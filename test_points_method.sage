def try_exec(s):
    try:
        exec(s)
        print(f"\x1b[32mSuccessful! ({s})\x1b[0m")
    except AttributeError:
        print(f"\x1b[31mAttributeError: {s}\x1b[0m")
    except AssertionError:
        print(f"\x1b[31mAssertionError: {s}\x1b[0m")
    except Exception as e:
        print(f"\x1b[31m[?] {type(e).__name__}: {s} (\x1b[34merror: {e}\x1b[31m)\x1b[0m")

from hyperelliptic_constructor import HyperellipticCurveSmoothModel

### Setup
K = GF(17)
R.<x> = K[]

for curve, order, order2 in [
    (x^5 + 5 * x + 1, 360, 90720),
    (x^6 + 5 * x + 1, 231, 86625),
    (3*x^6 + 5 * x + 1, 333, 86913),
]:
    print(f"\x1b[33mTesting curve y^2 = {curve}\x1b[0m")
    H = HyperellipticCurveSmoothModel(curve)

    J = H.jacobian()
    assert J.order() == order

    JH = J.point_homset()
    assert JH.order() == order

    ### unique constructor
    try_exec("assert H.jacobian() == H.jacobian()")
    try_exec("assert H.jacobian() is H.jacobian()")
    try_exec("assert J.point_homset() == J.point_homset()")
    try_exec("assert J.point_homset() is J.point_homset()")

    ### Check richcmp implementation
    P = J.random_element()
    try_exec("assert P == P")
    try_exec("assert not (P != P)")
    Z = J(0)
    try_exec("assert Z == 0")
    try_exec("assert not (Z != 0)")
    try_exec("assert Z != 1")
    try_exec("assert not (Z == 1)")

    ### Constructor fucks up
    try_exec("assert J(R.one(), R.zero()) == 0")
    try_exec("assert J(R.one(), 0) == 0")
    # this doesn't run
    try_exec("assert J(1, 0) == 0")
    try_exec("assert J() == 0")

    ### Arithmetic fucks up
    try_exec("assert 0 * P == 0")
    # this doesn't even run
    try_exec("assert P * 0 == 0")

    ### Ideally this should work
    from random import randrange
    from sage.groups.generic import discrete_log

    while (G := J.random_element()).order() != J.order():
        continue
    k = randrange(J.order())
    P = G * k
    try_exec("assert discrete_log(P, G, operation='+') == k")

    ### Some of these haven't been implemented
    try_exec("assert len(J.points()) == J.order()")
    # try_exec("assert J.points() == J.rational_points() == list(J)")
    # try_exec("assert J.points() == J.points()")
    try_exec("assert len(JH.points()) == JH.order()")
    # try_exec("assert JH.points() == JH.rational_points() == list(JH)")
    # try_exec("assert JH.points() == JH.points()")

    ### J and JH are two different things but are printed the same
    try_exec("assert type(J) != type(JH)")
    try_exec("assert str(J) != str(JH)")

    ### base ring is not changed correctly
    K2 = H.base_ring().extension(2)
    JH2 = J(K2)
    try_exec("assert JH2 != JH")
    # note this shouldn't be K2
    try_exec("assert JH2.base_ring() == K")
    try_exec("assert JH2.curve() == H")
    try_exec("assert JH2.extended_curve() == H.base_extend(K2)")
    try_exec("assert JH2.codomain() == J")
    try_exec("assert JH2.extended_codomain() == J.base_extend(K2)")

    ### counting points over an extension (using Weil's conjectures)
    try_exec("assert J.count_points(1) == order")
    try_exec("assert J.count_points(2) == [order, order2]")
    try_exec("assert J.count_points(10)[1::2] == JH2.count_points(5)")

    ### base change correctly
    J2 = J.change_ring(K2)
    try_exec("assert J2.order() == order2")
    JH2_ = J2.point_homset()
    try_exec("assert JH2_ != JH")
    try_exec("assert JH2_.base_ring() == K2")

    print("~" * 40)
