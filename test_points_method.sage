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
R.<x> = GF(17)[]
H = HyperellipticCurveSmoothModel(x^5 + 5 * x + 1)

J = H.jacobian()
assert J.order() == 360

JH = J.point_homset()
assert JH.order() == 360

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
try_exec("assert J(K.one(), K.zero()) == 0")
try_exec("assert J(K.one(), 0) == 0")
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
try_exec("J.points()")
try_exec("J.rational_points()")
try_exec("list(J)")
try_exec("JH.points()")
try_exec("JH.rational_points()")
# try_exec("list(JH)")  # hangs, #38566

### J and JH are two different things but are printed the same
try_exec("assert type(J) != type(JH)")
try_exec("assert str(J) != str(JH)")

### base ring is not changed correctly
K2 = H.base_ring().extension(2)
JH2 = J(K2)
try_exec("assert JH2 != JH")
try_exec("assert JH2.base_ring() == K2")
