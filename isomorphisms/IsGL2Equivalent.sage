from itertools import *

def equivalence(f1,f2,n):

    # Checking the ground field of the input polynomials
    assert f1.base_ring() == f2.base_ring(), "Input polynomials must be defined over the same field"

    # Base field
    FF = f1.base_ring()
    assert FF.is_field() & FF.is_finite(), "The base filed must be a finite field"


    # DEFINING SOME ADDITIONAL FUNCTIONS FIRST.
    # LIKELY WE ARE NOT GOING TO NEED THIS. HENCE WE POSTPONE.

    def adding_zeros(input, length):
        assert len(input) <= length

        return (input + length*[FF.zero()])[:length]

    def transform_polynomial(f,n,T):

        # S = FractionField(f.parent())

        a = T[0]
        b = T[1]
        c = T[2]
        d = T[3]

        return ((c*x+d)^n*f.substitute({x: (a*x+b)/(c*x+d)})).numerator()

    def compute_transformation(f1, f2, n, M):

        # ADD SOME CHECKS?

        output = []

        if M.left_kernel().dimension() == 1:
            # print(transform_polynomial(f1,n,M.left_kernel().basis()[0]))
            f3 = transform_polynomial(f1,n,M.left_kernel().basis()[0])
            if f2*f3.lc() == f3*f2.lc():
                output.append(M.left_kernel().basis()[0])

        return output




    # HERE STARTS THE ACTUAL FUNCTION


    # We assume that both degree n homogeneous polynomials defined by f1 and f2
    # have reduced zero locus on the projective line. Therefore, their factorizations
    # should have matching degrees. We now perform these checks.

    # First we check the point at infinity. Since the point at infinity is a rational
    # point, the multiplicity of the root at this point must be at most one.
    # This leads to the following constraint on the degrees of the input polynomials

    # Only working with polynomials of degree at least 3 (EXPLAIN)
    assert n >= 3, "The degree of the forms must be at least 3"
    assert f1.degree() >= 3 and f2.degree() >= 3, "Input polynomials must be of degree at least 3"


    # Check degrees of the input polynomials
    assert f1.degree() in [n-1,n], "Degree of the first input polynomial is not in {n-1, n}"
    assert f2.degree() in [n-1,n], "Degree of the second input polynomial is not in {n-1, n}"


    # Storing factorisations of f1 and f2 for later use
    # NOTE: f.factor() turns the polynomial into a monic one.
    Fact1 = f1.factor()
    if f1 == f2:
        Fact2 = Fact1  # Maybe we don't really care.
    else:
        Fact2 = f2.factor()

    # Checking that f1 and f2 are square free
    assert all(item[1] == 1 for item in Fact1), "The first input polynomial is not squarefree"
    assert all(item[1] == 1 for item in Fact2), "The second input polynomial is not squarefree"

    # Now we are ready to check factorizations taking into account the points at
    # infinity, if necessary, i.e. in the case when f1 and/or f2 have degree n-1.
    degrees_1 = [item[0].degree() for item in Fact1]
    degrees_2 = [item[0].degree() for item in Fact2]

    # Just in case checking that Fact1 and Fact2 are sorted according to the degree
    # (lowest degree first)
    assert all(a <= b for a, b in zip(degrees_1, degrees_1[1:]))
    assert all(a <= b for a, b in zip(degrees_2, degrees_2[1:]))

    # if deg(f1) = n-1 or deg(f2), we insert one extra factor that corresponds to
    # the point at infinity
    if f1.degree() == n-1:
        degrees_1.insert(0, 1)

    if f2.degree() == n-1:
        degrees_2.insert(0, 1)

    print("degrees of factors for f1 and f2 are", degrees_1, "and", degrees_2)

    # Checking that all factors have the same degrees.
    # (assumes that degrees are sorted in increasing order, we asserted this above)
    if degrees_1 != degrees_2:
        print("Factors of f1 and f2 have differnt degrees")
        return False

    # Here is the starting of the actual algorithm
    else:

        # Here we are going to compute the possible GL2 transformations between f1 and f2.

        output = []

        # Notation for a/the highest degree factor in f1 and its degree
        f = Fact1[-1][0]
        dmax = f.degree()
        # denoting by rf the root of f
        L.<rf> = FF.extension(f)

        # The algorithm distinguishes 3 cases: dmax >= 3, dmax = 2, and dmax = 1.

        if dmax >= 3:
            # looping over the factors of f2 of degree dmax
            for g in [item[0] for item in Fact2 if item[0].degree() == dmax]:
                # looping over the roots of g in L (all roots of g are automatically in L)
                # denoting by rg the roots of g
                for rg in [item[0] for item in g.roots(L)]:
                    M = matrix(FF, 4, dmax)
                    rows = [rg, L.one(), -rf*rg, -rf]
                    for i in range(4):
                        M[i] = adding_zeros(rows[i].polynomial().list(), dmax)

                    output.extend(compute_transformation(f1,f2,n,M))

        elif dmax == 2:
            h = Fact1[-2][0] # need a better name (=second to last factor)
            if h.degree() == 2:
                rh = h.roots(L)[0][0] # both f and h split in L

                # looping over pairs (g1,g2) of distinct factors of f2 of degree 2
                for g1, g2 in permutations([item[0] for item in Fact2 if item[0].degree() == dmax], int(2)):
                    # looping over pairs (rg1,rg2) of roots of these factors
                    for rg1, rg2 in product([item[0] for item in g1.roots(L)], [item[0] for item in g2.roots(L)]):
                        M = matrix(FF, 4, 4)
                        rows_f = [rg1, L.one(), -rf*rg1, -rf]
                        rows_h = [rg2, L.one(), -rh*rg2, -rh]
                        for i in range(4):
                            M[i] = adding_zeros(rows_f[i].polynomial().list(), dmax) + adding_zeros(rows_h[i].polynomial().list(), dmax)

                        output.extend(compute_transformation(f1,f2,n,M))


            else:
                assert h.degree() == 1
                rh = h.roots(L)[0][0] # both f and h split in L

                # Factor of degree 2
                # g1 = Fact2[-1][0]

                # Roots of the unique quadratic factor of f2
                roots_quadratic_f2 = [item[0] for item in Fact2[-1][0].roots(L)]

                # roots of linear factors of f2 (if degree(f2) = n-1, then we add infinity)
                roots_linear_f2 = [Fact2[i][0].roots()[0][0] for i in range(len(Fact2)-1)]
                if f2.degree() == n-1:
                    roots_linear_f2.append("infty")

                # Loop over roots_quadratic_f2 and roots_linear_f2
                for rg1, rg2 in product(roots_quadratic_f2, roots_linear_f2):

                    rows_f = [rg1, L.one(), -rf*rg1, -rf]

                    if rg2 == "infty":
                        rows_h = [L.one(), L.zero(), -rh, L.zero()]
                    else:
                        rows_h = [rg2, L.one(), -rh*rg2, -rh]

                    M = matrix(FF, 4, 3)
                    for i in range(4):
                        M[i] = adding_zeros(rows_f[i].polynomial().list(), dmax) + adding_zeros(rows_h[i].polynomial().list(), 1)

                    output.extend(compute_transformation(f1,f2,n,M))


        # This is the case when f1 and f2 have only linear factors
        else:

            # Let us fix the "last" 3 roots of f1 (all finite)
            r1 = Fact1[-1][0].roots()[0][0]
            r2 = Fact1[-2][0].roots()[0][0]
            r3 = Fact1[-3][0].roots()[0][0]

            triple_of_roots_of_f1 = [r1,r2,r3]

            roots_of_f2 = [item[0].roots()[0][0] for item in Fact2]
            if f2.degree() == n-1:
                roots_of_f2.append("infty")

            for rg1, rg2, rg3 in permutations(roots_of_f2, int(3)):
                M = matrix(FF, 4, 3)

                if rg1 == "infty":
                    rows_h1 = [L.one(), L.zero(), -r1, L.zero()]
                else:
                    rows_h1 = [rg1, L.one(), -r1*rg1, -r1]

                if rg2 == "infty":
                    rows_h2 = [L.one(), L.zero(), -r2, L.zero()]
                else:
                    rows_h2 = [rg2, L.one(), -r2*rg2, -r2]

                if rg3 == "infty":
                    rows_h3 = [L.one(), L.zero(), -r3, L.zero()]
                else:
                    rows_h3 = [rg3, L.one(), -r3*rg3, -r3]

                for i in range(4):
                    M[i] = adding_zeros(rows_h1[i].polynomial().list(), 1) + adding_zeros(rows_h2[i].polynomial().list(), 1) + adding_zeros(rows_h3[i].polynomial().list(), 1)

                output.extend(compute_transformation(f1,f2,n,M))

    return output
