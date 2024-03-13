reset()

load('IsGL2Equivalent.sage')

# degrees = [1, 1, 1, 1, 2]
R.<x> = GF(101)[]
f1 = 95*x^5 + 88*x^4 + 8*x^3 + 8*x^2 + 95*x + 1
f2 = 6*x^5 + 80*x^4 + 64*x^3 + 14*x^2 + 7*x + 17
print("Output of the function:", equivalence(f1,f2,6))

# # degrees = [2, 2, 2]
# R.<x> = GF(11)[]
# f1 = x^6 + 1
# f2 = 2*x^6 - 3*x^4 - 3*x^2 + 2
# print("Output of the function:", equivalence(f1,f2,6))

# # degrees = [1, 2, 3]
# R.<x> = GF(11)[]
# f1 = x^6 + x^4 + x^3 + x^2 + x + 1
# f2 = 2*x^6 + 4*x^5 + 4*x^4 - 3*x^2 - 4*x - 5
# print("Output of the function:", equivalence(f1,f2,6))


# # degrees = [1, 1, 1, 1, 2]
# R.<x> = GF(11)[]
# f1 = x^6 + 9*x^5 + 4*x^4 + 3*x^3 + 9*x^2 + 4*x + 2
# f2 = -2*x^6 - 2*x^5 + x^3 - 5*x^2 - 2*x - 1
# print("Output of the function:", equivalence(f1,f2,6))


# # degrees = [1, 1, 1, 1, 1, 1]
# R.<x> = GF(13)[]
# f1 = x^6 + 1
# f2 = 2*x^6 + 4*x^4 + 4*x^2 + 2
# print("Output of the function:", equivalence(f1,f2,6))


# # degrees = [1, 1, 1, 1, 1]
# R.<x> = GF(17)[]
# f1 = x^4 + 1
# f2 = x^5 + x
# print("Output of the function:", equivalence(f1,f2,5))


# # degrees = [1, 1, 1, 1, 1, 2]
# R.<x> = GF(11)[]
# f1 = (x^5+1)*(x^2 + 1)
# f2 = 9*x^6 + 5*x^4 + 4
# print("Output of the function:", equivalence(f1,f2,7))
