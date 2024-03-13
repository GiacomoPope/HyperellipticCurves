Context from Wouter

```
Maxim has finished the equivalence of binary forms (already over a month ago...), and all that remains to be done is to incorporate this into a function for finding isomorphisms between hyperelliptic curves, which is really not a big task: if

   H : y^2 = f(x)  and  H' : y^2 = f'(x),

then Maxim's function finds all matrices [a, b \\ c, d] such that (c*x + d)^6*f((a*x + b)/(c*x + d)) = lambda * f'(x) for some lambda. If lambda is a square, then there is a corresponding isomorphism, amounting to substituting sqrt(lambda)*y/(c*x + d)^3 for y and (a*x + b)/(c*x + d) for x. (If lambda is not a square, then it takes you to the quadratic twist.)
```

Context from Maxim

```
The file IsGL2Equivalent.sage contains the implementation of our function, which for now is simply called equivalence. We need a better name.

The second file examples.sage contains examples. It loads the first file and then computes an example. There are several examples there, you just need to comment out the current oue and uncomment the one you would like to run. For each example you see there also the degrees of factors of the polynomials, so that you can choose your example.

You can run this by starting sage in the terminal from the directory containing both files and then typing:

%runfile examples.sage
```
