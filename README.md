# factoring

This is a java project which factors integer values for small inputs.
It is optimized for integers below 48 bits.
It provides different implementations, junit tests for correctness and measuring performance.
The best implementation uses the Lehman factorization, which runs in O(n^1/3).
The main performance improvements is removing the expensive range checks in the inner loop for each multiplier k of the number to be factored. The expensive calculations of the square roots for the lower bound of the range were done only once and then reused. The check of the upper bound (calculation the sixth root) was removed completely in 98,5 % of the cases.

Unfortunately the project is not structured well at the moment. A lot of different implementations - some working some do not - of different algorithms lie in the same package. Maybe the most interesting algorithms are the Lehman and the Trial division algorithms, because they were pertty fast. Up to my knowledge the Lehman implementation - done together with Tillman Neuman - is the fastest I know. At least it is much faster then the java version of the YAFU implementation of the Lehman algorithms.

The trial division algorithm uses the fact that multiplications of (double) numbers are much faster the division. So instead of doing a trial division it is doing a trial multiplication with the inverse of a prime. This results in a much faster trial division algorithm, which can also be used in the Lehman algorithm.

It is planned to tidy this up. So far I would recommend to use the framework from Tilman Neumann https://github.com/TilmanNeumann/java-math-library since it anyhow uses the nice algorithms of this project and combines the best of both projects.

