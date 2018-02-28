# factoring

This is a java project which factors integers for small inputs.
It is optimized for integers below 48 bits.
It provides different implementations, junit tests for correctness and measuring performance.
The best implementation uses the Lehman factorization, which runs in O(n^1/3).
Since the squares were stored in an array the main work which is the range checks were fast.
