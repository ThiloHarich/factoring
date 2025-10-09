package de.harich.thilo.factoring.hart.calculator.aspects;

import de.harich.thilo.factoring.hart.calculator.DifferenceOfSquaresCalculator;

/**
 * The arrays enable Vector (MMX, SIMD, AVX) optimizations.
 * This alone on a Celeron N5100 with 4 kernels and Intel® SSE4.2 on 64 bits is
 * around 3 times faster than the calculation on primitive datatypes as in NativeDatatypeCalculator
 */
public class VectorCalculator
        implements DifferenceOfSquaresCalculator {

    // Only use when CPU supports fusedMultiplyAdd, otherwise Algorithm will slow down by a factor of ~ 500.
    protected boolean useFusedMultipleAdd = false;

//    protected long baseMultiplier = 3*3*5*7;
    static final long baseMultiplier = 3*3*5*7;


    // 21 BIT should be enough for number around 21*3 BIT
    protected static final int MULTIPLIERS_LIMIT_50_BIT = 1 << 21;

    static final double ROUND_UP_DOUBLE = 0.9999999665;

    // SIMD/AVX operations work on arrays, here it is. Let the JIT do the optimization
    protected double [] sqrtMultiplier = new double[MULTIPLIERS_LIMIT_50_BIT];

    public VectorCalculator() {
        for (int i=1; i<MULTIPLIERS_LIMIT_50_BIT; i++) {
            sqrtMultiplier[i] = Math.sqrt(i* baseMultiplier);
        }
    }

    public VectorCalculator(boolean useFusedMultipleAdd) {
        this();
        this.useFusedMultipleAdd = useFusedMultipleAdd;
    }

    @Override
    public long adjustX(long x, long multiplier, long n) {
        return x;
    }

    public long getFirstX(double sqrt4N, int i) {
        // we want to calculate ceil (sqrt (4N * multiplier * i));
        // calculating the sqrt here can be much slower than using the array on SIMD/AVX architectures.
        // we take a double Array sqrtMultiplier and convert to long (array?), it is not clear how JIT can do this
        // but since we get some reasonable speed advantage it can do it.
        // The rounding trick seems to not work for every CPU:
        if (useFusedMultipleAdd)
//            // instead of return (long) (sqrt4N * sqrtMultiplier[i] + ROUND_UP_DOUBLE);
//            // we use fusedMultiplyAdd is supported by AVX, so we use it together with a rounding trick
//            // Might save around 4 %
            return (long) Math.fma(sqrt4N, sqrtMultiplier[i], ROUND_UP_DOUBLE);
        return (long) Math.ceil(sqrt4N * sqrtMultiplier[i]);
    }



    @Override
    public long getMultiplier(int index) {
        return index * baseMultiplier;
    }

    @Override
    public void handleRuntimeException(long numberToFactorize, RuntimeException e) {
        System.out.println("Failed to factor N =" + numberToFactorize +
                ". " + numberToFactorize + " might have factors < below " + Math.ceil(Math.cbrt(numberToFactorize)) +
                " (n^1/3) -> add handling in findSmallFactor, the array for the stored square roots is too small " +
                "-> increase value of MULTIPLIERS_LIMIT_50_BIT, or the number is a square -> add handling in findSmallFactor");
        System.out.println("Exception is : " + e.getLocalizedMessage());
    }

    @Override
    public long findSmallFactor(long numberToFactorize) {
        // TODO this makes the algorithm slow, but why?
//        if (BigInteger.valueOf(numberToFactorize).isProbablePrime(10))
//            return numberToFactorize;
//        if (numberToFactorize < 4)
//            return numberToFactorize;
//        long sqrt = (long) Math.sqrt(numberToFactorize);
//        if (sqrt * sqrt == numberToFactorize)
//            return sqrt;
        return 1L;
    }

    @Override
    public long handleSolution(long greatestCommonDivisor, int multiplierIndex, long numberToFactorize) {
        return greatestCommonDivisor;
    }

    @Override
    public boolean hasNextMultiplier(int index) {
        return true;
    }

}
