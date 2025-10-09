package de.harich.thilo.factoring.hart;

import de.harich.thilo.factoring.hart.calculator.DifferenceOfSquaresCalculator;
import de.harich.thilo.math.BinaryGreatestCommonDivisorEngine;

/**
 * A clean code Approach for a Hart Factoring algorithm running in O(n^1/3),
 * if at least one Factor of size > n^(1/3) exist, and the number is no Square.
 * This is the Basic Algorithm. The details for good performance are in HartCalculator.
 */
public class HartFactorizationAlgorithm {

    // use decomposition not Inheritance,
    // the Algorithm itself is the same, only the way we calculate the values is different
    DifferenceOfSquaresCalculator calculator;
    BinaryGreatestCommonDivisorEngine greatestCommonDivisorEngine = new BinaryGreatestCommonDivisorEngine();

    public HartFactorizationAlgorithm(DifferenceOfSquaresCalculator calculator) {
        this.calculator = calculator;
    }

    public long findSingleFactor(long numberToFactorize) {
        // handle all easy cases
        if (numberToFactorize < 4)
            return numberToFactorize;
        long smallFactor = calculator.findSmallFactor(numberToFactorize);
        if (smallFactor > 1)
            return smallFactor;

        try {
            return findBigSingleFactor(numberToFactorize);
        } catch (final RuntimeException e) {
            calculator.handleRuntimeException(numberToFactorize, e);
            // no factor can be found
            return 1;
        }
    }

    private long findBigSingleFactor(long numberToFactorize) {
        final long fourN = 4 * numberToFactorize;
        final double sqrt4N = Math.sqrt(fourN);


        // this loop is pretty simple and has no branches except, when we find the solution,
        // this makes it a good candidate for parallel processing. In Java we can let the JIT compiler do the work
        // for us. We can still write clean Java code. In the SimdCalculator we have a simple array.
        // Intrinsic in Java or using AVX in C is not (much) faster, but really complex.

        // infinite loop with no range check speeds up by ~ 5%
//        for (int i = 1; i < ((Mod8GroupsMultiplierCalculator)calculator).analyzeUpToIndex; i++) {
//        for (int i=1; ;i++) {
        for (int i=1; calculator.hasNextMultiplier(i) ;i++) {

            // create a multiplier
            long multiplier = calculator.getMultiplier(i);

            // get an x = ceil sqrt (4 * n * multiplier)
            long x = calculator.getFirstX(sqrt4N, i);

            // get an x which might have a better chance to solve  y^2 = x^2 - 4 * multiplier * n
            x = calculator.adjustX(x, multiplier, numberToFactorize);

            // we look for y^2 = x^2 - 4 * multiplier * n
            // the number x^2 and multiplier * fourN might exceed 64 bits, but the difference is usually still correct
            // this is the reason we switch to long calculation here.
            long possibleSquare = x*x - multiplier * fourN;

            // the fastest way to check if a number is a square - if lower Bits are bits of a square - seems to be
            // to calculate the square root, which is usually only supported for double values. So another
            // conversion to double needs to be performed. Support for such conversions, or avoiding in AVX is complex.
            long y = (long) Math.sqrt(possibleSquare);
            if (y*y == possibleSquare) {
                long greatestCommonDivisor = greatestCommonDivisorEngine.greatestCommonDivisor(x + y, numberToFactorize);
                if (greatestCommonDivisor  > 1 && greatestCommonDivisor < numberToFactorize) {
//                    return greatestCommonDivisor;
                    long factor = calculator.handleSolution(greatestCommonDivisor, i, numberToFactorize);
                    if (factor > 1)
                        return factor;
                }
            }
        }
        return 1;
    }



    public DifferenceOfSquaresCalculator getCalculator() {
        return calculator;
    }

}
