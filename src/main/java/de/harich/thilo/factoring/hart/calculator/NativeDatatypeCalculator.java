package de.harich.thilo.factoring.hart.calculator;

/**
 * This is the really basic version
 */
public class NativeDatatypeCalculator
        implements DifferenceOfSquaresCalculator {

    // by multiplying a number n with m for x^2 - m*n = y^2 mod m for any x, meaning there are many solutions,
    // but the number is growing
    // the multiplier  3*3*5*7 increases performance by a factor of 2.5
    private int multiplier = 3*3*5*7;

    public NativeDatatypeCalculator(int multiplier) {
        this.multiplier = multiplier;
    }


    public long getFirstX(double sqrt4N, int i) {
        return (long) Math.ceil(sqrt4N * Math.sqrt(multiplier * i));
    }


    @Override
    public long adjustX(long x, long multiplier, long n) {
        return x;
    }


    @Override
    public long getMultiplier(int index) {
        return (long) index * multiplier;
    }

    @Override
    public void handleRuntimeException(long numberToFactorize, RuntimeException e) {
        System.out.println("Failed to factor N =" + numberToFactorize +
                ". Either it has factors < cbrt (N) needing trial division, or the number is a square");
    }

    @Override
    public long findSmallFactor(long numberToFactorize) {
        if (numberToFactorize < 4)
            return numberToFactorize;
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
