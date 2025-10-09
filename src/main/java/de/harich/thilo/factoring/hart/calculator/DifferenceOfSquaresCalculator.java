package de.harich.thilo.factoring.hart.calculator;

/**
 * provides operations on the Equation x^2 - y^2 = 4 + multiplier * n
 */
public interface DifferenceOfSquaresCalculator {


    /**
     * Has to calculate ceil (sqrt (4N * multiplier * multiplierIndex))
     * such that x^2 >  4 + multiplier * n = 4 + getMultiplier(multiplierIndex) * n
     * @param sqrt4N
     * @param multiplierIndex
     * @return
     */
    long getFirstX(double sqrt4N, int multiplierIndex);

    /**
     * Finds the next (bigger) x which might fulfill the Equation x^2 - 4 + multiplier * n = y^2
     *
     * @param x
     * @param multiplier
     * @param n
     * @return
     */
    public long adjustX(long x, long multiplier, long n);


    long getMultiplier(int index);

    /**
     * When we have not proper handling of lower numbers and squares in findSmallFactor
     * The loop
     * @param numberToFactorize
     * @param e
     */
    void handleRuntimeException(long numberToFactorize, RuntimeException e);

    /**
     * Since the Hart Algorithm is only working fast on numbers with factors > numberToFactorize^(1/3)
     * Handling of lower numbers and squares has to be handed in this method.
     * @param numberToFactorize
     * @return
     */
    long findSmallFactor(long numberToFactorize);

    long handleSolution(long greatestCommonDivisor, int multiplierIndex, long numberToFactorize);

    boolean hasNextMultiplier(int index);
}
