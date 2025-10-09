package de.harich.thilo.factoring.hart.calculator;

import de.harich.thilo.factoring.hart.calculator.aspects.VectorCalculator;
import de.harich.thilo.factoring.hart.multiplier.MultiplierFactory;

import java.util.Arrays;

/**
 * The adjustment of the x values alone makes the algorithm around 250% Faster on Celeron N5100
 * compared with the VectorCalculator.
 * Compared to the NativeDatatypeCalculator around 9 times faster (with the same multiplier) and
 * around 22 times faster when running NativeDatatypeCalculator with a multiplier of 1
 *
 * The adjustment has branches mod 8. But as long as the multipliers are generated in a linear way like:
 * multiplier = index * MULTIPLIER JIT can make it branch free, by unrolling the code.
 * This means on the other side using a list of good multipliers with no structure (mod 8) can hardly be optimized.
 */
public class AdjustXModPow2Calculator extends VectorCalculator {

    public static final int MODULUS_FOR_ADJUST_X = 8;

    protected long[] multiplier;

    public AdjustXModPow2Calculator() {
        super(false);

        int limit43Bit = 100000;

//        long[] smoothMultipliers = MultiplierFactory.createMultiplesOfPrimesBelow11(limit43Bit, MultiplierFactory::filterBigPrimesAndMultipleOf9);
        multiplier = MultiplierFactory.createMultiplesOf315(limit43Bit, MultiplierFactory::filterBigPrimesAndMultipleOf9);
//        long[] smoothMultipliers =  MultiplierFactory.createMultiplesOf315(limit43Bit, MultiplierFactory::filterBigPrimesAndMultipleOf9);
//        optimizeArraysMod8(smoothMultipliers);
    }

    @Override
    public long getMultiplier(int index) {
        return multiplier[index];
    }

    /**
     * Increases x to return the next possible solution on x for x^2 - 4 * multiplier * n = y^2 mod 64
     * Due to performance reasons we give back solutions for this equations modulo a
     * power of 2, since we can determine the solutions just by additions and binary
     * operations.
     * The branches depend on the value multiplier * n mod 8.
     * Since n is fixed, the branches only depend on multiplier mod 8
     * -> JIT works best if multiplier mod 8 is fixed, or we iterate over multiplier mod 8.
     * <p>
     * if multiplier is even x must be odd.
     * if multiplier*n == 3 mod 4 -> x = multiplier*n+1 mod 8
     * if multiplier*n == 1 mod 8 -> x = multiplier*n+1 mod 16 or -multiplier*n+1 mod 16
     * if multiplier*n == 5 mod 8 -> x = multiplier*n+1 mod 32 or -multiplier*n+1 mod 32
     *
     * @param x
     * @param multiplier
     * @param n
     * @return the improved x
     */
    @Override
    public long adjustX(long x, long multiplier, long n) {
        if ((multiplier&1)==0)
            return x | 1;
        final long kNPlus1 = multiplier*n+1;
        if ((kNPlus1 & 3) == 0)
        {
            return x + ((kNPlus1 - x) & 7);
        }
        else if ((kNPlus1 & 7) == 2) {
            final long adjust1 = ( kNPlus1 - x) & 15;
            final long adjust2 = (-kNPlus1 - x) & 15;

            return x + Math.min(adjust1, adjust2);
        }
        final long adjust1 = ( kNPlus1 - x) & 31;
        final long adjust2 = (-kNPlus1 - x) & 31;

        return x + (Math.min(adjust1, adjust2));
    }

    /**
     * After processing the multipliers, the multiplier (and sqrtMultiplier) array are structured in such a way,
     * that when iterating over it, the remainders of the multiplier are increasing by one modulus 8.
     * This helps to optimise the main loop.
     */
    protected void optimizeArraysMod8(long[] multipliers) {
        int length = (int) multipliers[0];
        long[][] multiplierPerRemainder = new long[MODULUS_FOR_ADJUST_X][2* (length / MODULUS_FOR_ADJUST_X)];
        int[] indexPerRemainder = fillRemainderLists(multipliers, length, multiplierPerRemainder);
        int minIndex = Arrays.stream(indexPerRemainder).min().getAsInt();


        // sorting each remainder list is not helping, but why?
//        for(int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++){
//            while (indexPerRemainder[remainder] < multiplierPerRemainder[remainder].length){
//                multiplierPerRemainder[remainder][indexPerRemainder[remainder]++] = Long.MAX_VALUE;
//            }
//            Arrays.sort(multiplierPerRemainder[remainder]);
//        }
        multiplier = new long[MODULUS_FOR_ADJUST_X * minIndex];
        sqrtMultiplier = new double[MODULUS_FOR_ADJUST_X * minIndex];

        for (int i= 0; i < minIndex; i++){
            for(int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++){
                int index = (i * MODULUS_FOR_ADJUST_X) + remainder;
                multiplier[index] = multiplierPerRemainder[remainder][i];
                sqrtMultiplier[index] = Math.sqrt(multiplier[index]);
            }
        }
    }

    static int[] fillRemainderLists(long[] multipliers, int length, long[][] multiplierPerRemainder) {
        int[] indexPerRemainder = {0,0,0,0,0,0,0,0};
        for (int index = 1; index < length; index++ ){
            long goodMultiplier = multipliers[index];
            int remainder = (int) (goodMultiplier % MODULUS_FOR_ADJUST_X);
            multiplierPerRemainder[remainder][indexPerRemainder[remainder]++] = goodMultiplier;
        }
        return indexPerRemainder;
    }
}
