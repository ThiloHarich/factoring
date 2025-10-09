package de.harich.thilo.factoring.hart.calculator.idea;

import de.harich.thilo.factoring.hart.calculator.AdjustXModPow2Calculator;

/**
 * For each of the 8 remainders of n * multiplier mod 8 we have a special array collecting the solutions.
 * Additionally we collect some of them in a group.
 * The idea is to always have same case in one Vector. This might help the JIT to always have the same instruction
 * within the different registers in the Vector. At least for SIMD this seems not to increase the speed, at least
 * for i * baseMultiplier.
 */
public class RemainderGroupArraysCalculator
        extends AdjustXModPow2Calculator {
//        extends SimdCalculator {

    static final long baseMultiplier = 3*3*5*7;

    // SIMD/AVX operations work on arrays, here it is. Let the JIT do the optimization
    final long [] multiplier = new long[MULTIPLIERS_LIMIT_50_BIT];
//    sqrtMultiplier = new double[][MULTIPLIERS_LIMIT_50_BIT];

//    int index = 1;
    int numGroups = 8;
    int numCasesAdjust = 8;
    int[] indexForRemainder = new int [] {0,1,2,3,4,5,6,7};
    public RemainderGroupArraysCalculator() {
        for (int i=0; i<MULTIPLIERS_LIMIT_50_BIT - numGroups * numCasesAdjust;) {
            for (int remainder = 0; remainder < numCasesAdjust; remainder++) {
//                if (remainder != 3)
//                for (int group = 0; group < (remainder == 3 ? numGroups*2 : numGroups) ; group++) {
                for (int group = 0; group < numGroups ; group++) {
                    long currentMultiplier = indexForRemainder[remainder] * baseMultiplier;
                    if (indexForRemainder[remainder] % numCasesAdjust != remainder) {
                        System.out.println("Fehler");
                    }
                    indexForRemainder[remainder] += 8;
                    multiplier[i] = currentMultiplier;
                    sqrtMultiplier[i++] = Math.sqrt(currentMultiplier);
                }
            }
        }
    }

    @Override
    public long getMultiplier(int index) {
        return multiplier[index];
    }
}
