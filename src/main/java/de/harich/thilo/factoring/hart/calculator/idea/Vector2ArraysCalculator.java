package de.harich.thilo.factoring.hart.calculator.idea;

import de.harich.thilo.factoring.hart.calculator.AdjustXModPow2Calculator;

/**
 * Introduce a second array for the multipliers. This gives us more freedom to do generate more complex multipliers.
 * We can avoid to create multipliers which create fewer solutions as other multipliers.
 * It might lose ~ 2% of the SimdCalculator
 */
public class Vector2ArraysCalculator
        extends AdjustXModPow2Calculator {
//        extends SimdCalculator {

    static final long baseMultiplier = 3*3*5*7;

    // SIMD/AVX operations work on arrays, here it is. Let the JIT do the optimization
    final long [] multiplier = new long[MULTIPLIERS_LIMIT_50_BIT];
//    sqrtMultiplier = new double[][MULTIPLIERS_LIMIT_50_BIT];

    int index = 1;
    public Vector2ArraysCalculator() {
        for (int i=1; i<MULTIPLIERS_LIMIT_50_BIT; i++) {
                multiplier[index] = i * baseMultiplier;
                sqrtMultiplier[index++] = Math.sqrt(i * baseMultiplier);
        }

    }

    @Override
    public long getMultiplier(int index) {
        return multiplier[index];
    }
}
