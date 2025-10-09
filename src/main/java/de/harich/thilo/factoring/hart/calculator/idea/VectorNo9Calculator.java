package de.harich.thilo.factoring.hart.calculator.idea;

/**
 * Introduce a second array for the multipliers. This gives us more freedom to do generate more complex multipliers.
 * We can avoid to create multipliers which create fewer solutions as other multipliers.
 * We leave out another 9 works when we have no dependency of the loop execution on the generated numbers,
 * this will be faster, but if we do some operations depending on the generated multiplier,
 * like optimizing x in AdjustXModPow2Calculator, the JIT can no longer do its work properly, and the algorithm
 * will slow down, so this is not working together with the optimization on the x values.
 *
 */
public class VectorNo9Calculator
        extends Vector2ArraysCalculator {


    // SIMD/AVX operations work on arrays, here it is. Let the JIT do the optimization
    final long [] multiplier = new long[MULTIPLIERS_LIMIT_50_BIT];
//    sqrtMultiplier = new double[][MULTIPLIERS_LIMIT_50_BIT];

    int index = 1;
    public VectorNo9Calculator() {
        for (int i=1; i<MULTIPLIERS_LIMIT_50_BIT; i++) {
            if (i % 9 != 0)
            {
                multiplier[index] = i * this.baseMultiplier;
                sqrtMultiplier[index++] = Math.sqrt(i * this.baseMultiplier);
            }
        }

    }

    @Override
    public long getMultiplier(int index) {
        return multiplier[index];
    }
}
