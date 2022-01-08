package factoring.sieve.triple;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import org.junit.Test;

import java.math.BigInteger;

import static java.lang.Math.ceil;
import static java.lang.Math.sqrt;
import static org.junit.Assert.assertTrue;

public class SmoothNumbersTest {

    FactorAlgorithm sieve = new QuadroLookupSieve(18);

    @Test
    public void test87463() {
        sieve.factor(BigInteger.valueOf(87463));
    }

    @Test
    public void test13_17() {
        sieve.factor(BigInteger.valueOf(13 * 17));
    }

    @Test
    public void test13_23() {
        sieve.factor(BigInteger.valueOf(29 * 23));
    }

    @Test
    public void test8462977() {
        sieve.factor(BigInteger.valueOf(8462977));
    }

    @Test
    public void test769311827() {
        sieve.factor(BigInteger.valueOf(33637 * 22871));
    }

    @Test
    public void test8000009_4000037() {
        final long n = 8_000_0009l * 4_000_037l;
        sieve.factor(BigInteger.valueOf(n));
    }

    @Test
    public void test28Bit() {
        //        final long n = 150000001l * 300000007l;
        final long n = 150000001l * 280000027l;
        sieve.factor(BigInteger.valueOf(n));
    }

    @Test
    public void testSmoothBitset()
    {
        QuadroLookupSieve quadroSieve = new QuadroLookupSieve(18);
        int number = 288;
        do {
            int z = (int) ceil(sqrt(number));
            final int yRange = (int) Math.pow(number, .25);
            double numberPowQuarter = Math.pow(number, .25);

            int numberSqrtCeil = (int) ceil(sqrt(number));
            BitSetWithOffset smoothAscendingFromLower = quadroSieve.getSmoothNumbersInRange(number, numberSqrtCeil, yRange);
            for (int yDiff = smoothAscendingFromLower.ysMinyOffset.nextSetBit(0); yDiff >= 0; yDiff = smoothAscendingFromLower.ysMinyOffset.nextSetBit(yDiff + 1)) {
                int y = smoothAscendingFromLower.yOffset + yDiff;
                int lowerProductSmooth = z - y;
                int higherProductSmooth = z + y;
                int smoothNumberResult = lowerProductSmooth * higherProductSmooth;

                assertTrue("lower number " + lowerProductSmooth + " is not smooth for number : " + number , quadroSieve.smoothNumber.get(lowerProductSmooth));
                assertTrue(quadroSieve.smoothNumber.get(higherProductSmooth));

                final double variance = 4 * numberPowQuarter * yRange;
                assertTrue("upper range does not fit for number : " + number, smoothNumberResult <= number + variance);
                assertTrue("lower range does not fit for number : " + number, smoothNumberResult >= number - variance);
            }
            number++;
        }while (number < 310);
     }

}
