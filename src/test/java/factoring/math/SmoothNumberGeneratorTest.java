package factoring.math;

import factoring.sieve.triple.BitSetWithOffset;
import org.junit.Test;

import static java.lang.Math.*;
import static org.junit.Assert.assertTrue;

public class SmoothNumberGeneratorTest {


    @Test
    public void testIsSmoothByFourFactors(){
//        int number = 87463;
//        long number = 3363 * 2287;
        final long number = 8_000_0009l * 4_000_037l;
        SmoothNumberGenerator smoothNumberGenerator = SmoothNumberGenerator.ofMaxSize(2 * number);
        int smoothNumberCount = 0;
        int smoothRightCount = 0;

        int numberSqrt = (int) ceil(sqrt(number));
        final double factorBaseSize = smoothNumberGenerator.factorBaseSize(numberSqrt);
        System.out.println("factor base size : " + factorBaseSize);
        int rangeOrig = (int) (3* factorBaseSize);
        System.out.println("work per loop: " + rangeOrig / 64);

        double numberPowQuarter = pow(number, .25);
        final int variance = (int) (3 * numberPowQuarter * rangeOrig);
        System.out.println(" expected size of numbers : " + variance + " = " + (variance/ numberSqrt) + " * size of QS");
        int zDiff = 0;
        int z = numberSqrt + zDiff;
        int range = rangeOrig;
        // z^2 - t
        for (; smoothNumberCount < factorBaseSize;) {
            BitSetWithOffset smoothNumbersSet = smoothNumberGenerator.getSmoothNumbersBy4Products(number, z, range);
            for (int yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(0); yDiff >= 0; yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(yDiff + 1)) {
                long y = smoothNumbersSet.yOffset + yDiff;
                long lowerProductSmooth = z - y;
                long higherProductSmooth = z + y;
                long smoothNumber = lowerProductSmooth * higherProductSmooth;
                assertTrue("lower number " + lowerProductSmooth + " is not smooth for number : " + number, smoothNumberGenerator.smoothNumber.get((int) lowerProductSmooth));
                assertTrue(smoothNumberGenerator.smoothNumber.get((int) higherProductSmooth));

                int right = (int) (smoothNumber - number);
                if (smoothNumberGenerator.smoothNumber.get(abs(right))) {
                    System.out.println("found smooth number right : " + right);
                    smoothRightCount++;
                }

                if (smoothNumber == 320007290814900l){
                    System.out.println();
                }
                assertTrue("range does not fit for number : " + smoothNumber + " size of diff to sqrt n " + right, abs(right) < numberSqrt);
                assertTrue("generated number " + smoothNumber + " distance to number " + (smoothNumber - number) + " greater than .5 * sqrt(n) ", abs(right) < 1 * numberSqrt);
                smoothNumberCount++;
            }
            range = (int) (rangeOrig / sqrt(++zDiff));
            z = numberSqrt + zDiff;
        }
        System.out.println("we have created " +smoothNumberCount + " numbers within " + zDiff + " loops. all smooth numbers smaller then .5 * sqrt(n)");
        System.out.println("we have created " +smoothRightCount + " smooth numbers on the right");
        System.out.println("last range " + range);

    }

    @Test
    public void testIsSmoothByTwoProducts(){
//        int number = 87463;
//        long number = 3363 * 2287;
        final long number = 8_000_0009l * 4_000_037l;
        SmoothNumberGenerator smoothNumberGenerator = SmoothNumberGenerator.ofMaxSize(2 * number);
        int smoothNumberCount = 0;
        int smoothRightCount = 0;

        int numberSqrt = (int) ceil(sqrt(number));
        final double factorBaseSize = smoothNumberGenerator.factorBaseSize(number);
        System.out.println("factor base size : " + factorBaseSize);
        int range = (int) (3* factorBaseSize);
        System.out.println("work per loop: " + range / 64);

        double numberPowQuarter = pow(number, .25);
        final int variance = (int) (3 * numberPowQuarter * range);
        System.out.println(" expected size of numbers : " + variance + " = " + (variance/ numberSqrt) + " * size of QS");
        int zDiff = 0;
        int z = numberSqrt + zDiff;
        int rangeForZ = range;
        for (; smoothNumberCount < factorBaseSize;) {
            BitSetWithOffset smoothNumbersSet = smoothNumberGenerator.getSmoothNumbersBy2Products(number, z, rangeForZ);
            for (int yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(0); yDiff >= 0; yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(yDiff + 1)) {
                long y = smoothNumbersSet.yOffset + yDiff;
                long lowerProductSmooth = z - y;
                long higherProductSmooth = z + y;
                long smoothNumber = lowerProductSmooth * higherProductSmooth;
                assertTrue("lower number " + lowerProductSmooth + " is not smooth for number : " + number, smoothNumberGenerator.smoothNumber.get((int) lowerProductSmooth));
                assertTrue(smoothNumberGenerator.smoothNumber.get((int) higherProductSmooth));

                int right = (int) (smoothNumber - number);
                if (smoothNumberGenerator.smoothNumber.get(abs(right))) {
                    System.out.println("found smooth number right : " + right);
                    smoothRightCount++;
                }

                if (smoothNumber == 320007290814900l){
                    System.out.println();
                }
                assertTrue("range does not fit for number : " + smoothNumber + " size of diff to sqrt n " + right, abs(right) < numberSqrt);
                assertTrue("generated number " + smoothNumber + " distance to number " + (smoothNumber - number) + " greater than .5 * sqrt(n) ", abs(right) < 1 * numberSqrt);
                smoothNumberCount++;
            }
            rangeForZ = (int) (range / sqrt(++zDiff));
            z = numberSqrt + zDiff;
        }
        System.out.println("we have created " +smoothNumberCount + " numbers within " + zDiff + " loops. all smooth numbers smaller then .5 * sqrt(n)");
        System.out.println("we have created " +smoothRightCount + " smooth numbers on the right");
        System.out.println("last range " + rangeForZ);

    }

    @Test
    public void testIsSmoothByMultiplyingN(){
//        int number = 87463;
//        long number = 3363 * 2287;
        final long numberOrig = 8_000_0009l * 4_000_037l;
        SmoothNumberGenerator smoothNumberGenerator = SmoothNumberGenerator.ofMaxSize(10 * numberOrig);
        int smoothNumberCount = 0;
        int smoothRightCount = 0;

        final double factorBaseSize = smoothNumberGenerator.factorBaseSize(numberOrig);
        System.out.println("factor base size : " + factorBaseSize);
        int range = (int) (2* factorBaseSize);
        System.out.println("work per loop: " + range / 64);

        double numberPowQuarter = pow(numberOrig, .25);
        int variance = (int) (3 * numberPowQuarter * range);
        int numberSqrt = (int) sqrt(numberOrig);
        System.out.println(" expected size of numbers : " + variance + " = " + (variance/ numberSqrt) + " * size sqrt(n) = " +numberSqrt);

        long number = 0;
        int k = 0;

        for (; smoothNumberCount < factorBaseSize;) {
            number += numberOrig;
            k++;
            System.out.println(" k = " + k);
            numberPowQuarter = pow(number, .25);
            variance = (int) (3 * numberPowQuarter * range);
            numberSqrt = (int) ceil(sqrt(number));
            BitSetWithOffset smoothNumbersSet = smoothNumberGenerator.getSmoothNumbersBy2Products(number, numberSqrt, range);
            for (int yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(0); yDiff >= 0; yDiff = smoothNumbersSet.ysMinyOffset.nextSetBit(yDiff + 1)) {
                long y = smoothNumbersSet.yOffset + yDiff;
                long lowerProductSmooth = numberSqrt - y;
                long higherProductSmooth = numberSqrt + y;
                long smoothNumber = lowerProductSmooth * higherProductSmooth;
                assertTrue("lower number " + lowerProductSmooth + " is not smooth for number : " + number, smoothNumberGenerator.smoothNumber.get((int) lowerProductSmooth));
                assertTrue(smoothNumberGenerator.smoothNumber.get((int) higherProductSmooth));

                int right = (int) (smoothNumber - number);
                if (smoothNumberGenerator.smoothNumber.get(abs(right))) {
                    System.out.println("found smooth number right : " + right);
                    smoothRightCount++;
                }

                if (smoothNumber == 320000931906359l){
                    System.out.println();
                }
                assertTrue("range does not fit for number : " + smoothNumber + " size of diff to sqrt n " + right, abs(right) < numberSqrt);
                assertTrue("generated number " + smoothNumber + " distance to number " + (smoothNumber - number) + " greater than .5 * sqrt(n) ", abs(right) < 1 * numberSqrt);
                smoothNumberCount++;
            }
        }
        System.out.println("we have created " +smoothNumberCount + " numbers within " + k + " loops. all smooth numbers smaller then .5 * sqrt(n)");
        System.out.println("we have created " +smoothRightCount + " smooth numbers on the right");

    }
}
