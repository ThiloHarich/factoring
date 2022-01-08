package factoring.math;

import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.SortedMultiset_BottomUp;
import factoring.sieve.triple.BitSetWithOffset;
import factoring.trial.TDiv31Barrett;

import java.math.BigInteger;
import java.util.BitSet;

import static java.lang.Math.*;
import static java.lang.Math.log;

public class SmoothNumberGenerator {

    int smoothBound;
    BitSet smoothNumber;
    BitSet smoothNumberReverse;
    private int smoothBits;
    TDiv31Barrett smallFactoriser = new TDiv31Barrett();


    public static SmoothNumberGenerator ofMaxSize(long maxSize){
        return new SmoothNumberGenerator(maxSize);
    }

    public SmoothNumberGenerator(long maxSize) {
        this.smoothBound = (int) sqrt(2*maxSize);
        initSmoothNumbers();
    }

    public double factorBaseSize(long n) {
        return 1 * sqrt(exp(sqrt(log(n) * log(log(n)))));
    }

    /**
     *
     * @param number
     * @param numberSqrt
     * @param yRange
     * @return
     */
    BitSetWithOffset getSmoothNumbersBy4Products(long number, int numberSqrt, int yRange) {
        // if number is smaller than lookup table just return the Bitset
        // do a loop with some higher numbers
        long z = numberSqrt;

        int diffNumberToSquare = (int) (z * z - number);
        double yOptDouble = sqrt(diffNumberToSquare);   // <= sqrt(2*sqrt(number)) = sqrt(2)*number^1/4
        // we try to write number = z^2 - y^2
        int yOptInt = (int) round(yOptDouble);
        // (y + yRange)^2 = y^2 + 2y * yRange + yRange^2 ~  2 * sqrt(2)*number^1/4 * yRange < 3 * *number^1/4 * yRange
        int yRangeAdjust = yRange > yOptInt ? yOptInt : yRange;
        BitSet ysMinyOffset = new BitSet();

        // we search ner sqrt(number). Around a number2 were we know that a product of non-smooth numbers is near
        // the original number.
        long number2 = z - yOptInt;
//        for (int number2 = z - yOptInt - yRange; number2 < z - yOptInt + yRange; number2++) {
            int numberSqrt2 = (int) ceil(sqrt(number2));
            final double factorBaseSize = factorBaseSize(number2);
//            int yRange2 =  (int) (2* factorBaseSize);
            BitSetWithOffset smoothAscendingFromLower = getSmoothNumbersBy2Products(number2, numberSqrt2, yRange);
            for (int yDiff = smoothAscendingFromLower.ysMinyOffset.nextSetBit(0); yDiff >= 0; yDiff = smoothAscendingFromLower.ysMinyOffset.nextSetBit(yDiff + 1)) {
                // operate on index i here
                if (yDiff == Integer.MAX_VALUE) {
                    break; // or (i+1) would overflow
                }
                int y = smoothAscendingFromLower.yOffset + yDiff;

                long smoothLowerNumber = (numberSqrt2 - y) * (numberSqrt2 + y);
                long complementNumber = (int) (numberSqrt + (numberSqrt - smoothLowerNumber));
                int nDivLower = (int) (number / smoothLowerNumber);
                if (complementNumber != nDivLower)
                    System.out.println("number to far away from sqrt?");
                int[] rep = potentiallySmoothFermat(nDivLower);

                if (rep != null) {
                    System.out.println("found smooth (x-t) " + smoothLowerNumber + " and (x+t) " + complementNumber);
                    int rest = (int) (smoothLowerNumber * complementNumber - number2);
                    System.out.println("smooth number  : " + smoothLowerNumber * complementNumber + " and rest " + rest);
                }
//            }
        }
        int yOffset = (int) (z + yOptInt - yRange);
        BitSetWithOffset result = new BitSetWithOffset(ysMinyOffset, yOffset);
        return result;
    }

    /**
     * we want number = z^2 - y^2 = (z - y) * (z + y)
     * z^2 = {0,1,4} - number mod 8
     * With a mod 2^k argument we find some valid z'
     * if we have some z
     *
     * @return
     */
    int[] potentiallySmoothFermat(int number) {
        // do some mod arguments here. reuse Smoothumbers class!?
        int numberMod8 = number % 8;
        int higherSquareSqrt = (int) ceil(sqrt(number));
        int z = higherSquareSqrt;

        int zRange = (int) exp(sqrt(log(higherSquareSqrt) * log(log(higherSquareSqrt))));
        do {
            int zMod8 = z % 8;
            int zSquaredMod8 = (z * z) % 8;

            int diffNumberToSquare = z * z - number;
            double yOptDouble = sqrt(diffNumberToSquare);   // <= sqrt(2*sqrt(number)) = sqrt(2)*number^1/4
            // we try to write number = x^2 - y^2
            int yOptInt = (int) round(yOptDouble);
            if (yOptInt * yOptInt == diffNumberToSquare) {
                int lowerNumber = z - yOptInt;
                int higherNumber = z + yOptInt;
                if (smoothNumber.get(lowerNumber)) {
                    if (smoothNumber.get(higherNumber))
                        return new int[]{lowerNumber, higherNumber};
                }
            }
            z++;
        } while (z < higherSquareSqrt + zRange);
        return null;
    }
    /**
     * given a number we will return a BitSet of numbers y' and a yOffset
     * stored in BitSetWithOffset
     * with y = y' + yOffset
     *
     * z^2 - y^2 = (z - y) * (z + y) = number +/- 3 * number^1/4 * (sqrt(numberSqrt - sqrt(number)) * yRange
     * with
     * (z - y) and (z + y) being smooth
     * z = ceil(sqrt(number)) + l
     * ((z+l) + (y+k)) * ((z+l) z - (y+k)) = (z+l)^2 - (y+k)^2) = z^2 + 2zl + l^2 - y^2 -2ky - k^2
     * = z^2 + 2zl + l^2 - y^2 -2ky - k^2
     * y^2 = 2zl -> y = n^1/4* l^1/2
     *
     * @param number
     * @return
     */
    // TODO provide an Iterator here
    BitSetWithOffset getSmoothNumbersBy2Products(long number, int numberSqrt, int yRange) {
        // if number is smaller than lookup table just return the Bitset
        // do a loop with some higher numbers
        int z = numberSqrt;

        int diffNumberToSquare = (int) (z * z - number);
        double yOptDouble = sqrt(diffNumberToSquare);   // <= sqrt(2*sqrt(number)) = sqrt(2)*number^1/4
        // we try to write number = x^2 - y^2
        int yOptInt = (int) round(yOptDouble);
        // (y + yRange)^2 = y^2 + 2y * yRange + yRange^2 ~  2 * sqrt(2)*number^1/4 * yRange < 3 * *number^1/4 * yRange
        int yRangeAdjust = yRange > yOptInt ? yOptInt : yRange;
        int yBegin = yOptInt - yRangeAdjust;
        // Z - Y_opt + range -> Z - y_opt  + range
        int lowerProduct = z - yBegin;
        // z + y_opt - range -> z + y_opt - range
        int higherProduct = z + yBegin;
//        int smoothRange = Math.max(64, 2 * yRange);
        BitSet smoothAscendingFromHigher = smoothNumber.get(higherProduct, higherProduct + 2*yRangeAdjust);
        BitSet smoothDescendingFromLower = smoothNumberReverse.get(smoothBound - lowerProduct, smoothBound - (lowerProduct - 2*yRangeAdjust));
        smoothAscendingFromHigher.and(smoothDescendingFromLower);
        return new BitSetWithOffset(smoothAscendingFromHigher, yBegin);
    }

    private void initSmoothNumbers() {
        long start = System.currentTimeMillis();


        smoothNumber = new BitSet(smoothBound);
        smoothNumberReverse = new BitSet(smoothBound);
        smoothNumber.set(1);
        smoothNumberReverse.set(smoothBound - 1);
        int smoothCount = 0;
        for (long j = 2; j < smoothBound; j++) {
            double baseSize = factorBaseSize(j * j);
            SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
//                boolean isSmooth = smallFactoriser.factor((int) j, baseSize, factors);
            // since we store the high factors for small numbers use the maximal factor base
            if (j % 500000 == 0)
                System.out.println("initializing : " + j + " " + ((j* 100.0) /smoothBound) + "% done");
            int maxFactorIndex = smallFactoriser.findMaxFactorIndex((int) j, (int) ceil(baseSize)) + 1;
            if (maxFactorIndex < baseSize) {
                smoothNumber.set((int) j);
                smoothNumberReverse.set((int) (smoothBound - j));
            }
        }
        System.out.println("initialization ready ");
    }
}
