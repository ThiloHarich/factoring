package de.harich.thilo.factoring.hart.multiplier;

import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.function.Predicate;

/**
 * a class, which can generate multipliers, which can be uses as input for the HartFactorisationAlgorithm.
 * The generated Multipliers should result in fast factorizations.
 * The generation splits in two parts, the generated numbers itself, and a filter, which eliminates numbers, which
 * do not lead to (many) factorizations. We also tried sorting the multipliers, but this was not helpful.
 * This creation process does not has to be fast, since this is done only once.
 *
 */
public class MultiplierFactory {

    public static long[] createMultiplesOf315(int limit, Predicate<Integer> filter){
        long baseMultiplier315 = 3 * 3 * 5 * 7 * 11;
//        long baseMultiplier315 = 3 * 3 * 5 * 7;
//        long baseMultiplier315 = 3 * 3 * 5;
        long[] multiplier = new long[limit];
        int index = 1;
        for (int i=1; i<limit; i++) {
//            if (filter.test(i))
            multiplier[index++] = i* baseMultiplier315;
        }
        multiplier[0] = index;
        return multiplier;
    }

    /**
     * Removing some numbers helps speeding up the Hart Algorithm.
     * It is not obvious which numbers to not use by the Algorithm.
     * One assumption is that smooth numbers, which have many ways to be split in two sets (dominator und numerator)
     * should work well, since a number to be factorized can be approximated by such a fraction.
     * A close approximation leads to a factorization.
     *
     * @param i
     * @return
     */
    public static boolean filterBigPrimesAndMultipleOf9(Integer i){

        // multiples of 9 do not work out good. Found by looking at how often multipliers i (with i % 9 == 0)
        // factorize random semiprimes
        boolean isAMultipleOf9 = i % 9 == 0;
        if (isAMultipleOf9)
            return false;

        // a prime by definition only consists out of one factor and such does not increase the number of
        // possible fractions combined with the base multiplier.

        boolean isABigPrime = i > 100 && BigInteger.valueOf(i).isProbablePrime(10);
        if (isABigPrime)
            return false;

        // the following filters do not help

        // small numbers i % 8 == 0 do not lead to many factorisations. But when removing such multipliers,
        // the algorithm does not speed up. such numbers might be useful for bigger multipliers
//        boolean isAMultipleOf8 = i % 8 == 0 ;
//        if (isAMultipleOf8)
//            return false;

        // squares only have two multipliers, but are seldom.
        int sqrt = (int) Math.sqrt(i);
        boolean isABigSquare = i > 100 && sqrt * sqrt == i && i % 2 == 0;

//        if (isABigSquare)
//            return false;

        return true;
    }

    /**
     * Numbers of the form 3*3*5 * i, 3*3*5*7 * i and 3*3*5*7*11 * i, a little advantage over numbers of the form 3*3*5 * i
     */
    public static long[] createMultiplesOfPrimesBelow11(int limit, Predicate<Integer> filter){
        final long SMALL_MULTIPLIER = 3*3*5;
        final long MID_MULTIPLIER =   3*3*5*7;
        final long BIG_MULTIPLIER =   3*3*5*7*11;

        long[] goodMultipliers = new long[limit];

        int index = 1;
        for (int i=1; i<limit / 3; i++) {
            if (filter.test(i)) {
                goodMultipliers[index++] = i * BIG_MULTIPLIER;

                boolean noCollision = i % 11 != 0;
                if (noCollision) {
                    goodMultipliers[index++] = i * MID_MULTIPLIER;

                    noCollision = i % 7 != 0;
                    if (noCollision) {
                        goodMultipliers[index++] = i * SMALL_MULTIPLIER;
                    }
                }
            }
        }
        // store the length in the first cell; depending on the filter we apply the size might shrink
        goodMultipliers[0] = index;
        return goodMultipliers;
    }

    public static long [] fromFile(String pathToMultiplierFile) {

        long[] multiplier = new long[0];
        if (Files.exists(Path.of(pathToMultiplierFile))) {
            List<String> lines;
            int i = 0;
            try {
                lines = Files.readAllLines(Path.of(pathToMultiplierFile));
                final String line = lines.get(0);
                final String[] numbers = line.split(",");
                multiplier = new long[Integer.parseInt(numbers[i++])];
                for (; i < numbers.length; i++) {
                    final String number = numbers[i];
                    if (!number.isEmpty()) {
                        final int k = Integer.parseInt(number);
                        double quality = 0;
                        if (BigInteger.valueOf(k).isProbablePrime(5)){
                            quality = (i + 0.0) / k;
                        }
                        multiplier[i-1] = k;
//                        sqrtMultiplier[i-1] = Math.sqrt(k);
                    }
                }
            } catch (final IOException e) {
                e.printStackTrace();
            }
        }
        return multiplier;
    }
}
