package de.harich.thilo.factoring.hart.multiplier;

import de.harich.thilo.factoring.hart.HartFactorizationAlgorithm;
import de.harich.thilo.factoring.hart.util.MultiplierFactorizationCount;
import de.harich.thilo.factoring.hart.util.MultiplierFactorizationQuality;
import factoring.primes.Primes;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static de.harich.thilo.factoring.hart.multiplier.FindGoodMultiplierCalculator.MODULUS_FOR_ADJUST_X;

public class GoodMultiplierMod8Finder {

    private static final String MULTIPLIER_FILE = "./multipliers.txt";
    private final boolean readFromFile = false;
    private FindGoodMultiplierCalculator findGoodMultiplierCalculator;

    public int[] findMultipliers (int maxMultiplier, int bitsOfNumbersToFactorize, int numberSemiprimes){
        int runtimeOfHart = 1 << (bitsOfNumbersToFactorize / 3);
        long[] semiprimes = Primes.makeNumberList(bitsOfNumbersToFactorize, runtimeOfHart, numberSemiprimes, readFromFile);
//        long[] semiprimes = Primes.makeSemiPrimesList(bitsOfNumbersToFactorize, numberSemiprimes, readFromFile);

        findGoodMultiplierCalculator = new FindGoodMultiplierCalculator(maxMultiplier);
        HartFactorizationAlgorithm factorizer = new HartFactorizationAlgorithm(findGoodMultiplierCalculator);

        for (final Long semiprime : semiprimes) {
             factorizer.findSingleFactor(semiprime);
        }
        for (int i = 5000; i < 5100; i++) {
            System.out.println(findGoodMultiplierCalculator.multiplierFactorizationCount[i]);
        }

        Arrays.sort(findGoodMultiplierCalculator.multiplierFactorizationCount);
        int index = 0;
        int goodIndex = 0;
        int[] goodMultipliers = new int[findGoodMultiplierCalculator.multiplierFactorizationCount.length];
        for (MultiplierFactorizationCount multiplierFactorizationCount : findGoodMultiplierCalculator.multiplierFactorizationCount){
//            if (multiplierFactorizationCount.multiplier > 100 && Primes.isProbablePrime(multiplierFactorizationCount.multiplier)){
            if (multiplierFactorizationCount.count < 100){
                index = index;
            }
            if (index < multiplierFactorizationCount.multiplier * 5){
                goodMultipliers[goodIndex++] = multiplierFactorizationCount.multiplier;
            }
            index++;
        }
//        MultiplierFactorizationQuality[][] goodMultiplierQuality = new MultiplierFactorizationQuality[8][maxMultiplier/ MODULUS_FOR_ADJUST_X];
//        // find the pos, where all remainders have a Multiplier with count > 0
//        int goodMultipliersForAll = (maxMultiplier/ MODULUS_FOR_ADJUST_X)-1;
//        int goodMultiplierForOne = 0;
//        for (int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++){
//            goodMultiplierQuality[remainder] = findGoodMultiplierQuality(remainder);
//            int goodMultiplerIndex = (maxMultiplier/ MODULUS_FOR_ADJUST_X)-1;
//            while (goodMultiplierQuality[remainder][goodMultipliersForAll].quality == 0){
//                goodMultipliersForAll--;
//            }
//            while (goodMultiplierQuality[remainder][goodMultiplerIndex].quality == 0){
//                goodMultiplerIndex--;
//            }
//            goodMultiplierForOne = Math.max(goodMultiplierForOne, goodMultiplerIndex);
//        }
//        MultiplierFactorizationCount[][] goodMultiplierCounts = new MultiplierFactorizationCount[8][maxMultiplier/ MODULUS_FOR_ADJUST_X];
//        // find the pos, where all remainders have a Multiplier with count > 0
//        goodMultipliersForAll = (maxMultiplier/ MODULUS_FOR_ADJUST_X)-1;
//        goodMultiplierForOne = 0;
//        for (int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++){
//            goodMultiplierCounts[remainder] = findGoodMultiplierCounts(remainder);
//            int goodMultiplerIndex = (maxMultiplier/ MODULUS_FOR_ADJUST_X)-1;
//            while (goodMultiplierCounts[remainder][goodMultipliersForAll].count == 0){
//                goodMultipliersForAll--;
//            }
//            while (goodMultiplierCounts[remainder][goodMultiplerIndex].count == 0){
//                goodMultiplerIndex--;
//            }
//            goodMultiplierForOne = Math.max(goodMultiplierForOne, goodMultiplerIndex);
//        }
//        // merge the lists, such that each MODULUS_FOR_ADJUST_X we have the same remainder, and have a good Multiplier
//        int indexInResult = 0;
//        int[] goodMultipliers = new int[MODULUS_FOR_ADJUST_X * (goodMultipliersForAll + 1)];
//        for (int indexLevel = 0; indexLevel <= goodMultipliersForAll; indexLevel++){
//            for (int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++){
//                goodMultipliers[indexInResult++] = goodMultiplierCounts[remainder][indexLevel].multiplier;
//            }
//        }
        storeMultipliers(goodMultipliers);
        return goodMultipliers;
    }

    public void storeMultipliers(int[] multipliers){
        String delimiter = ",";
        String factorSequence = multipliers.length + delimiter;
        for (final int multiplier : multipliers) {
            factorSequence += multiplier + delimiter;
        }
        final Path path = Paths.get(MULTIPLIER_FILE);
        final byte[] strToBytes = factorSequence.getBytes();

        try {
            Files.write(path, strToBytes);
        } catch (IOException e) {
            //
        }
        System.out.println("wrote " + multipliers.length + " factors ");
    }

//    private MultiplierFactorizationCount[] findGoodMultiplierCounts(int remainder) {
//        final List<MultiplierFactorizationCount> multiplierFactorizationCounts = new ArrayList<>();
//        for (int i = 0; i < findGoodMultiplierCalculator.analyzeUpToIndex / MODULUS_FOR_ADJUST_X; i++) {
//            int multiplier = MODULUS_FOR_ADJUST_X * i + remainder;
//            final int count = findGoodMultiplierCalculator.multiplierFactorizationCount[remainder][i];
//
//            final MultiplierFactorizationCount multiplierFactorizationCount = new MultiplierFactorizationCount(multiplier, count);
//            multiplierFactorizationCounts.add(multiplierFactorizationCount);
////            if (i < 10)
////                System.out.println(multiplier + "\t" + count + "   " + multiplier + " = "
////                            + factorization.getFactorisationAsString((long) multiplierFactorizationCount.multiplier)
////                );
//        }
//        final MultiplierFactorizationCount[] factorizationCounts = new MultiplierFactorizationCount[multiplierFactorizationCounts.size()];
//        multiplierFactorizationCounts.toArray(factorizationCounts);
//        Arrays.sort(factorizationCounts);
//        return factorizationCounts;
//    }

    private MultiplierFactorizationQuality[] findGoodMultiplierQuality(int remainder) {
        final List<MultiplierFactorizationQuality> multiplierFactorizationQuality = new ArrayList<>();
        for (int i = 0; i < findGoodMultiplierCalculator.analyzeUpToIndex / MODULUS_FOR_ADJUST_X; i++) {
            int multiplier = MODULUS_FOR_ADJUST_X * i + remainder;
            final double quality = findGoodMultiplierCalculator.multiplierFactorizationQuality[remainder][i];

            final MultiplierFactorizationQuality multiplierFactorizationCount = new MultiplierFactorizationQuality(multiplier, quality);
            multiplierFactorizationQuality.add(multiplierFactorizationCount);
//            if (i < 10)
//                System.out.println(multiplier + "\t" + count + "   " + multiplier + " = "
//                            + factorization.getFactorisationAsString((long) multiplierFactorizationCount.multiplier)
//                );
        }
        final MultiplierFactorizationQuality[] factorizationQualities = new MultiplierFactorizationQuality[multiplierFactorizationQuality.size()];
        multiplierFactorizationQuality.toArray(factorizationQualities);
        Arrays.sort(factorizationQualities);
        return factorizationQualities;
    }

    public static void main(String[] args) {
        GoodMultiplierMod8Finder finder = new GoodMultiplierMod8Finder();

//        finder.findMultipliers(maxMultiplier, 12, 20);
        int bitsOfNumbersToFactorize = 42;
        int maxMultiplier = (int) Math.cbrt(1L << bitsOfNumbersToFactorize);
        maxMultiplier = 200000;
        int numSemiprimes = maxMultiplier * 1;
        finder.findMultipliers(maxMultiplier, bitsOfNumbersToFactorize, numSemiprimes);
    }
}
