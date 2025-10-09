package de.harich.thilo.factoring.hart.multiplier;

import de.harich.thilo.factoring.hart.calculator.AdjustXModPow2Calculator;
import de.harich.thilo.factoring.hart.util.MultiplierFactorizationCount;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FindGoodMultiplierCalculator extends AdjustXModPow2Calculator {

    final int analyzeUpToIndex;
//    public int[][] multiplierFactorizationCount;
    public MultiplierFactorizationCount[] multiplierFactorizationCount;
    public double[][] multiplierFactorizationQuality;
    Set<Integer> mustHaveMultiplier = new HashSet<>();
    int numHitsPerNumber = 0;
    int lowestSuccessfulIndex = Integer.MIN_VALUE;
    List<Integer> successfulMultipliers = new ArrayList<>();


    public FindGoodMultiplierCalculator(int analyzeUpToIndex) {
        for (int i = 1; i < multiplier.length; i++) {
            sqrtMultiplier[i] = Math.sqrt(getMultiplier(i));
        }
        this.analyzeUpToIndex = analyzeUpToIndex;
//        multiplierFactorizationCount = new int[MODULUS_FOR_ADJUST_X][analyzeUpToIndex / MODULUS_FOR_ADJUST_X];
        multiplierFactorizationCount = new MultiplierFactorizationCount[analyzeUpToIndex];
        for (int i = 0; i < analyzeUpToIndex; i++) {
            multiplierFactorizationCount[i] = new MultiplierFactorizationCount(i);
        }
//        multiplierFactorizationQuality = new double[MODULUS_FOR_ADJUST_X][analyzeUpToIndex / MODULUS_FOR_ADJUST_X];
    }

    @Override
    public long handleSolution(long greatestCommonDivisor, int multiplierIndex, long numberToFactorize) {
//        multiplierFactorizationCount[multiplierIndex % MODULUS_FOR_ADJUST_X][multiplierIndex / MODULUS_FOR_ADJUST_X]++;
            multiplierFactorizationCount[multiplierIndex].count++;

        return greatestCommonDivisor;
//        if (multiplierIndex >= analyzeUpToIndex) {
//            if (numHitsPerNumber == 0){
//                System.out.println("no Multiplier foud in range");
//            }
//            if (numHitsPerNumber <= 3){
//                System.out.println(numberToFactorize + " hard number, only " + (numHitsPerNumber -1)  + " multipliers successful, first factorized by : " + lowestSuccessfulIndex);
//                mustHaveMultiplier.add(lowestSuccessfulIndex);
//
//            }
//            for (int i = 0; i < successfulMultipliers.size(); i++) {
//                Integer successfulMultiplier = successfulMultipliers.get(i);
//                double weight = (1.0 / (numHitsPerNumber - 0.9)) // less Factorisations -> more weight, only one -> weight 10
//                        * (1.0 / (i+1));
//                multiplierFactorizationQuality[successfulMultiplier % MODULUS_FOR_ADJUST_X][successfulMultiplier / MODULUS_FOR_ADJUST_X]
//                        += weight;                      // lower Index -> more weight
//            }
//
//            return greatestCommonDivisor;
//        }
//        if (numHitsPerNumber == 0){
//            lowestSuccessfulIndex = multiplierIndex;
//        }
//        numHitsPerNumber++;
//        successfulMultipliers.add(multiplierIndex);
//        multiplierFactorizationCount[multiplierIndex % MODULUS_FOR_ADJUST_X][multiplierIndex / MODULUS_FOR_ADJUST_X]++;
//
//
//        return 1;
    }
    @Override
    public boolean hasNextMultiplier(int index) {
//        return true;
        return index <= analyzeUpToIndex;
    }

    @Override
    public long getMultiplier(int index) {
        return index * 315L;
    }

    @Override
    public long findSmallFactor(long numberToFactorize) {
        numHitsPerNumber = 0;
        lowestSuccessfulIndex = Integer.MIN_VALUE;
        successfulMultipliers = new ArrayList<>();
        return 1;
    }
}
