package de.harich.thilo.factoring;

import de.harich.thilo.factoring.hart.HartFactorizationAlgorithm;
import de.harich.thilo.factoring.hart.calculator.AdjustXModPow2Calculator;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Factorization {

    HartFactorizationAlgorithm algorithm = new HartFactorizationAlgorithm(new AdjustXModPow2Calculator());

    public String getFactorisationAsString(Long number){
        List<Long> sortedFactors = getSortedPrimeFactors(number);
        String factorisation = "";

        long lastFactor = sortedFactors.get(0);
        int exponent = 1;
        
        for (int factorIndex = 1; factorIndex < sortedFactors.size(); ){
            long factor = sortedFactors.get(factorIndex);
            if (factor == lastFactor){
                exponent++;
            }
            else{
                factorisation += lastFactor + ((exponent == 1) ? "" : ("^" + exponent + " "));
                exponent = 1;
            }
            lastFactor = factor;
            factorIndex++;
        }
        factorisation += lastFactor + ((exponent == 1) ? "" : ("^" + exponent + " "));

        return factorisation.trim();
    }

    public List<Long> getSortedPrimeFactors(Long number) {
        List<Long> factors = getPrimeFactors(number);
        factors.sort(Comparator.comparingLong(Long::longValue));
        return factors;
    }

    private List<Long> getPrimeFactors(Long number) {
        if (number < 4)
            return new ArrayList<>(List.of(number));
        List<Long> primeFactors = new ArrayList<>();
        do{
            long factor = algorithm.findSingleFactor(number);
            if (BigInteger.valueOf(factor).isProbablePrime(10))
                primeFactors.add(factor);
            else{
                primeFactors.addAll(getPrimeFactors(factor));
            }
            number /= factor;
        }
        while (number > 1);
        return primeFactors;
    }
}
