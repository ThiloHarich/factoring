package factoring.fermat.lehman;

import com.google.common.collect.TreeMultiset;

/**
 * Created by Thilo Harich on 26.03.2018.
 */
public class AbstractFactorFinder {

    // TODO define an intrface
    LehmanReverseFact impl;
//LehmanNoSqrtFact impl;

    public TreeMultiset<Long> factorLong(long n) {
        // if we have a prime return an empty set
        TreeMultiset<Long> factorsEven = TreeMultiset.create();
        while ((n & 1) == 0)
        {
            factorsEven.add(2l);
            n = n >> 1;
        }
        if (n == 1) {
            return factorsEven;
        }
        TreeMultiset<Long> primeFactors = TreeMultiset.create();
        // find one factor and decomposite this factor and n/factor
        long factor1 = impl.findPrimeFactors(n, primeFactors);
        // if we do not find a divisor just return it
        if (factor1 == n){
            factorsEven.add(n);
            return factorsEven;
        }
        // also divide out the prime factorsEven
        long factor2 = n/factor1;
        for (long factor : primeFactors) {
            factor2 /= factor;
        }
        TreeMultiset<Long> subFactors1 = factorLong(factor1);
        TreeMultiset<Long> subFactors2 = factorLong(factor2);
        factorsEven.addAll(subFactors1);
        factorsEven.addAll(subFactors2);
        factorsEven.addAll(primeFactors);
        return factorsEven;
    }
}
