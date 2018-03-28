package factoring;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by Thilo Harich on 26.03.2018.
 */
public interface FactorizationOfLongs {


    /**
     * prints out the factorization in a nice ways starting with the lowest factors.
     * @param n
     * @return
     */
    default String printFactorization(long n) {
        final TreeMultiset<Long> factors = factorization(n);

        final List<String> s = new ArrayList<String>();
        for (final Multiset.Entry<Long> entry : factors.entrySet()) {
            final int exponent = entry.getCount();
            String part = "" + entry.getElement();
            part += exponent == 1 ? "" : "^" + exponent;
            s.add(part);
        }
        return String.join(" * ", s);
    }

    default long findSingleFactor(long n) {
        return getImpl(n).findFactors(n, null);
    }

    /**
     * This method returns a complete factorization of n.
     * It uses the implementation returned by {@link #getImpl(long)} and calls
     * {@link FactorFinderLong#findFactors(long, Collection)}. This will return a factor.
     * This factor does not have to be a prime factor, and has to be factorized again by
     * findFactors().
     * @param n
     * @return
     */
    default TreeMultiset<Long> factorization(long n) {
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
        long factor1 = getImpl(n).findFactors(n, primeFactors);
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
        TreeMultiset<Long> subFactors1 = factorization(factor1);
        TreeMultiset<Long> subFactors2 = factorization(factor2);
        factorsEven.addAll(subFactors1);
        factorsEven.addAll(subFactors2);
        factorsEven.addAll(primeFactors);
        return factorsEven;
    }

    /**
     * return an implementation which perform good for the number n.
     *
     * @param n
     * @return
     */
    FactorFinderLong getImpl(long n);

}
