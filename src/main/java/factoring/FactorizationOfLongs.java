package factoring;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * This is an Interface which gives an {@link #factorization(long)} of a long number.
 * It calls {@link SingleLongFactorFinder#findFactors(long, Collection)} at the Implementation returned by {@link #getImpl(long)}.
 * Created by Thilo Harich on 26.03.2018.
 */
public interface FactorizationOfLongs extends SingleLongFactorFinder {


    /**
     * If the {@link SingleLongFactorFinder#findFactors(long, Collection)} of the underlying {@link #getImpl(long)}
     * always retuns prime Factors this method returns true.
     * It then fills the prime factors of the given collection only with prime factors or always retuns a prime
     * factor if no prime factors are passed over.
     * @return
     */
    default boolean returnsOnlyPrimeFactors(){
        return true;
    }

    /**
     * returns a full prime factorization of the number.
     * The factorization is given as a TreeMultiset of the prime factors.
     * If the underlying algorithm returns only prime factors, we do not have to factorize the factors we have found.
     *
     * @param n
     * @return
     */
    default TreeMultiset<Long> factorization(long n) {
        TreeMultiset<Long> allFactors;
        if (returnsOnlyPrimeFactors())
            allFactors = factorizationByPrimes(n);
        else
            allFactors = factorizationByFactors(n);
        return allFactors;
    }


    /**
     * prints out the factorizationByFactors in a nice ways starting with the lowest factors.
     * @param n
     * @return
     */
    default String printFactorization(long n) {
        final TreeMultiset<Long> factors = factorizationByFactors(n);

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
     * This method returns a complete factorizationByFactors of n.
     * It uses the implementation returned by {@link #getImpl(long)} and calls
     * {@link SingleLongFactorFinder#findFactors(long, Collection)}. This will return a factor.
     * This factor does not have to be a prime factor, and has to be factorized again by
     * findFactors().
     * 
     * @deprecated ("due to performance reasons this should no be used. It can be called to check correctness of an algorithm or if there is no algorithm available
     * where {@link #returnsOnlyPrimeFactors} returns false.")
     * @see #factorizationByPrimes if possible.
     * @param n
     * @return
     */
    default TreeMultiset<Long> factorizationByFactors(long n) {
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
        TreeMultiset<Long> subFactors1 = factorizationByFactors(factor1);
        TreeMultiset<Long> subFactors2 = factorizationByFactors(factor2);
        factorsEven.addAll(subFactors1);
        factorsEven.addAll(subFactors2);
        factorsEven.addAll(primeFactors);
        return factorsEven;
    }

    default TreeMultiset<Long> factorizationByPrimes(long n) {
        TreeMultiset<Long> primeFactors = TreeMultiset.create();
        while ((n & 1) == 0)
        {
            primeFactors.add(2l);
            n = n >> 1;
        }
        if (n == 1) {
            return primeFactors;
        }
        // try to find all prime factors
        long remainder = getImpl(n).findFactors(n, primeFactors);
        // if we do not find a trivial divisor add it; this should only be the case if n
        // without the powers of 2 is a prime
        if (remainder != 1){
            primeFactors.add(n);
        }
        return primeFactors;
    }

    /**
     * return an implementation which perform good for the number n.
     *
     * @param n
     * @return
     */
    default SingleLongFactorFinder getImpl(long n){
        return this;
    }

}
