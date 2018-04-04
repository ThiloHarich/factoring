package factoring;

import java.util.Collection;

/**
 * Created by Thilo Harich on 27.03.2018.
 */
public interface FactorFinderLong {

    /**
     * This is an interface which is optimized for Integer values lower then 64 bits, which fit in a Long
     * value. To be able to get the full factorizationByFactors fast, there is a parameter primeFactors where
     * the implementation can store factors which are prime.
     * @param n the number to be factorized.
     * @param primeFactors if this collection is given possible prime factors should be stored here.
     * @return a factor of n if there is one. If n is prime n will be returned.
     */
    long findFactors (long n, Collection<Long> primeFactors);

}
