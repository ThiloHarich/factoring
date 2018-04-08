package factoring;

import java.util.Collection;

/**
 * This tries to find only one factor of a long value n.
 * It is a simplified version of de.tilman_neumann.math.factor.SingleFactorFinder
 * Created by Thilo Harich on 27.03.2018.
 */
public interface SingleLongFactorFinder {

    /**
     * This is an interface which is optimized for Integer values lower then 64 bits, which fit in a Long
     * value. To be able to get the full factorizationByFactors fast, there is a parameter primeFactors where
     * the implementation can store factors which are prime.
     * @param n the number to be factorized.
     * @param primeFactors if this collection is given possible prime factors should be stored here.
     * @return a factor of n if there is one. If n is prime n will be returned.
     */
    long findFactors (long n, Collection<Long> primeFactors);

    /**
     * Find a single factor of the given N, which is composite and odd.
     * It is a simplified version of de.tilman_neumann.math.factor.SingleFactorFinder
     * It is like calling findFactors with no prime Factor list
     * @param n
     * @return factor
     */

    default long findSingleFactor(long n) {
        return findFactors(n, null);
    }
}
