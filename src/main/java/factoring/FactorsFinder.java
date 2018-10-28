package factoring;

import java.util.Collection;

/**
 * This tries to find (prime) factors of a long value n.
 * Is is not guaranteed to factor the number completely. If a {@link #setMaxFactor(int) maxFactor} is set,
 * it has to provide that a factor by  {@link #findFactors(long, Collection)} maxFactor}below maxFactor if such a factor exists.
 * This is an interface which is optimized for Integer values lower then 64 bits, which fit in a Long
 * value. 
 * Created by Thilo Harich on 27.03.2018.
 * Deprecated use the common interface {@link FactorFinder}
 */
@Deprecated
public interface FactorsFinder extends FactorFinder{

    /**
     * To be able to get the full factorization fast, there is a parameter primeFactors where
     * the implementation should store all the prime factors. If {@link #setMaxFactor(int)} is called,
     * only factors lower then the maxPrimeFactor will be added. 
     * @param n the number to be factorized.
     * @param primeFactors if this collection is given possible prime factors should be stored here.
     * @return n divided by the prime factors of n added to the primeFactors collection. If n is prime n will be returned. If 1 is returned the number is factorized completely.
     */
    long findFactors (long n, Collection<Long> primeFactors);
    
    default boolean findsPrimes() {
    		return true;
    }
 }
