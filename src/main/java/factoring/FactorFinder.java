package factoring;

import java.util.Collection;

/**
 * This is a generic interface for classes which can find some factors of a number up to a certain factor.
 * This maximal factor can be given by {@link #setMaxFactor(int)}.
 * The implementation can either return prime or composite factors of the number.
 * It can also choose if it only gives back one factor per call to {@link #findFactors(long, Collection)} as a return value,
 * or stores every factor it finds in the collection.
 *
 * @author thiloharich
 *
 */
public interface FactorFinder {

	/**
	 * Sets the maximalfactor the FactorFinder should look for factors.
	 * @param maxTrialFactor
	 */
	default void setMaxFactor(int maxPrimeFactor)
	{
	}

	/**
	 * Gives back at least one factor not exceeding {@link #setMaxFactor(int) maxFactor} of the number n, if there is any.
	 * If {@link #findsOneFactor()} is true, the factor will be returned by the return of the function. In this case the factors
	 * parameter will not be used. If {@link #findsOneFactor()} is false the factors will be stored in the factors collection
	 * then the return is n divided by the factors of n added to the factors collection.
	 * If {@link #findsPrimes()} the implementation should return only prime factors (either as return value or in the factors collection). If {@link #setMaxFactor(int)} is called,
	 * only factors lower then the maxPrimeFactor will be added.
	 *
	 * @param n the number to be factorized.
	 * @param factors Only If {@link #findsOneFactor()} is false the factors will be stored in this collection. In this case factors can not be null.
	 * @return If {@link #findsOneFactor()} is true it will return a factor of the number n. If {@link #findsPrimes()}  it is a prime factor.
	 * If {@link #findsOneFactor()} is false, the return is n divided by the factors of n added to the factors collection.
	 * If n is prime n will be returned. If 1 is returned the number is factorized completely.
	 */
	long findFactors (long n, Collection<Long> factors);

	/**
	 * Indicates that the implementation will only return prime factors of the number in {@link #findFactors(long, Collection)}.
	 * @return
	 */
	default boolean findsPrimes() {
		return true;
	}

	/**
	 * If the implementation will only return one factor as return value of {@link #findFactors(long, Collection)} and not store it in the factors collection, this
	 * function returns true.
	 * If this returns false, in {@link #findFactors(long, Collection)} the factors will be stored in the second parameter. The return of the function will just be the
	 * remainder of n divided by the factors stored in the collection.
	 * @return
	 */
	default boolean findsOneFactor() {
		return false;
	}



}
