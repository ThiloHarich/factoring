package factoring;

import java.util.Collection;

/**
 * Created by Thilo Harich on 14.01.2018.
 */
public abstract class FindPrimeFact implements Factorizer {

	/**
	 * We take the number n call {@link #findPrimeFactors(long, Collection)}.
	 * The call to this methods adds all prime factors to the factors.
	 * The factor returned by this method has to be a prime as well.
	 * @param n
	 * @param factors
	 * @return
	 */
	@Override
	public Collection<Long> storeFactors(long n, Collection<Long> factors) {
		// first make the number odd
		while ((n & 1) == 0)
		{
			factors.add(2l);
			n = n >> 1;
		}
		if (n>1) {
			final long remainder = findPrimeFactors(n, factors);
			if (remainder != 1)
				factors.add(remainder);
		}
		return factors;
	}

	/**
	 * each call either finds one ore more prime factors of n and adds it to the factors.
	 * There might already be some primes stored in the collection.
	 * The returned number if n divided by the factors. So by checking if the returned number
	 * is n indicates that the call has not found any factors.
	 * In order to have a prime factor decomposition of the number you have to ensure that the
	 * @param n
	 * @param factors
	 * @return a divisor of n.
	 */
	public abstract long findPrimeFactors(long n, Collection<Long> factors);


}
