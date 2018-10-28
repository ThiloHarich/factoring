package factoring;

/**
 * TODO Seems to be unused
 * This tries to find only one prime factor of a long value n.
 * Is is not guaranteed to factor the number completely. If a {@link #setMaxFactor(int) maxFactor} is set,
 * it has to provide a factor by {@link #findPrimeFactor(long)} below maxFactor, if such a factor exists.
 * Example : SmallestPrimeFactorTable
 * @Deprecated("use the common interface{@link FactorFinder}")
 * Created by Thilo Harich on 27.03.2018.
 */
@Deprecated
public interface SinglePrimeFactorFinder extends SingleLongFactorFinder{


	/**
	 * Find a single prime factor of the given number.
	 * @param n
	 * @return a factor of the number n. If no factor can be found (i.e. the number is prime) it returns n.
	 */
	@Override
	long findPrimeFactor(long n);
}
