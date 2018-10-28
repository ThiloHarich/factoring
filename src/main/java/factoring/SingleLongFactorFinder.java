package factoring;

import java.util.Collection;

/**
 * This tries to find only one factor of a long number n.
 * It is a simplified version of de.tilman_neumann.math.factor.SingleFactorFinder
 * Is is not guaranteed to factor the number completely. If a {@link #setMaxFactor(int) maxFactor} is set,
 * it has to provide a factor by {@link #findFactors(long, Collection)} maxFactor} below maxFactor, if such a factor exists.
 * Example : Pollard Rho, which always just gives one not necessary prime factor
 * @Deprecated("use the common interface{@link FactorFinder}")
 * Created by Thilo Harich on 27.03.2018.
 */
@Deprecated
public interface SingleLongFactorFinder extends FactorFinder{


	/**
	 * Find a single factor of the given n, which is might be composite but odd.
	 * @param n
	 * @return a factor of the number n. If no factor can be found (i.e. the number is prime) it returns n.
	 */
	long findPrimeFactor(long n);
}
