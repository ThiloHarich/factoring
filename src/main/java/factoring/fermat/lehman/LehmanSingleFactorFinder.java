package factoring.fermat.lehman;

import java.math.BigInteger;

import com.google.common.collect.TreeMultiset;
import de.tilman_neumann.math.factor.SingleFactorFinder;
import de.tilman_neumann.types.SortedMultiset;
import de.tilman_neumann.types.SortedMultiset_BottomUp;
import factoring.FactorFinderLong;
import factoring.FactorizationOfLongs;

/**
 *
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanSingleFactorFinder implements SingleFactorFinder, FactorizationOfLongs {
	private final FactorFinderLong impl;

	boolean factorizationByPrimes = false;
//	int bits;
//	float maxFactorMultiplier;

	// 3 cases
	// all factors + factors below n^1/3 -> use trial divison with small primes first
	// factors around n^1/3 -> use trial divison with high primes first -> worst case
	// factors above n^1/3 -> use fermat with high primes first

	/**
	 * Configure the algorithms to be used, depending on the size of the number to be factorized and the
	 * size of the factors in the factorization. In the Input is n and the second biggest number in the
	 * factorization is n^a, use a as the factorExponent.
	 * If the numbers are random, or are ascending
	 * the second biggest number of the factorization is small compared to n. In this case you should use
	 * 0 as factorExponent (or any number below .33). If you know the numbers n to be factorized are semiprimes -
	 * both with roughly the same size - the number have size n^1/2. This is
	 * the exponent of both primes are ~ .5. Use .5 as factorExponent here. Only if you know that the factorizationByFactors
	 * has two or three primes, where one prime has size around n^1/3, use 1/3 as factorExponent here.<br>
	 * For factorExponent < 0.33f we do trial division beginning with small primes is used first, then we use the lehman algorithm.<br>
	 * For factorExponent > 0.37f the lehman algorithm optimized for big factors is applied first.<br>
	 * For  0.33 <= factorExponent <= .37 trial division with factors in this range is done first.<br>
	 * @param bitsOfNumber the maximal number of bits of the number to be factorized
	 * @param factorExponent the exponent of the second biggest factor of the factorizationByFactors. If you factorize n
	 *                          this is factor should have size n ^ factorExponent.
	 */
	public LehmanSingleFactorFinder(int bitsOfNumber, float factorExponent) {
//		this.bitsOfNumber = bitsOfNumber;
//		this.maxFactorMultiplier = maxFactorMultiplier;
		if (factorExponent < .33) {
			impl = new LehmanFactorFinder(bitsOfNumber, 1f);
			factorizationByPrimes = true;
		}
		else {
			// TODO find the right boarder here
			if (factorExponent <= .37)
				// there is still place for improvement here
				impl = new LehmanReverseFact(bitsOfNumber, 3f);
			else
				// we use the same algorithm as for small numbers but configure it to first look for big numbers.
				// in this case we will get no primes out of the algorithm, but since the numbers are big we do not care
				impl = new LehmanFactorFinder(bitsOfNumber, 3f);
		}
	}

	@Override
	public BigInteger findSingleFactor(BigInteger n) {
		//        LehmanFactorFinder impl = new LehmanFactorFinder(41, 1.001f);
		final long factor = getImpl(n.longValue()).findFactors(n.longValue(), null);
		return BigInteger.valueOf(factor);
	}

	/**
	 * This is just for be compatible with SingleFactorFinder.
	 * has the drawback to convert BigInteger to long and the conversion from
	 * TreeMultiset to SortedMultiset_BottomUp
	 * TODO what is the performance drawback
	 * @param n
	 * @return
	 */
	@Override
	public SortedMultiset<BigInteger> factor(BigInteger n) {
		TreeMultiset<Long> allFactors;
		allFactors = factorization(n.longValue());
		SortedMultiset result = new SortedMultiset_BottomUp(allFactors);
		return result;
	}

	@Override
	public String getName() {
		return LehmanSingleFactorFinder.class.getSimpleName();
	}

	@Override
	public boolean returnsOnlyPrimeFactors() {
		return factorizationByPrimes;
	}

	@Override
	public FactorFinderLong getImpl(long n) {
		// TODO return the best possible impl here,
		return impl;
	}
}
