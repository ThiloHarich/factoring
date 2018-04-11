package factoring.fermat.lehman;

import java.util.Collection;

import factoring.SingleLongFactorFinder;
import factoring.FactorizationOfLongs;

/**
 * This algorithm combines different fast Implementation of the lehman factoring algorithm.
 * It uses a fast implementaion of a trial division algorithm as well.
 * Usually it is 3 times faster then the YAFU implementation by... 
 * It works for numbers up to 41 bits.
 * For numbers up to around 34 bits it faster then the java implementation of the SQUAFU algorithm provided 
 * by til neuman. But you have to know about the size of the factors in the factorization. Where as in the SQUAFU
 * implementation usually the running time is independent of the size of the numbers.
 * For small factors (below n^1/3 if the number to factorize is n) the algorithm basically uses trial division.
 * Here it is much faster then the SQUAFU algorithm. If n consists of only two factors in the area of n^1/2 it 
 * uses the lehman algorithm first preferably high numbers > 3*n^1/3, then switching to trial division.
 * For numbers in between a mixed approach is going to be applied.
 * For random numbers it should always be faster then the SQUAFU up
 * 
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanFactorization implements FactorizationOfLongs {
	private final SingleLongFactorFinder impl;

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
	public LehmanFactorization(int bitsOfNumber, float factorExponent) {
//		this.bitsOfNumber = bitsOfNumber;
//		this.maxFactorMultiplier = maxFactorMultiplier;
		if (factorExponent < .33) {
			impl = new LehmanLongFactorFinder(bitsOfNumber, 1f);
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
				impl = new LehmanLongFactorFinder(bitsOfNumber, 3f);
		}
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		return getImpl(n).findFactors(n, primeFactors);
	}

	//	@Override
//	public BigInteger findSingleFactor(BigInteger n) {
//		//        LehmanLongFactorFinder impl = new LehmanLongFactorFinder(41, 1.001f);
//		final long factor = getImpl(n.longValue()).findFactors(n.longValue(), null);
//		return BigInteger.valueOf(factor);
//	}

//	/**
//	 * This is just for be compatible with SingleFactorFinder.
//	 * has the drawback to convert BigInteger to long and the conversion from
//	 * TreeMultiset to SortedMultiset_BottomUp
//	 * TODO what is the performance drawback
//	 * @param n
//	 * @return
//	 */
//	@Override
//	public SortedMultiset<BigInteger> factor(BigInteger n) {
//		TreeMultiset<Long> allFactors;
//		allFactors = factorization(n.longValue());
//		SortedMultiset result = new SortedMultiset_BottomUp(allFactors);
//		return result;
//	}


	public String toString() {
		return impl.toString();
	}

	@Override
	public boolean returnsOnlyPrimeFactors() {
		return factorizationByPrimes;
	}

	@Override
	public SingleLongFactorFinder getImpl(long n) {
		// TODO return the best possible impl here,
		return impl;
	}
}
