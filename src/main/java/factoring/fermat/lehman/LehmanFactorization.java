package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;

/**
 * This algorithm combines different fast implementations of the lehman factoring algorithm.
 * It uses a fast implementaion of a trial division algorithm as well.
 * Usually it is 3 times faster then the YAFU implementation by Warren D. Smith
 * It works for numbers up to 41 bits.
 * For those numbers it is usually faster then the java implementation of the SQUAFU algorithm provided by Til Neuman.
 * For random number it is always faster.
 * The performance is strongly related to the size of the factors.
 * Since we have to combine trial division for numbers below a*n^1/3 and the lehman phase for higher numbers,
 * this algorithm also needs information on the size of the (second highest) factor(s) to run fast.
 * In other algorithms like SQUAFU usually the running time is independent of the size of the numbers.
 * For small factors (below n^1/3 if the number to factorize is n) this algorithm basically uses trial division.
 * Here it is much faster then the SQUAFU algorithm. If n consists of only two factors in the area of n^1/2 it
 * uses the lehman algorithm first. It is also faster then SQUAFU.
 * For numbers in between we choose the breakpoint such that we always use trial divison.
 * Here is the worst case for this algorithm, here SQUAFU is usually faster for bigger numbers.
 *
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanFactorization implements FactorizationOfLongs, FactorFinder {
	private final FactorFinder impl;

	boolean factorizationByPrimes = false;
	//	int bits;
	//	float maxFactorMultiplier;

	// 3 cases
	// all factors + factors below n^1/3 -> use trial divison with small primes first
	// factors around n^1/3 -> use trial divison with high primes first -> worst case
	// factors above n^1/3 -> use fermat with high primes first

	/**
	 * Configure the algorithms to be used, depending on the size of the number to be factorized and the
	 * size of the factors in the factorization. If the Input is n and the the second biggest number in the
	 * factorization is n^a, use a as the factorExponent.
	 *
	 * If the second biggest factor is n^1/3 (a = 1/3) or lower, then there can at least be one prime factor greater then n^1/3.
	 * In this case the algorithm will processed numbers up to n^1/3 with trial division by {@link LehmanFactorFinder} first.
	 * Then the factors above n^1/3 were analyzed with a simple lehman process. The found factor must be a prime factor in this case.
	 *
	 * For a > .0372 the factors above 3.2 * n^1/3 will be inspected first with the lehman phase in
	 * {@link LehmanFactorFinderMod12}. This uses also a mod 3 argument and is the fastest in this package.
	 * If the number is not fully factorized trial division is applied for primes below 3.2 * n^1/3.
	 *
	 *
	 * If 0.33 <= a < 0.37, use the factors below n^a will be factored out by trial division
	 * in {@link LehmanFactorFinderMod12}.
	 * Then we inspect the higher factors, in case there are any.
	 * This case is the worst case behavior, since the prime factors can not be found by one of the algorithms trial division and lehman
	 * exclusive and fast. It tries to find all factors below n^a with trial division first.
	 *
	 * If the numbers
	 * are random, or
	 * the second highest factors is lower then n^1/3, or
	 * you want to factorize a sequence of numbers
	 * then you should use
	 * 0 as factorExponent (or any number below .33).
	 *
	 * If you know the numbers n to be factorized are semiprimes -
	 * both with roughly the same size - the number have size n^1/2. Use .5 as factorExponent here.
	 *
	 * Only if you know that the factorization
	 * has two primes, where the lowest prime n^a has size around between n^1/3 (~ n^0.33) and n^0.37, use this a here.
	 *
	 *
	 *
	 * {@link LehmanReverseFact} is applied then and does trial division on numbers
	 * @param bitsOfNumber the maximal number of bits of the number to be factorized
	 * @param factorExponent the exponent of the second biggest factor of the factorizationByFactors. If you factorize n
	 *                          this is factor should have size n ^ factorExponent.
	 */
	public LehmanFactorization(int bitsOfNumber, float factorExponent) {
		//		this.bitsOfNumber = bitsOfNumber;
		//		this.maxFactorMultiplier = maxFactorMultiplier;
		if (factorExponent < .33) {
			impl = new LehmanFactorFinder(bitsOfNumber, 1f, true);
			factorizationByPrimes = true;
		}
		else {
			// TODO find the right boarder here
			// 2^(40/3)*3,2 = 2^(40 *(1/3+ log_2(3.2)/40) = 2^(40 *(1/3+ .042)
			final float mult = (float) ((factorExponent - .28) * bitsOfNumber);
			if (factorExponent <= .38)
				// there is still place for improvement here
				impl = new LehmanFactorFinderMod12(bitsOfNumber, mult, true);
			else
				// we use the same algorithm as for small numbers but configure it to first look for big numbers.
				// in this case we will get no primes out of the algorithm, but since the numbers are big we do not care
				impl = new LehmanFactorFinderMod12(bitsOfNumber, 3.2f, false);
		}
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		return getImpl(n).findFactors(n, primeFactors);
	}

	//	@Override
	//	public BigInteger findSingleFactor(BigInteger n) {
	//		//        LehmanFactorFinder impl = new LehmanFactorFinder(41, 1.001f);
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


	@Override
	public String toString() {
		return impl.toString();
	}

	@Override
	public boolean findsPrimes(){
		return factorizationByPrimes;
	}

}
