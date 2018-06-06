package factoring.fermat.lehman;

import java.util.Collection;

import factoring.SingleLongFactorFinder;
import factoring.FactorizationOfLongs;

/**
 * This algorithm combines different fast Implementation of the lehman factoring algorithm.
 * It uses a fast implementaion of a trial division algorithm as well.
 * Usually it is 3 times faster then the YAFU implementation by Warren D. Smith
 * It works for numbers up to 41 bits.
 * For those numbers it is usually faster then the java implementation of the SQUAFU algorithm provided
 * For random number it is always faster.
 * The performance is strongly related to the size of the factors.
 * Since we have to combine trial division for numbers below a*n^1/3 and the lehman phase for higher numbers,
 * this algorithm also needs information on the size of the (second highest) factor(s) to run fast.
 * In other algorithms like SQUAFU
 *  usually the running time is independent of the size of the numbers.
 * For small factors (below n^1/3 if the number to factorize is n) this algorithm basically uses trial division.
 * Here it is much faster then the SQUAFU algorithm. If n consists of only two factors in the area of n^1/2 it 
 * uses the lehman algorithm first. It is also faster then SQUAFU.
 * For numbers in between we choose the breakpoint such that we always use trial divison.
 * Here is the worst case for this algorithm, here SQUAFU is usually faster for bigger numbers.
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
	 * size of the factors in the factorization. In the Input is n and the lowest number in the
	 * factorization is n^a, use a as the factorExponent.
     * For a < .3333 all prime factors up to n^1/3
     * were processed with trial division by {@link LehmanFactorFinder}
     * Here the lehman phase is used rather seldom and we use the simple version without (increased)
     * multiplier and mod 3 argument.
     *
     * For a > .0372 the factors above 3.2 * n^1/3 will be inspected first with the lehman phase in
     * {@link LehmanFactorFinderMod12}. Then trial division is done for primes below 3.2 * n^1/3.
     * This uses also a mod 3 argument.
     *
     * If 0.33 <= a < 0.37, use the factors below n^a will be factored out by trial division
     * in {@link LehmanFactorFinderMod12}.
     * Then we inspect the higher factors, in case there are any.
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
     * has two primes, where the lowest prim n^a has size around between n^1/3 (~ n^0.33) and n^0.37, use this a here.
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
            float mult = (float) ((factorExponent - .28) * bitsOfNumber);
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
