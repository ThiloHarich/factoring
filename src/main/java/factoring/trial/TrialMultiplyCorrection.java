package factoring.trial;

import java.util.Collection;

import factoring.FactorizationOfLongs;
import factoring.primes.Primes;

/**
 * This implementation is generating a list of all primes up to a limit.
 * Beside of storing the prime itself, it also store the reciprocal value.
 * When checking a number if it is divisible by the prime, we will not divide
 * the number by the prime, we will multiply by the inverse, since this is faster.
 * Due to precision we have to check for a given range near an Integer.
 * And then do a Long division.
 * This implementation is around two and a half times faster then a version  based on long numbers.
 * Since Double only has 52 bis for the remainder, this can only work for numbers below 2^52.
 * We can only factorize numbers up to maxFactor^2
 * When calling it with bigger numbers only prime factors below
 * maxFactor were added to the factors. {@link #findFactors(long, Collection)} then might return a
 * composite number.
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialMultiplyCorrection implements FactorizationOfLongs {

	// Size of number is ~ 2^52
	private static final int DISCRIMINATOR_BITS = 10; // experimental result
	private static final double DISCRIMINATOR = 1.0/(1<<DISCRIMINATOR_BITS);
	private int maxFactor = 65535;
	double[] primesInv;
	int[] primes;

	public TrialMultiplyCorrection(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		//		initPrimes();
		final Primes primesGen = Primes.initPrimesEratosthenes(maxFactor);
		primes = primesGen.primes;
		primesInv = primesGen.primesInv;
	}


	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		final int Nbits = 64-Long.numberOfLeadingZeros(n);
		final int multiplicationWorksBits = Nbits - 53 + DISCRIMINATOR_BITS;
		int primeIndex = 0;
		if (multiplicationWorksBits>0) {
			// for the smallest primes we must do standard trial division
			final int multiplicationWorks = 1<<multiplicationWorksBits;
			for (; primes[primeIndex] < multiplicationWorks; primeIndex++) {
				if (n%primes[primeIndex]==0) {
					return primes[primeIndex];
				}
			}
		}
		final int primeLimit = (int) Math.min(Math.sqrt(n), maxFactor);
		for (; primes[primeIndex] <= primeLimit; primeIndex++) {
			// round the number. Casting to long is faster then rounding the double number itself, but we
			// have to prevent some cases were the number is not correctly rounded by adding a small number
			long nDivPrime = (long) (n*primesInv[primeIndex] + DISCRIMINATOR);
			// TODO if we want to return all factors this must we a while loop and some special handling
			if (nDivPrime * primes[primeIndex] == n) {
				if (primeFactors == null)
					return primes[primeIndex];
				primeFactors.add((long) primes[primeIndex]);
				n = nDivPrime;
				// if the remaining number is lower then the maximal factor we look for it must be a prime
				nDivPrime = (long) (n*primesInv[primeIndex] + DISCRIMINATOR);
			}
			// if the remaining part is less then the maximal prime, it must be a prime power
		}
		return -1;
	}


	@Override
	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}
}
