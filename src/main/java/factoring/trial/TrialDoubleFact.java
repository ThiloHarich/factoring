package factoring.trial;

import java.util.Collection;

import factoring.FactorizationOfLongs;

/**
 * This implementation is based on double values instead of integer/long values.
 * Instead of dividing n the number to factorize by a possible factor,
 * we multiply n by the inverse of the possible factor.
 * This only provides performance benefits, if we will check multiple factors,
 * and store the reciprocal of all possible factors.
 * So we are generating a list of all primes up to a limit.
 * Beside of storing the prime itself, we also store the reciprocal value.
 * When checking a number if it is divisible by the prime, we will not divide
 * the number by the prime, we will multiply by the inverse, since this is faster.
 * Due to precision we have to check for a given range near an Integer.
 * And then do a Long division.
 * This implementation is around two times faster then a version  based on long numbers.
 * Since Double only has 52 bis for the remainder, this can only work for numbers below 2^52.
 * We can only factorize numbers up to maxFactor^2
 * When calling it with bigger numbers only prime factors below
 * maxFactor were added to the factors. {@link #findFactors(long, Collection)} then might return a
 * composite number.
 *
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialDoubleFact implements FactorizationOfLongs {

	// The number of values to be printed
	private static final int PRINT_NUM = 20000;
	// for printing we need a value a little bit above 1
	public static final double PRINT_CONST = 1.0000001;
	private int maxFactor = 65535;
	double[] primesInv;
	int[] primes;

	/**
	 * finds the prime factors up to maxFactor by the sieve of eratosthenes.
	 * Not optimized, since this is only called once when initializing.
	 */
	void initPrimesEratosthenes()
	{
		final double logMaxFactor = Math.log(maxFactor);
		final int maxPrimeIndex = (int) ((maxFactor) / (logMaxFactor - 1.1)) + 4;
		primesInv = new double [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		primes = new int [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		int primeIndex = 0;
		final boolean [] noPrimes = new boolean [maxFactor+1];
		for (int i = 2; i <= Math.sqrt(maxFactor); i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
			for (int j = i * i; j <= maxFactor; j += i) {
				noPrimes[j] = true;
			}
		}
		for (int i = (int) (Math.sqrt(maxFactor)+1); i <= maxFactor; i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
		}
		for (int i=primeIndex; i < primes.length; i++) {
			primes[i] = Integer.MAX_VALUE;
		}

		System.out.println("Prime table built max factor '" + maxFactor + "'       bytes used : " + primeIndex * 12);
	}


	public TrialDoubleFact(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		//		initPrimes();
		initPrimesEratosthenes();
	}


	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		double nD = n;
		for (int primeIndex = 1; primes[primeIndex] <= maxFactor; primeIndex++) {
			double nDivPrime = nD * primesInv[primeIndex];
			if (primes[primeIndex] == 0)
				System.out.println();
			// TODO choose the precision factor with respect to the maxFactor, if we are close to 52 bits
			long nDivPrimeLong;
			while (((nDivPrimeLong= (long)nDivPrime) - nDivPrime) >=  -0.001 &&
					nDivPrime - nDivPrimeLong < 0.001 &&
					//			while (Math.abs(nDivPrime - (((long)nDivPrime))) <=  0.001 &&
					nD > 1.0 &&
					nD % primes[primeIndex] == 0) {
				if (primeFactors == null)
					return primes[primeIndex];
				primeFactors.add((long) primes[primeIndex]);
				nD = Math.round(nDivPrime);
				//				// if the remainder n is lower then the maximal prime factor and it can not be split it must also
				//				// be prime factor
				//				if (n < maxFactor && n*n > maxFactor) {
				//					primeFactors.add(n);
				//					return 1;
				//				}

				nDivPrime = nD*primesInv[primeIndex];
			}
		}
		return (long) nD;
	}


	@Override
	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}
}
