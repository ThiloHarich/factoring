package factoring.trial;

import java.util.Collection;

/**
 * This implementation is generating a list of all primes up to a limit.
 * Beside of storing the prime itself, it also store the reciprocal value.
 * When checking a number if it is dividable by the prime, we will not divide
 * the number by the prime, we will multiply by the inverse, since this is faster.
 * Due to precision we have to check for a given range near an Integer.
 * And then do a Long division.
 * This implementation is around two times faster then a version  based on long numbers.
 * Since Double only has 52 bis for the remainder, this can only work for numbers below 2^52.
 * We can only factorize numbers up to maxFactor^2
 * When calling it with bigger numbers only prime factors below
 * maxFactor were added to the factors. {@link #findPrimeFactors(long, Collection, double, double)} then might return a
 * composite number.
 **
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialRangeFact {

	// The number of values to be printed
	private static final int PRINT_NUM = 20000;
	// for printing we need a value a little bit above 1
	public static final double PRINT_CONST = 1.0000001;
	private int maxFactor = 65535;
	double[] primesInv;
	int[] primes;
	int maxFactorIndex = 0;



	/**
	 * finds the prime factors up to maxFactor by the sieve of eratosthenes.
	 * Not optimized, since this is only called once when initializing.
	 */
	void initPrimesEratosthenes()
	{
		final double logMaxFactor = Math.log(maxFactor);
		final int maxPrimeIndex = (int) ((maxFactor) / (logMaxFactor - 1.1)) + 1;
		primesInv = new double [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		primes = new int [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		final boolean [] noPrimes = new boolean [maxFactor];
		for (int i = 2; i <= Math.sqrt(maxFactor); i++) {
			if (!noPrimes[i]) {
				primes[maxFactorIndex] = i;
				primesInv[maxFactorIndex++] = 1.0 / i;
			}
			for (int j = i * i; j < maxFactor; j += i) {
				noPrimes[j] = true;
			}
		}
		for (int i = (int) (Math.sqrt(maxFactor)+1); i < maxFactor; i++) {
			if (!noPrimes[i]) {
				primes[maxFactorIndex] = i;
				primesInv[maxFactorIndex++] = 1.0 / i;
			}
		}
		for (int i = maxFactorIndex; i < primes.length; i++) {
			primes[i] = Integer.MAX_VALUE;
		}

		System.out.println("Prime table built max factor '" + maxFactor + "'       bytes used : " + maxFactorIndex * 12);
	}


	public TrialRangeFact(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		//		initPrimes();
		initPrimesEratosthenes();
	}


	/**
	 * since
	 * f(x) = x * ln (x)
	 * f(x*b) = x*b * ln (x*b) = x*b * (ln (x) + ln(b)) = f(x)*b + x*b*ln(b) < =  f(x)*b * (1 +  b*ln(b)/ln(x))
	 * b = begin - begin * Math.
	 * begin - 1 =
	 * primes [(int)(maxFactorIndex*begin)] < primes [maxFactorIndex]
	 * @param n
	 * @param primeFactors
	 * @param begin
	 * @param end
	 * @return
	 */
	public long findPrimeFactors(long n, Collection<Long> primeFactors, double begin, double end) {

		// adjust the begin and end value such that it begins/ends just below n^1/3
		// this ensures we always have prime factors in the lehman phase as well
//		double b = end == 1 ? begin : end;
//		double correctFakt = b - b * Math.log(b) * 3 / Math.log(n);
//		begin = begin == 0 ? 0 : correctFakt;
//		end = end == 1 ? 1 : correctFakt;
		for (int primeIndex = (int)(maxFactorIndex*begin); primeIndex < maxFactorIndex*end; primeIndex++) {
			double nDivPrime = n*primesInv[primeIndex];
			// TODO choose the precision factor with respect to the maxFactor!?
			if (primes[primeIndex] == 0)
				System.out.println();
			while (Math.abs(Math.round(nDivPrime) - nDivPrime) < 0.01 && n > 1 && n % primes[primeIndex] == 0) {
				if (primeFactors == null)
					return (long) primes[primeIndex];
				primeFactors.add((long) primes[primeIndex]);
				n = Math.round(nDivPrime);
				nDivPrime = n*primesInv[primeIndex];
			}
		}
		return n;
	}

	public long findPrimeFactor(long n, Collection<Long> factors, int primeIndex) {
		double nDivPrime = n*primesInv[primeIndex];
		// TODO choose the precision factor with respect to the maxFactor!?
		if (primes[primeIndex] == 0)
            System.out.println();
		while (Math.abs(Math.round(nDivPrime) - nDivPrime) < 0.01 && n > 1 && n % primes[primeIndex] == 0) {
            factors.add((long) primes[primeIndex]);
            n = Math.round(nDivPrime);
            nDivPrime = n*primesInv[primeIndex];
        }
		return n;
	}


	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}
}
