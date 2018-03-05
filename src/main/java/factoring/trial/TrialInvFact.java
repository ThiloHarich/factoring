package factoring.trial;

import java.util.Collection;

import factoring.FindPrimeFact;

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
 * maxFactor were added to the factors. {@link #findPrimeFactors(long, Collection)} then might return a
 * composite number.
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialInvFact extends FindPrimeFact {

	// The number of values to be printed
	private static final int PRINT_NUM = 20000;
	// for printing we need a value a little bit above 1
	public static final double PRINT_CONST = 1.0000001;
	private int maxFactor = 65535;
	double[] primesInv;
	int[] primes;

	void initPrimes()
	{
		final double logMaxFactor = Math.log(maxFactor);
		final int maxPrimeIndex = (int) ((maxFactor) / (logMaxFactor - 1.1)) + 100;
		primesInv = new double [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		primes = new int [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		// TODO we do not need 2
		primesInv[0]= 1.0 / 2.0;
		primesInv[1]= 1.0 / 3.0;
		primesInv[2]= 1.0 / 5.0;
		primes[0]=2;
		primes[1]=3;
		primes[2]=5;
		int k = 3;
		for(int i=7; i<maxFactor; i+=2){
			boolean isPime = true;
			for(int j = 0; primes[j]* primes[j] <= i && isPime; j++){
				final double nDivPrime = i*primesInv[j];
				final long nDivPrimeLong = (long) (i*primesInv[j] + 0.0001);
				//				if (Math.round(nDivPrime) == nDivPrime && i % primes[j]==0) {
				if (Math.abs(nDivPrimeLong - nDivPrime) < 0.0001 && i % primes[j]==0) {
					//					if (Math.abs(Math.round(nDivPrime)-nDivPrime) < 0.001 && i % primes[j]==0) {
					isPime = false;
				}
			}
			if (isPime) {
				primesInv[k] = 1.0 / i;
				primes[k] = i;
				k++;
			}
		}
		//		assert(k==6542);
		//		primesInv[k] = 65535; //sentinel
		System.out.printf("Prime table[0..%d] built: ", k);
		for(int i=0; i<Math.min(PRINT_NUM,maxPrimeIndex) ; i++)
		{
			System.out.printf("%d,", (int)(PRINT_CONST / primesInv[i]));
		}
		if (maxPrimeIndex > PRINT_NUM)
			System.out.printf(" ,..., %f,%f%n", PRINT_CONST / primesInv[k-2], PRINT_CONST / primesInv[k-1]);
		else
			System.out.println();
	}

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
		int primeIndex = 0;
		boolean [] noPrimes = new boolean [maxFactor];
		for (int i = 2; i <= Math.sqrt(maxFactor); i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
			for (int j = i * i; j < maxFactor; j += i) {
				noPrimes[j] = true;
			}
		}
		for (int i = (int) (Math.sqrt(maxFactor)+1); i < maxFactor; i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
		}
		System.out.printf("Prime table[0..%d] built: ", primeIndex);
		for(int i=0; i<Math.min(PRINT_NUM,maxPrimeIndex) ; i++)
		{
			System.out.printf("%d,", (int)(PRINT_CONST / primesInv[i]));
		}
		if (maxPrimeIndex > PRINT_NUM)
			System.out.printf(" ,..., %f,%f%n", PRINT_CONST / primesInv[primeIndex-2], PRINT_CONST / primesInv[primeIndex-1]);
		else
			System.out.println();
		for (;primeIndex < maxPrimeIndex; primeIndex++)
		{
			primes[primeIndex] = Integer.MAX_VALUE;
//			primesInv[primeIndex++] = 1.0 / i;
		}
	}


	public TrialInvFact(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
//		initPrimes();
		initPrimesEratosthenes();
	}


	@Override
	public long findPrimeFactors(long n, Collection<Long> factors) {
		for (int primeIndex = 0; primes[primeIndex] <= maxFactor; primeIndex++) {
			double nDivPrime = n*primesInv[primeIndex];
			// TODO choose the precision factor with respect to the maxFactor!?
			while (Math.abs(Math.round(nDivPrime) - nDivPrime) < 0.01 && n > 1 && n % primes[primeIndex] == 0) {
				factors.add((long) primes[primeIndex]);
				n = Math.round(nDivPrime);
				nDivPrime = n*primesInv[primeIndex];
			}
		}
		return n;
	}


	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}
}
