package factoring.trial;

import java.util.Collection;

import factoring.FindPrimeFact;

/**
 * This implementation is generating a list of all primes up to a limit.
 * Instead of storig the prime itself, it will store the reziprok value.
 * When checking a number if it is dividable by the prime, we will not divide
 * the number by the prime, we will multiply by the inverse, since this is faster.
 * number is divideable. Here we us a limit of 65536=2^16.
 * We can only factorize numbers up to 2^32.
 * When calling it with bigger numbers only prime factors below
 * 2^16 were added to the factors. {@link #findPrimeFactors(long, Collection)} then might return a
 * composite number. This strictly means it can not be applied for factoring big numbers,
 * but it can be applied for finding the prime factors below 2^16 and applying other algorithms
 * to factor the returned remainder. This is exactly what the lehman and hart algorithms need.
 * We have choosen 2^16 because when factoring long numbers by the lehman method, they have to be
 * lower n =  2^48 = (2^16)^3, which means we can find all primesInv below n^1/3 with this algorithm.
 *
 * Create a new instance, unless you want to have control over the search interval by using {@link #setMaxFactor(int)}.
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialInvFact extends FindPrimeFact {

	private static final int PRINT_NUM = 20000;
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
			System.out.printf("%d,", (int)(1.0000001 / primesInv[i]));
		}
		if (maxPrimeIndex > PRINT_NUM)
			System.out.printf(" ,..., %f,%f%n", 1.001 / primesInv[k-2], 1.0 / primesInv[k-1]);
		else
			System.out.println();
	}


	public TrialInvFact(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		initPrimes();
	}


	@Override
	public long findPrimeFactors(long n, Collection<Long> factors) {
		// TODO fill the end of the array with Integer.MaxValue
		//		final double maxFactorInv = 1.0 / maxFactor;
		for (int primeIndex = 0; primes[primeIndex] <= maxFactor && primesInv[primeIndex] > 0; primeIndex++) {
			double nDivPrime = n*primesInv[primeIndex];
			//			long nDivPrimeLong = (long) (nDivPrime + .0001);

			//			while (Math.abs(nDivPrimeLong - nDivPrime) < 0.0001 && n > 1 && n % primes[primeIndex]==0) {
			while (Math.round(nDivPrime) == nDivPrime && n > 1 && n % primes[primeIndex] == 0) {
				factors.add(Math.round(1.0 / primesInv[primeIndex]));
				n = Math.round(nDivPrime);
				nDivPrime = n*primesInv[primeIndex];
				//				nDivPrimeLong = (long) (n*primesInv[primeIndex]+ 0.0001);
			}
		}
		return n;
	}


	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}
}
