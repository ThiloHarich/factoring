package factoring.trial;

import java.util.Collection;

import factoring.FindPrimeFact;

/**
 * This implementation is generating a list of all primes up to a limit and will then check if the
 * number is divideable. Here we us a limit of 65536=2^16.
 * We can only factorize numbers up to 2^32.
 * When calling it with bigger numbers only prime factors below
 * 2^16 were added to the factors. {@link #findPrimeFactors(long, Collection)} then might return a
 * composite number. This strictly means it can not be applied for factoring big numbers,
 * but it can be applied for finding the prime factors below 2^16 and applying other algorithms
 * to factor the returned remainder. This is exactly what the lehman and hart algorithms need.
 * We have choosen 2^16 because when factoring long numbers by the lehman method, they have to be
 * lower n =  2^48 = (2^16)^3, which means we can find all primes below n^1/3 with this algorithm.
 *
 * Create a new instance, unless you want to have control over the search interval by using {@link #setMaxFactor(int)}.
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialPrimesDynamicFact extends FindPrimeFact {

	private int maxFactor = 65535;
	int[] primes;
	// we keep the index so that we can continue with the factoring after applying
	// other algorithms like lehman, not so nice
	int primeIndex = 0;

	void initPrimes()
	{
		final double logMaxFactor = Math.log(maxFactor);
		final int maxPrimeIndex = (int) ((maxFactor) / (logMaxFactor - 1.1));
		primes = new int [maxPrimeIndex]; //the 6542 primes up to 65536=2^16, then sentinel 65535 at end
		// TODO we do not need 2
		primes[0]=2;
		primes[1]=3;
		primes[2]=5;
		int k = 3;
		for(int i=7; i<maxPrimeIndex; i+=2){
			boolean isPime = true;
			for(int j=0; primes[j]* primes[j] <= i && isPime; j++){
				if(i%primes[j]==0)
					isPime = false;
			}
			if (isPime) {
				primes[k] = i;
				k++;
			}
		}
		//		assert(k==6542);
		//		primes[k] = 65535; //sentinel
		System.out.printf("Prime table[0..%d] built: ", k);
		for(int i=0; i<20; i++)
		{
			System.out.printf("%d,", primes[i]);
		}
		System.out.printf("%d,...,%d,(%d)%n", primes[20],primes[maxPrimeIndex-2],primes[maxPrimeIndex-1]);
	}


	public TrialPrimesDynamicFact(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		initPrimes();
	}


	@Override
	public long findPrimeFactors(long n, Collection<Long> factors) {
		for (; primes[primeIndex] < maxFactor; primeIndex++) {
			while (n % primes[primeIndex] == 0) {
				factors.add((long)primes[primeIndex]);
				n /= primes[primeIndex];
			}
		}
		return n;
	}
}
