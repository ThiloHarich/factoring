package factoring.trial;

import java.util.Collection;

import factoring.FactorizationOfLongs;

/**
 * This implementation is generating a array of all number up to a limit.
 * It stores the lowest prime factor of the number.
 * We only consider numbers which have a factor > 11.
 * This is we only have to store numbers with
 * 1,13, 169=13*13 mod 210
 *
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class SmallestPrimeFactorTable implements FactorizationOfLongs {

	// The number of values to be printed
	private static final int PRINT_NUM = 20000;
	// for printing we need a value a little bit above 1
	public static final double PRINT_CONST = 1.0000001;
	private int maxFactor = 65535;
	int[] smallestFactor;


	/**
	 * finds the prime factors up to maxFactor by the sieve of eratosthenes.
	 * Not optimized, since this is only called once when initializing.
	 * We might sieve in buckets to reduce the amount of memory to be used.
	 */
	void initPrimesEratosthenes()
	{
		smallestFactor = new int [maxFactor];
		for (int i = 2; i <= Math.sqrt(maxFactor); i++) {
			for (int j = 2*i; j < maxFactor; j += i) {
				if (smallestFactor[i] == 0) {
					smallestFactor[j] = i;
				}
			}
		}
		System.out.println("Smalles factor table built max factor '" + maxFactor + "'       bytes used : " + maxFactor * 4);
	}


	public SmallestPrimeFactorTable(int maxFactor) {
		//        if (maxFactor > 65535)
		//            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
		this.maxFactor = maxFactor;
		//		initPrimes();
		initPrimesEratosthenes();
	}

	@Override
	public void setMaxFactor(int maxTrialFactor) {
		maxFactor = maxTrialFactor;
	}


	public long findPrimeFactor(long n) {
		if (n > maxFactor)
			return n;
		final int factor = smallestFactor[(int) n];
		if (factor == 0)
			return n;
		return factor;
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		while(true) {
			final long primeFactor = findPrimeFactor(n);
			primeFactors.add(primeFactor);
			if (primeFactor == n) {
				return n/primeFactor;
			}
			n = n/primeFactor;
		}
	}
}
