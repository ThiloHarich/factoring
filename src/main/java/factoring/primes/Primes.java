package factoring.primes;

public class Primes {

	public int[] primes;
	public double[] primesInv;

	/**
	 * finds the prime factors up to maxFactor by the sieve of eratosthenes.
	 * Not optimized, since this is only called once when initializing.
	 * @param primesInv
	 * @param primes
	 * @param maxFactor
	 */
	public static Primes initPrimesEratosthenes(int maxFactor)
	{
		final Primes primes = new Primes();
		final double logMaxFactor = Math.log(maxFactor);
		final int maxPrimeIndex = (int) Math.max(100, ((maxFactor) / (logMaxFactor - 1.1)) + 4);
		primes.primesInv = new double [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		primes.primes = new int [maxPrimeIndex]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end
		int primeIndex = 0;
		final boolean [] noPrimes = new boolean [maxFactor];
		for (int i = 2; i <= Math.sqrt(maxFactor); i++) {
			if (!noPrimes[i]) {
				primes.primes[primeIndex] = i;
				primes.primesInv[primeIndex++] = 1.0 / i;
			}
			for (int j = i * i; j < maxFactor; j += i) {
				noPrimes[j] = true;
			}
		}
		for (int i = (int) (Math.sqrt(maxFactor)+1); i < maxFactor; i++) {
			if (!noPrimes[i]) {
				primes.primes[primeIndex] = i;
				primes.primesInv[primeIndex++] = 1.0 / i;
			}
		}
		for (int i=primeIndex; i < primes.primes.length; i++) {
			primes.primes[i] = Integer.MAX_VALUE;
		}
		return primes;
	}

}
