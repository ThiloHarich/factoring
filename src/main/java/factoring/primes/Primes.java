package factoring.primes;

import java.io.*;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;

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

	public static long[] makeSemiPrimesList(int bits, int numPrimes, boolean readFromFile) {
		long[] semiPrimes = new long[numPrimes];
		final String file = "semiPrimes" + bits + ".dat";
		final Path path = Paths.get(file);

		if (Files.exists(path) && readFromFile)
		try {
			ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream((file)));
			semiPrimes = (long[]) inputStream.readObject();

//			List<String> semiPrimeList = Files.readAllLines(path);
//			System.out.println("found " + semiPrimes.length + " semi primes in file "+ path);
			return semiPrimes;
		} catch (final IOException | ClassNotFoundException e) {
			e.printStackTrace();
		}
		System.out.println("no semi primes file "+ path + " will create at least " + numPrimes + " semi primes");
		long start = System.currentTimeMillis();

		Random rnd = new Random();
		for (int i=0; i< numPrimes;)
		{
//			final int smallFactorBitsMin = (int) Math.ceil(bits * .37);
			final int smallFactorBitsMin = (int) (bits * .5);
			final int smallFactorBitsMax = (int) (bits * .5);
//			final int smallFactorBitsMax = 20;
//			final int smallFactorBitsMin = smallFactorBitsMax;
//			final int smallFactorBitsMin = 14;
			for(int bitsFirst = smallFactorBitsMin; bitsFirst <= smallFactorBitsMax && i< numPrimes; bitsFirst++) {
//						final int smallFactorBits = (bits / 3) - 1;

//			rnd = new Random();
				final BigInteger fact1 = BigInteger.probablePrime(bitsFirst, rnd);
				final int bigFactorBits = bits - bitsFirst;
//			rnd = new Random();
				final BigInteger fact2 = BigInteger.probablePrime(bigFactorBits, rnd);
				semiPrimes[i] = fact1.longValue() * fact2.longValue();
				i++;
				if (i % 10000 == 0)
				System.out.println((100.0 * i / numPrimes) + "% of semi prime creation done. Created " + i + " semi primes");
			}
		}
//		String semiprimeList = "";
//		for (Long semiPrime : semiPrimes) {
//			if (semiPrime != null)
//			semiprimeList += semiPrime + System.lineSeparator();
//		}
		long endCreation = System.currentTimeMillis();
//		final byte[] strToBytes = semiprimeList.getBytes();


		try {
			ObjectOutputStream outputStream = new ObjectOutputStream(new FileOutputStream(file));
			outputStream.writeObject(semiPrimes);
//			Files.write(path, semiPrimes);
		} catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long endWrite = System.currentTimeMillis();
		System.out.println("wrote " + semiPrimes.length + " semi primes. Took " + (endCreation - start)/ 1000.0 + " sec to create numbers.");
		System.out.println("wrote " + semiPrimes.length + " semi primes. Took " + (endWrite - endCreation)/ 1000.0 + " sec to write numbers.");

		return semiPrimes;
	}
}
