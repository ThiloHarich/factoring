package factoring;

import java.math.BigInteger;
import java.util.Random;

import org.junit.Test;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.factor.lehman.Lehman_Fast;
import factoring.fermat.lehman.Lehman_Fast30;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard {

	final static int bits = 50;
	final static int numPrimes = 1845;
	final static int loop = 20;
	final static int smallFactorBits = bits / 2;
	static long[] semiprimes;

	public static void main(String[] args) {
		//		final FactorAlgorithmBase factorizer2 = new factoring.fermat.lehman.Lehman_Fast(false);
		//		final FactorAlgorithmBase factorizer2 = new SquFoF31();
		//		final FactorizationOfLongs factorizer1 = new TrialDoubleFact(1 << (bits/2));
		//		final FactorAlgorithmBase factorizer1 = new PollardRhoBrentMontgomery63();
		//		final FactorAlgorithmBase factorizer2 = new Lehman_Fast24_4(true);
		//		final FactorizationOfLongs factorizer1 = new PollardRhoBrentDouble53();
		final FactorAlgorithmBase factorizer2 = new Lehman_Fast30(true);
		//		final FactorAlgorithmBase factorizer1 = new Lehman_Fast6(true);
		//		final FactorizationOfLongs factorizer2 = new TrialInvFact(1 << (bits/2));
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorFinder(bits, 2.f, false);
		//		final FactorizationOfLongs factorizer2 = new PollardRhoBrentDouble52();
		final FactorAlgorithmBase factorizer1 = new Lehman_Fast(true);
		semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);
		test2(factorizer1);

		//		System.out.println("loop 6k      first : " + factorizer1.loop_6_1);
		//		System.out.println("loop 6      ground : " + factorizer1.loop_ground);
		//		System.out.println("loop 6k     second : " + factorizer1.loop_6_2);
		//		System.out.println("loop 6k + 3 ground : " + factorizer1.loop_3);

		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, loop);

		test2(factorizer2);
		//		System.out.println("loop 6k      first : " + factorizer2.loop_6_1);
		//		System.out.println("loop 6      ground : " + factorizer2.loop_ground);
		//		System.out.println("loop 6k     second : " + factorizer2.loop_6_2);
		//		System.out.println("loop 6k + 3 ground : " + factorizer2.loop_3);
	}

	public static void test2(FactorAlgorithmBase factorizer) {
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
	}

	@Test
	public static void test2(FactorizationOfLongs factorizer){
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
	}

	protected static long findFactors(final FactorAlgorithmBase factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findSingleFactor(BigInteger.valueOf(semiprime));
				//                factorizer1.factor(BigInteger.valueOf(semiprime));
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-45s", factorizer1.getClass().getSimpleName());
		System.out.println(name + " :    \t" +  time);
		return time;
	}

	protected static long findFactors(final FactorizationOfLongs factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				//				factorizer1.findFactors(semiprime, null);
				factorizer1.factorization(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-50s", factorizer1);
		System.out.println(name + " :    \t" + time);

		return time;
	}

	public static long[] makeSemiPrimesList(int bits, int smallFactorBits, int numPrimes) {
		final long[] semiPrimes = new long[numPrimes];
		for (int i=0; i< numPrimes; i++)
		{
			final Random rnd = new Random();
			final BigInteger fact1 = BigInteger.probablePrime(smallFactorBits, rnd);
			final BigInteger fact2 = BigInteger.probablePrime(bits - smallFactorBits, rnd);
			semiPrimes[i] = fact1.longValue() * fact2.longValue();
		}

		return semiPrimes;
	}




}
