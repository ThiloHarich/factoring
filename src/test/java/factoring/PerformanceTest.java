package factoring;

import java.math.BigInteger;
import java.util.Random;

import org.junit.Test;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.factor.squfof.SquFoF63;
import factoring.fermat.lehman.LehmanFactorFinder;
import factoring.fermat.lehman.LehmanFactorFinderRange;
import factoring.fermat.lehman.LehmanSimple;
import factoring.fermat.lehman.Lehman_FastOrig;
import factoring.rho.PollardRhoBrentDouble52;
import factoring.trial.TrialDoubleFact;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceTest {

	@Test
	public void testPerfHard(){
		// based on the the second biggest factor
		// 15 Bits TrialDoubleFact is the fastest
		// 31 Bits SquFoF31
		// 4x Bits Lehman

		final int bits = 40;
		final int numPrimes = 945;
		final int loop = 40;
		final int smallFactorBits = bits / 2;
		final long[] semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		//		final SquFoF63 factorizer1 = new SquFoF63();
		//		final Lehman factorizer3 = new Lehman(2);
		FactorizationOfLongs factorizer;
		//		factorizer = new LehmanFactorFinder(bits, 2.5f, false);
		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		//		final FactorizationOfLongs factorizer3 = new TrialInvFact(1 << (40/3));

		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		final long time1 = findFactors(factorizer, semiprimes, loop, 1l);
		findFactors(factorizer, semiprimes, loop, time1);
		//		findFactors(factorizer3, semiprimes, loop, time1);

		final long time2 = findFactors(factorizer, semiprimes, loop, time1);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time2);
		//		findFactors(factorizer3, semiprimes, loop, time2);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		final long time3 = findFactors(factorizer, semiprimes, loop, time2);
		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		//		factorizer = new LehmanFactorFinder(bits, 2.5f, false);

		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);

		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		factorizer = new LehmanFactorFinder(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time3);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);
		//		findFactors(factorizer3, semiprimes, loop, time3);

		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		//		factorizer = new LehmanFactorFinder(bits, 2.5f, false);

		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);

	}
	@Test
	public void testPerfHard2(){
		// based on the the second biggest factor
		// 15 Bits TrialDoubleFact is the fastest
		// 31 Bits SquFoF31
		// 4x Bits Lehman

		final int bits = 40;
		final int numPrimes = 945;
		final int loop = 40;
		final int smallFactorBits = bits / 2;
		final long[] semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		//		final SquFoF63 factorizer1 = new SquFoF63();
		//		final Lehman factorizer3 = new Lehman(2);
		FactorizationOfLongs factorizer;
		factorizer = new LehmanFactorFinder(bits, 2.5f, false);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		//		final FactorizationOfLongs factorizer3 = new TrialInvFact(1 << (40/3));

		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		final long time1 = findFactors(factorizer, semiprimes, loop, 1l);
		findFactors(factorizer, semiprimes, loop, time1);
		//		findFactors(factorizer3, semiprimes, loop, time1);

		final long time2 = findFactors(factorizer, semiprimes, loop, time1);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time2);
		//		findFactors(factorizer3, semiprimes, loop, time2);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		final long time3 = findFactors(factorizer, semiprimes, loop, time2);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		factorizer = new LehmanFactorFinder(bits, 2.5f, false);

		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);

		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		//		factorizer = new LehmanFactorFinder(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time3);
		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);
		//		findFactors(factorizer3, semiprimes, loop, time3);

		//		factorizer = new LehmanFactorFinderRange(bits, 2.5f, false);
		factorizer = new LehmanFactorFinder(bits, 2.5f, false);

		findFactors(factorizer, semiprimes, loop, time3);
		findFactors(factorizer, semiprimes, loop, time3);

	}

	@Test
	public void testPerfOneThird(){
		final int bits = 50;
		final int numPrimes = 16;
		final int loop = 500;
		final int smallFactorBits = bits / 3+2;
		final long[] semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		final SquFoF63 factorizer1 = new SquFoF63();
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorization(bits, .25f);
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorization(bits, .25f);
		final FactorizationOfLongs factorizer2 = new TrialDoubleFact(1 << smallFactorBits);
		factorizer2.setMaxFactor(1 << smallFactorBits);
		//        final FactorizationOfLongs factorizer2 = new LehmanYafuFact(5.8f);
		//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, .34f);
		final FactorizationOfLongs factorizer3 = new PollardRhoBrentDouble52();

		final long time1 = findFactors(factorizer1, semiprimes, loop, 1l);
		findFactors(factorizer2, semiprimes, loop, time1);
		findFactors(factorizer3, semiprimes, loop, time1);

		final long time2 = findFactors(factorizer1, semiprimes, loop, time1);
		findFactors(factorizer2, semiprimes, loop, time2);
		findFactors(factorizer3, semiprimes, loop, time2);

		final long time3 = findFactors(factorizer1, semiprimes, loop, time2);
		findFactors(factorizer2, semiprimes, loop, time3);
		findFactors(factorizer3, semiprimes, loop, time3);

		final long time4 = findFactors(factorizer1, semiprimes, loop, time3);
		findFactors(factorizer2, semiprimes, loop, time4);
		findFactors(factorizer3, semiprimes, loop, time4);

		final long time5 = findFactors(factorizer1, semiprimes, loop, time3);
		findFactors(factorizer2, semiprimes, loop, time4);
		findFactors(factorizer3, semiprimes, loop, time4);


	}

	protected long findFactors(final FactorizationOfLongs factorizer1, final long[] semiprimes, int loop, long time1) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				//				factorizer1.findFactors(semiprime, null);
				factorizer1.factorization(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-50s", factorizer1);
		System.out.println(name + " :    \t" + (0.0 + time)/time1 + "\t" + time);

		return time;
	}

	protected long findFactors(final FactorAlgorithmBase factorizer1, final long[] semiprimes, int loop, long time1) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findSingleFactor(BigInteger.valueOf(semiprime));
				//                factorizer1.factor(BigInteger.valueOf(semiprime));
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-45s", factorizer1.getClass().getSimpleName());
		System.out.println(name + " :    \t" +  (0.0 + time)/time1);
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


	@Test
	public void testPerfRandom()
	{
		final int bits = 25;
		//		final int bits = 35;
		final int range = 55000;


		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits/2);
		//		final Factorizer factorizer1 = new LehmanSmallRangeFact(bits, 1);
		//		final FactorizationOfLongs factorizer1 = new LehmanNoSqrtFact(bits,1f);
		//		final FactorizationOfLongs factorizer3 = new LehmanFactorFinderMod12(bits, 1f, true);
		//		Factorizer factorizer1 = new HartFact();
		//		Factorizer factorizer2 = new FermatResiduesRec();
		//				final FactorizationOfLongs factorizer2 = new TrialInvFact(1 << bits/2);
		//		Factorizer factorizer2 = new FermatFact();
		//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, 1f);
		final FactorizationOfLongs factorizer1 = new LehmanSimple(false);
		//		final FactorizationOfLongs factorizer2 = new LehmanPowFactorization(bits, 0f);
		//		final Factorizer factorizer2 = new TrialWithPrimesFact();
		//		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1, false);
		//		final FactorizationOfLongs factorizer3 = new PollardRhoBrentMultiGcd();
		final FactorAlgorithmBase factorizer2 = new Lehman_FastOrig(false);
		//		final FactorizationOfLongs factorizer3 = new PollardRhoBrentDouble();

		//		((TrialFactMod)factorizer1).setLimit(1 << 16);

		final long time1 = getFactors(factorizer1, bits, range, 1l);
		getFactors(factorizer2, bits, range, time1);
		//		getFactors(factorizer3, bits, range, time1);

		//		assertEquals(factors, factors2);

		final long time2 = getFactors(factorizer1, bits,range, time1);
		getFactors(factorizer2, bits, range, time2);
		//		getFactors(factorizer3, bits, range, time2);

		//		assertEquals(factors3, factors4);


		final long time3 = getFactors(factorizer1, bits,range, time1);
		getFactors(factorizer2, bits, range, time3);
		//		getFactors(factorizer3, bits, range, time3);


		final long time4 = getFactors(factorizer1, bits,range, time1);
		getFactors(factorizer2, bits, range, time4);
		//		getFactors(factorizer3, bits, range, time4);

		//		assertEquals(factors5, factors6);
	}


	public long getFactors(FactorizationOfLongs factorizer, int bits, int range, long time1) {

		// warmup
		final int factors = 0;
		final long begin = (1l << bits) +1584;
		final long start = System.nanoTime();
		for (long i = begin; i < begin + range; i++)
		{
			factorizer.factorization(i).size();
			//			factorizer.factorizationByFactors(i).size();
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-25s", factorizer.getClass().getSimpleName());
		System.out.println(name + " :    \t" + (0.0 + time)/time1);
		return time;
	}
	public long getFactors(FactorAlgorithm factorizer, int bits, int range, long time1) {

		// warmup
		final int factors = 0;
		final long begin = (1l << bits) +1584;
		final long start = System.nanoTime();
		for (long i = begin; i < begin + range; i++)
		{
			factorizer.factor(BigInteger.valueOf(i)).size();
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-20s", factorizer.getClass().getSimpleName());
		System.out.println(name + " :    \t" + (0.0 + time)/time1);
		return time;
	}

}
