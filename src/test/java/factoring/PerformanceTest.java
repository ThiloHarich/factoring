package factoring;

import java.math.BigInteger;
import java.util.Random;

import de.tilman_neumann.math.factor.FactorAlgorithmBase;
import factoring.fermat.lehman.*;
import factoring.fermat.lehman.playground.LehmanFactorFinderMod36;
import org.junit.Test;

import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceTest {

	@Test
	public void testPerfHard(){
		final int bits = 40;
		final int numPrimes = 1340;
		int loop = 20;
        int smallFactorBits = bits / 2;
        final long[] semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		final FactorizationOfLongs factorizer1 = new LehmanFactorFinderMod36(bits, 5);
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));
        final FactorizationOfLongs factorizer2 = new LehmanFactorFinderMod12(bits,5f, false);
//        final FactorizationOfLongs factorizer2 = new LehmanYafuFact(5.8f);
//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, .34f);
		final FactorAlgorithmBase factorizer3 = new CombinedFactorAlgorithm(1);

		long time1 = findFactors(factorizer1, semiprimes, loop, 1l);
		findFactors(factorizer2, semiprimes, loop, time1);
		findFactors(factorizer3, semiprimes, loop, time1);

        long time2 = findFactors(factorizer1, semiprimes, loop, time1);
        findFactors(factorizer2, semiprimes, loop, time2);
        findFactors(factorizer3, semiprimes, loop, time2);

        long time3 = findFactors(factorizer1, semiprimes, loop, time2);
        findFactors(factorizer2, semiprimes, loop, time3);
        findFactors(factorizer3, semiprimes, loop, time3);

        long time4 = findFactors(factorizer1, semiprimes, loop, time3);
        findFactors(factorizer2, semiprimes, loop, time4);
        findFactors(factorizer3, semiprimes, loop, time4);


	}

	protected long findFactors(final FactorizationOfLongs factorizer1, final long[] semiprimes, int loop, long time1) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
                factorizer1.findSingleFactor(semiprime);
//                factorizer1.factorization(semiprime);
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
	private long[] makeSemiPrimesList(int bits, int smallFactorBits, int numPrimes) {
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
		final int bits = 30;
		//		final int bits = 35;
		final int range = 93000;


		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits/2);
//		final Factorizer factorizer1 = new LehmanSmallRangeFact(bits, 1);
//		final FactorizationOfLongs factorizer1 = new LehmanNoSqrtFact(bits,1f);
//		final FactorizationOfLongs factorizer3 = new LehmanFactorFinderMod12(bits, 1f, true);
		//		Factorizer factorizer1 = new HartFact();
		//		Factorizer factorizer2 = new FermatResiduesRec();
//				final FactorizationOfLongs factorizer2 = new TrialInvFact(1 << bits/2);
		//		Factorizer factorizer2 = new FermatFact();
//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, 1f);
        final FactorizationOfLongs factorizer1 = new LehmanFactorization(bits, 0f);
//		final FactorizationOfLongs factorizer2 = new LehmanPowFactorization(bits, 0f);
		//		final Factorizer factorizer2 = new TrialWithPrimesFact();
		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1);
		final FactorizationOfLongs factorizer3 = new LehmanYafuFact(1f);

		//		((TrialFactMod)factorizer1).setLimit(1 << 16);

		final long time1 = getFactors(factorizer1, bits, range, 1l);
		getFactors(factorizer2, bits, range, time1);
		getFactors(factorizer3, bits, range, time1);

		//		assertEquals(factors, factors2);

		final long time2 = getFactors(factorizer1, bits,range, time1);
		getFactors(factorizer2, bits, range, time2);
		getFactors(factorizer3, bits, range, time2);

		//		assertEquals(factors3, factors4);


        final long time3 = getFactors(factorizer1, bits,range, time1);
        getFactors(factorizer2, bits, range, time3);
        getFactors(factorizer3, bits, range, time3);


        final long time4 = getFactors(factorizer1, bits,range, time1);
        getFactors(factorizer2, bits, range, time4);
        getFactors(factorizer3, bits, range, time4);

		//		assertEquals(factors5, factors6);
	}


    public long getFactors(FactorizationOfLongs factorizer, int bits, int range, long time1) {

        // warmup
        final int factors = 0;
        final long begin = (1l << bits) +1584;
        final long start = System.nanoTime();
        for (long i = begin; i < begin + range; i++)
        {
			factorizer.factorizationByPrimes(i).size();
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
