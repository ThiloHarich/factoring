package factoring;

import java.math.BigInteger;
import java.util.Random;

import factoring.fermat.lehman.*;
import org.junit.Test;

import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
import de.tilman_neumann.math.factor.FactorAlgorithm;
import de.tilman_neumann.math.factor.SingleFactorFinder;

public class PerformanceTest {

	@Test
	public void testPerfHard(){
		final int bits = 34;
		final int numPrimes = 2000;
		int loop = 200;
		final long[] semiprimes = makeSemiPrimesList(bits, bits/2, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		final SingleLongFactorFinder factorizer3 = new LehmanFactorization(bits, .5f);
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));
		final SingleLongFactorFinder factorizer2 = new LehmanFactorFinderStep8Mod3(bits,3.2f);
//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, .34f);
		final SingleFactorFinder factorizer1 = new CombinedFactorAlgorithm(1);

		findFactors(factorizer1, semiprimes, loop);
		findFactors(factorizer2, semiprimes, loop);
		findFactors(factorizer3, semiprimes, loop);

		findFactors(factorizer1, semiprimes, loop);
		findFactors(factorizer2, semiprimes, loop);
		findFactors(factorizer3, semiprimes, loop);

		findFactors(factorizer1, semiprimes, loop);
		findFactors(factorizer2, semiprimes, loop);
		findFactors(factorizer3, semiprimes, loop);

		findFactors(factorizer1, semiprimes, loop);
		findFactors(factorizer2, semiprimes, loop);
		findFactors(factorizer3, semiprimes, loop);


	}

	protected void findFactors(final SingleFactorFinder factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findSingleFactor(BigInteger.valueOf(semiprime));
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-45s", factorizer1.getName());
		System.out.println(name + " :    \t" + (time));
	}

	protected void findFactors(final SingleLongFactorFinder factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findSingleFactor(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-25s", factorizer1.getClass().getSimpleName());
		System.out.println(name + " :    \t" + (time));
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
		final int range = 130000;


		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits/2);
//		final Factorizer factorizer1 = new LehmanSmallRangeFact(bits, 1);
				final FactorizationOfLongs factorizer1 = new LehmanNoSqrtFact(bits,1f);
		//		Factorizer factorizer1 = new HartFact();
		//		Factorizer factorizer2 = new FermatResiduesRec();
		//		final Factorizer factorizer2 = new TrialInvFact(1 << bits/2);
		//		Factorizer factorizer2 = new FermatFact();
//		final SingleFactorFinder factorizer1 = new LehmanFactorization(bits, 1f);
		final FactorizationOfLongs factorizer3 = new LehmanFactorization(bits, 0f);
		final FactorizationOfLongs factorizer2 = new LehmanPowFactorization(bits, 0f);
		//		final Factorizer factorizer2 = new TrialWithPrimesFact();
//		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1);
//				final Factorizer factorizer1 = new LehmanYafuFact(2.8f);

		//		((TrialFactMod)factorizer1).setLimit(1 << 16);

		final int factors = getFactors(factorizer1, bits, range);
		final int factors2 =  getFactors(factorizer2, bits, range);
				final int factors13 =  getFactors(factorizer3, bits, range);

		//		assertEquals(factors, factors2);

		final int factors3 = getFactors(factorizer1, bits,range);
		final int factors4 =  getFactors(factorizer2, bits, range);
				final int factors14 =  getFactors(factorizer3, bits, range);

		//		assertEquals(factors3, factors4);

		final int factors5 = getFactors(factorizer1, bits, range);
		final int factors6 =  getFactors(factorizer2, bits, range);
				final int factors15 =  getFactors(factorizer3, bits, range);

		final int factors7 = getFactors(factorizer1, bits, range);
		final int factors8 =  getFactors(factorizer2, bits, range);
				final int factors9 =  getFactors(factorizer3, bits, range);

		//		assertEquals(factors5, factors6);
	}


    public int getFactors(FactorizationOfLongs factorizer, int bits, int range) {

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
        System.out.println(name + " :    \t" + (time));
        return factors;
    }
    public int getFactors(FactorAlgorithm factorizer, int bits, int range) {

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
        System.out.println(name + " :    \t" + (time));
        return factors;
    }

}
