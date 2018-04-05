package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;

import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
import de.tilman_neumann.math.factor.FactorAlgorithm;
import factoring.fermat.lehman.LehmanReverseFact;
import org.junit.Test;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import factoring.fermat.FermatFact;
import factoring.fermat.lehman.LehmanNoSqrtFact;
import factoring.fermat.lehman.LehmanSmallRangeFact;
import factoring.fermat.lehman.LehmanYafuFact;

public class CorrectnessTest {

	@Test
	public void test10037_4339 ()
	{
		final int p = 10037;
		final int q = 4339;
		final List<Long> factors = factors(p, q);
		final ErrorShiftFact fact = new ErrorShiftFact(true);
		final long factor = fact.findFactor(p, q);

		assertTrue(factors.contains(factor));
	}

	List<Integer> qs ()
	{
		return Arrays.asList(20023, 20029, 24043, 24029, 30029, 45007, 60029, 60017, 60013, 90023, 90007);
	}

	@Test
	public void testCorrect() {
		final int bits = 8;
		final int bitsMax = 32;

		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits + 4);
		//		Factorizer factorizer1 = new Fermat24();
		//		Factorizer factorizer1 = new LehmanBigFact(bitsMax, 1);
		//		final Factorizer factorizer2 = new LehmanMod16Fact(bitsMax);
		//		final Factorizer factorizer2 = new LehmanApproxFact();
		final Factorizer factorizer2 = new LehmanNoSqrtFact(bitsMax, 1.01f);
		final Factorizer factorizer1 = new LehmanYafuFact(1);
		//		final Factorizer factorizer2 = new TrialInvFact(1 << bits + 4);

		//		for (int i = 65538; i < 1 << (bits + 1); i++)
		long begin = (1L << bits) +1;  // = 2^4 * 3^2 * 5
		begin = 315l	; // * 23
		// 29*23 * 53
		// 29*53 * 23 ->
		while (begin < Long.MAX_VALUE / 1000)
		{
			for (long i = begin; i < begin + begin/8; i++) {
				final Collection<Long> factors = factorizer1.findAllPrimeFactors(i);
				System.out.println(i + ": " + factorizer1.printFactors(i));
				//			Collection<Integer> factors = factorizer1.findAllPrimeFactors(i);
				final Collection<Long> factors2 = factorizer2.findAllPrimeFactors(i);
				System.out.println(i + ": " + factorizer2.printFactors(i));

				assertEquals("Test failed for " + i, factors.size(), factors2.size());
				final Multiset<Long> factorsSet = TreeMultiset.create();
				factorsSet.addAll(factors2);

				for (final long factor :
					factors) {
					assertTrue("Test failed for " + i, factorsSet.contains(factor));
				}
			}
			begin = begin + begin;
		}

	}
	@Test
	public void testPerfHard(){
		final int bits = 40;
		final int numPrimes = 1000;
		int loop = 30;
		final long[] semiprimes = makeSemiPrimesList(bits, bits/3+1, numPrimes);

		System.out.println("finished making hard numbers");
		final long start = System.currentTimeMillis();
		final Factorizer factorizer2 = new LehmanNoSqrtFact(bits, 3f);
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));
//		final Factorizer factorizer1 = new LehmanYafuFact(2.8f);
		final Factorizer factorizer1 = new LehmanReverseFact(bits, 1.5f);
		final FactorAlgorithm factorizer3 = new CombinedFactorAlgorithm(1);

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

	protected void findFactors(final FactorAlgorithm factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.factor(BigInteger.valueOf(semiprime)).size();
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-25s", factorizer1.getClass().getSimpleName());
		System.out.println(name + " :    \t" + (time));
	}

	protected void findFactors(final Factorizer factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findAllPrimeFactors(semiprime);
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
		final int bits = 40;
		//		final int bits = 35;
		final int range = 40000;


		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits/2);
//		final Factorizer factorizer1 = new LehmanSmallRangeFact(bits, 1);
		//		final Factorizer factorizer2 = new LehmanRange1Fact();
		//		Factorizer factorizer1 = new HartFact();
		//		Factorizer factorizer2 = new FermatResiduesRec();
		//		final Factorizer factorizer2 = new TrialInvFact(1 << bits/2);
		//		Factorizer factorizer2 = new FermatFact();
		final Factorizer factorizer2 = new LehmanNoSqrtFact(bits, 1);
		final Factorizer factorizer1 = new LehmanNoSqrtFact(bits, 0.2f);
		//		final Factorizer factorizer2 = new TrialWithPrimesFact();
						final FactorAlgorithm factorizer3 = new CombinedFactorAlgorithm(1);
		//		final Factorizer factorizer1 = new LehmanYafuFact();

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

	public int getFactors(Factorizer factorizer, int bits, int range) {

		// warmup
		int factors = 0;
		final long begin = (1l << bits) +1584;
		final long start = System.nanoTime();
		for (long i = begin; i < begin + range; i++)
		{
			factors += factorizer.findAllPrimeFactors(i).size();
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-20s", factorizer);
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


	private void checkFactors(int p) {
		final int operations = 0;
		final int fermatOperations = 0;
		for (final int q : qs())
		{
			final List<Long> factors = factors(p, q);
			//			final Factorizer fact = new ErrorShiftFact(true);
			final FermatFact fact = new FermatFact();
			//			final ErrorShift2DFact fact = new ErrorShift2DFact();
			//		final ErrorYIncShiftFact fact = new ErrorYIncShiftFact();
			fact.findFactors(p*q, factors);
			final long factor = factors.get(0);
			//			assertTrue(factors.contains(factor));
			if (!factors.contains(factor))
				System.err.println("Factor not found");
			//			operations += fact.operations;
			//			fermatOperations += fact.getOperationsFermat();
		}
		System.out.println("Overall Speedup " + (fermatOperations + 0.0)/operations);
	}

	@Test
	public void test9721 ()
	{
		final int p = 9721;
		checkFactors(p);
	}

	@Test
	public void test9719 ()
	{
		final int p = 9719;
		checkFactors(p);
	}
	@Test
	public void test11003 ()
	{
		final int p = 11003;
		checkFactors(p);
	}
	@Test
	public void test11007 ()
	{
		final int p = 10037;
		checkFactors(p);
	}


	private List<Long> factors(long p, long q) {
		return Arrays.asList(p,q);
	}

}