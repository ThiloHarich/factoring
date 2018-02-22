package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import factoring.fermat.FermatFact;
import factoring.fermat.lehman.LehmanNoSqrtFact;
import factoring.fermat.lehman.LehmanNoSqrtFact2;

public class ErrorShiftTest {

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
		final int bits = 16;

		//		Factorizer factorizer1 = new LehmanResidueFact();
		//		Factorizer factorizer1 = new Fermat24();
		//		Factorizer factorizer1 = new HartFloorFact();
		final Factorizer factorizer2 = new LehmanNoSqrtFact2();
		final Factorizer factorizer1 = new LehmanNoSqrtFact();
		//		Factorizer factorizer1 = new LehmanYafuFact();
		//		Factorizer factorizer2 = new Lehman8kFirstFact();

		//		for (int i = 65538; i < 1 << (bits + 1); i++)
		final int exp = 16;
		long begin = (1L << exp);  // = 2^4 * 3^2 * 5
		begin = 3337L	; // * 23
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
	public void testPerf()
	{
		final int bits = 42;
		//		int bits = 25;


		//		Factorizer factorizer1 = new Lehman8kFirstFact();
		//		Factorizer factorizer2 = new LehmanResidueFact();
		//		final Factorizer factorizer2 = new LehmanRange1Fact();
		//		Factorizer factorizer1 = new HartFact();
		//		Factorizer factorizer2 = new FermatResiduesRec();
		//		Factorizer factorizer2 = new FermatResiduesSieve();
		//		Factorizer factorizer2 = new FermatFact();
		final Factorizer factorizer1 = new LehmanNoSqrtFact();
		//		Factorizer factorizer2 = new LehmanYafuFact();
		final Factorizer factorizer2 = new LehmanNoSqrtFact2();
		//		Factorizer factorizer1 = new LehmanSquaresFact();

		//		((TrialFactMod)factorizer1).setLimit(1 << 16);

		final int factors = getFactors(factorizer1, bits);
		final int factors2 =  getFactors(factorizer2, bits);

		assertEquals(factors, factors2);

		final int factors3 = getFactors(factorizer1, bits);
		final int factors4 =  getFactors(factorizer2, bits);

		//		assertEquals(factors3, factors4);

		final int factors5 = getFactors(factorizer1, bits);
		final int factors6 =  getFactors(factorizer2, bits);
		final int factors7 = getFactors(factorizer1, bits);
		final int factors8 =  getFactors(factorizer2, bits);

		//		assertEquals(factors5, factors6);
	}

	public int getFactors(Factorizer factorizer, int bits) {

		// warmup
		final int factors = 0;
		final long begin = (1l << bits) +1584;
		final int range = 8000;
		final long start = System.nanoTime();
		for (long i = begin; i < begin + range; i++)
		{
			factorizer.findAllPrimeFactors(i).size();
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
