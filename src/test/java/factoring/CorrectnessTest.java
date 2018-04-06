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

	private List<Long> factors(long p, long q) {
		return Arrays.asList(p,q);
	}

}