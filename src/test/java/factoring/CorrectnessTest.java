package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.util.SortedMultiset;
import factoring.fermat.lehman.Lehman_TillSimple;
import factoring.fermat.lehman.Lehman_TillSimple3;
import factoring.rho.PollardRhoBrentDouble52;
import factoring.shift.ErrorShiftFact;
import factoring.trial.variant.TrialFact;

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
		final int bits = 40;

		long begin = (1L << bits) +1;
		begin = 5640012124823l	; // * 23
		//		final FactorAlgorithm factorizer1 = new Lehman_TillSimple4();
		final FactorAlgorithmBase factorizer1 = new Lehman_TillSimple(1);

		//		Factorizer factorizer1 = new Fermat24();
		//		Factorizer factorizer1 = new LehmanBigFact(bitsMax, 1);
		//		final Factorizer factorizer2 = new LehmanMod16Fact(bitsMax);
		//		final Factorizer factorizer2 = new LehmanApproxFact();
		final FactorAlgorithmBase factorizer2 = new Lehman_TillSimple3(true);
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinderRange(bits, 2f, true);
		//		final FactorizationOfLongs factorizer2 = new TrialDoubleFact(1 << (bits/2));
		//		final FactorizationOfLongs factorizer1 = new PollardRhoBrentParallel();
		//		final FactorizationOfLongs factorizer1 = new PollardRho((int) (Math.sqrt(begin)));
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2f, true);
		//		final FactorizationOfLongs factorizer1 = new TrialInvFact(1 << (bits/2));
		//		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1, false);


		while (begin < Long.MAX_VALUE / 1000)
		{
			for (long i = begin; i < begin + begin/1; i++) {
				//				final Collection<Long> factors = factorizer1.factorization(i);
				//				System.out.println(i + ": " + factorizer1.printFactorization(i));
				final SortedMultiset<BigInteger> factors = factorizer1.factor(BigInteger.valueOf(i));
				System.out.println(i + ": " + factors);
				//			Collection<Integer> factors = factorizer1.findAllPrimeFactors(i);
				//				final Collection<Long> factors2 = factorizer2.factorization(i);
				//				System.out.println(i + ": " + factorizer2.printFactorization(i));
				final SortedMultiset<BigInteger> factors2 = factorizer2.factor(BigInteger.valueOf(i));
				System.out.println(i + ": " + factors2);

				if (factors.size()!=factors2.size())
					System.out.println();
				//				assertEquals("Test failed for " + i, factors.size(), factors2.totalCount());
				assertEquals("Test failed for " + i, factors.totalCount(), factors2.totalCount());
				//				final Multiset<Long> factorsSet = TreeMultiset.create();
				//								factorsSet.addAll(factors2);
				//
				//								for (final long factor :
				//									factors) {
				//									assertTrue("Test failed for " + i, factorsSet.contains(factor));
				//								}
			}
			begin = begin + begin;
		}

	}

	@Test
	public void testCorrectThird() {
		final int bits = 41;

		final int numPrimes = 1634;
		final int loop = 20;
		final int smallFactorBits = bits / 3;
		final long[] semiprimes = PerformanceTest.makeSemiPrimesList(bits, smallFactorBits, numPrimes);


		long begin = (1L << bits) +1;
		final FactorizationOfLongs factorizer2 = new PollardRhoBrentDouble52();
		factorizer2.setMaxFactor(1 << smallFactorBits);
		final FactorizationOfLongs factorizer1 = new TrialFact();
		//		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1, false);


		for (int i = 0; i < loop; i++) {

		}
		{
			for (final long k : semiprimes) {
				final Collection<Long> factors = factorizer1.factorization(k);
				System.out.println(k + ": " + factorizer1.printFactorization(k));
				//			Collection<Integer> factors = factorizer1.findAllPrimeFactors(i);
				final Collection<Long> factors2 = factorizer2.factorization(k);
				System.out.println(k + ": " + factorizer2.printFactorization(k));
				//				final SortedMultiset<BigInteger> factors2 = factorizer2.factor(BigInteger.valueOf(i));
				//				System.out.println(i + ": " + factors2);

				if (factors.size()!=factors2.size())
					System.out.println();
				assertEquals("Test failed for " + k, factors.size(), factors2.size());
				final Multiset<Long> factorsSet = TreeMultiset.create();
				//				factorsSet.addAll(factors2);
				//
				//				for (final long factor :
				//					factors) {
				//					assertTrue("Test failed for " + i, factorsSet.contains(factor));
				//				}
			}
			begin = begin + begin;
		}

	}

	private List<Long> factors(long p, long q) {
		return Arrays.asList(p,q);
	}

}