package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import factoring.primes.Primes;
import factoring.trial.TDiv23InverseFMA;
import factoring.trial.TDiv31Barrett;
import org.junit.Test;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.util.SortedMultiset;
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

		long begin = (1L << bits) +7;
				begin = 9l;
//				begin = 73891306919159L;
//				begin = 24389;
//				begin = 5989;
//				begin = 515193651703l;
		//		final LehmanFactorFinder factorizer1 = new LehmanFactorFinder(50, 1, false);
		//		final FactorAlgorithm factorizer2 = new SquFoF31();
		//		final FactorAlgorithm factorizer1 = new LehmanMultiplier6_5_7(true);
		//		final FactorAlgorithm factorizer2 = new Lehman_CustomKOrderTh(false);
		//		final FactorAlgorithm factorizer1 = new Hart_FastAdjustMap(false);
//		final FactorAlgorithm factorizer1 = new Hart315Primes(false);
		//		final FactorAlgorithm factorizer2 = new HartMod8(true);
//		final FactorAlgorithm factorizer2 = new Hart_FastT(false);
//				final FactorAlgorithm factorizer2 = new SmoothNumbersSieve();
		//		final FactorAlgorithm factorizer1 = new LehmanMidRange7(0,1);
		//		final FactorAlgorithm factorizer1 = new factoring.hart.Hart_TDiv_Race();
//				final FactorAlgorithm factorizer2 = new LehmanHart2();
		//		final FactorAlgorithm factorizer2 = new TrialMultiplyUnrol(1 << (bits/2));
		//		Factorizer factorizer1 = new Fermat24();
		//		Factorizer factorizer1 = new LehmanBigFact(bitsMax, 1);
		//		final Factorizer factorizer2 = new LehmanMod16Fact(bitsMax);
		//		final FactorAlgorithm factorizer2 = new de.tilman_neumann.jml.factor.lehman.Lehman_Fast(false);
		//		final FactorAlgorithm factorizer1 = new Lehman_FastOrig(false);
		//		final FactorAlgorithm factorizer2 = new LehmanMultiplier6_5_7_11(true);
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinderRange(bits, 2f, true);
		//		final FactorizationOfLongs factorizer2 = new TrialDoubleFact(1 << (bits/2));
		//		final FactorizationOfLongs factorizer2 = new TrialInvFact(1 << (bits/2));
		//		final FactorizationOfLongs factorizer2 = new TrialFloatFact(1 << (bits/2));
		//		final FactorizationOfLongs factorizer1 = new PollardRhoBrentDouble52();
		//		final FactorizationOfLongs factorizer1 = new PollardRho((int) (Math.sqrt(begin)));
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2f, true);
		final FactorAlgorithm factorizer1 = new TDiv31Barrett();
		final FactorAlgorithm factorizer2 = new TDiv23InverseFMA();
		//		final FactorAlgorithm factorizer2 = new CombinedFactorAlgorithm(1, false);

		//		while (begin < Long.MAX_VALUE / 1000)
		//		{
		for (long i = begin; ; i+=1) {
			//				final Collection<Long> factors = factorizer1.factorization(i);
			//				System.out.println(i + ": " + factorizer1.printFactorization(i));
			final SortedMultiset<BigInteger> factors = factorizer1.factor(BigInteger.valueOf(i));
			System.out.println(i + ": " + factors);
			//				final Collection<Long> factors = factorizer1.factorization(i);
			//			final Collection<Long> factors2 = factorizer2.factorization(i);
			//				System.out.println(i + ": " + factorizer2.printFactorization(i));
			final SortedMultiset<BigInteger> factors2 = factorizer2.factor(BigInteger.valueOf(i));
			//																				final SortedMultiset<BigInteger> factors2 = factorizer2.findAllPrimeFactors(i);
			//				System.out.println(i + ": " + factorizer1.printFactorization(i));
			System.out.println(i + ": " + factors2);
			//
							if (factors.totalCount()!=factors2.totalCount())
								System.out.println();
			//			assertEquals("Test failed for " + i, factors.totalCount(), factors2.size());
//			assertEquals("Test failed for " + i, factors.totalCount(), factors2.totalCount());
			//				final Multiset<Long> factorsSet = TreeMultiset.create();
			//								factorsSet.addAll(factors2);
			//
			//								for (final long factor :
			//									factors) {
			//									assertTrue("Test failed for " + i, factorsSet.contains(factor));
			//								}
		}
		//			begin = begin + begin;
		//		}

	}

	@Test
	public void testCorrectThird() {
		final int bits = 41;

		final int numPrimes = 1634;
		final int loop = 20;
		final int smallFactorBits = bits / 3;
		boolean readFromFile = false;
		final long[] semiprimes = Primes.makeSemiPrimesList(bits, numPrimes, readFromFile);


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