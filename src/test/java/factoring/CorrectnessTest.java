package factoring;

import static org.junit.Assert.assertEquals;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;

import factoring.fermat.lehman.LehmanNoSqrtFact;
import factoring.fermat.lehman.LehmanSingleFactorFinder;

public class CorrectnessTest {
	@Test
	public void testCorrect() {
		final int bits = 8;
		final int bitsMax = 32;

		//		final Factorizer factorizer1 = new TrialPrimesDynamicFact(1 << bits + 4);
		//		Factorizer factorizer1 = new Fermat24();
		//		Factorizer factorizer1 = new LehmanBigFact(bitsMax, 1);
		//		final Factorizer factorizer2 = new LehmanMod16Fact(bitsMax);
				final Factorizer factorizer2 = new LehmanNoSqrtFact(bitsMax, 1);
		final FactorizationOfLongs factorizer1 = new LehmanSingleFactorFinder(bitsMax, 0f);
//		final Factorizer factorizer1 = new LehmanYafuFact(1.5f);
		//		final Factorizer factorizer2 = new TrialInvFact(1 << bits + 4);

		//		for (int i = 65538; i < 1 << (bits + 1); i++)
		long begin = (1L << bitsMax) +1;  // = 2^4 * 3^2 * 5
		begin = 37 * 37	; // * 23
		// 29*23 * 53
		// 29*53 * 23 ->
		while (begin < Long.MAX_VALUE / 1000)
		{
			for (long i = begin; i < begin + begin/8; i++) {
				final TreeMultiset<Long> factors = factorizer1.factorizationByPrimes(i);
				System.out.println(i + ": " + getPrintString(factors));

//				final Collection<Long> factors = factorizer1.findAllPrimeFactors(i);
//				System.out.println(i + ": " + factorizer1.printFactors(i));
				Collection<Long> factors2 = factorizer2.findAllPrimeFactors(i);
                System.out.println(i + ": " + factorizer2.printFactors(i));
//				SortedMultiset<BigInteger> factors2 = factorizer2.factor(BigInteger.valueOf(i));
//				final Set<Map.Entry<BigInteger, Integer>> factors2Set = factors2.entrySet();
//				System.out.println(i + ": " + getPrintString(factors2Set));

                assertEquals("Test failed for " + i, factors.size(), factors2.size());
//                assertEquals("Test failed for " + i, factors.size(), factors2.totalCount());
				final Multiset<Long> factorsSet = TreeMultiset.create();
//				factorsSet.addAll(factors2);

//				for (final long factor :
//					factors) {
//					assertTrue("Test failed for " + i, factorsSet.contains(factor));
//				}
			}
			begin = begin + begin;
		}

	}

	private String getPrintString(Set<Map.Entry<BigInteger, Integer>> factors2) {
		final List<String> s = new ArrayList<String>();
		for (final Map.Entry<BigInteger, Integer> entry : factors2) {
			final int exponent = entry.getValue();
			String part = "" + entry.getKey();
			part += exponent == 1 ? "" : "^" + exponent;
			s.add(part);
		}
		return String.join(" * ", s);
	}

	private String getPrintString(TreeMultiset<Long> factors) {
		final List<String> s = new ArrayList<String>();
		for (final Multiset.Entry<Long> entry : factors.entrySet()) {
			final int exponent = entry.getCount();
			String part = "" + entry.getElement();
			part += exponent == 1 ? "" : "^" + exponent;
			s.add(part);
		}
		return String.join(" * ", s);
	}


}
