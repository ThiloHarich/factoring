package factoring.math;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.TreeMultiset;

import factoring.FactorFinder;
import factoring.fermat.lehman.LehmanFactorization;
import factoring.trial.variant.TrialFact;

public class IsPrimeTest {

	@Test
	public void testBelow1000000()
	{
		final int range = 1000000;
		final int begin = 364289;
		for (int n = begin; n < begin + range ; n+= 2)
		{
			final boolean isPrime = PrimeMath.isPrime41Bit(n);
			final FactorFinder trial = new  TrialFact();
			final long factor = trial.findFactors(n, null);
			assertEquals("" + n, factor == n, isPrime);
		}
	}

	@Test
	public void testBig()
	{
		final int range = 100;
		for (long n = 4294967311l; n< 4294967311l + range ; n+= 2)
		{
			final boolean isPrime = PrimeMath.isPrime41Bit(n);
			final FactorFinder trial = new LehmanFactorization(32, 0f);
			final long factor = trial.findFactors(n, TreeMultiset.create());
			assertEquals("" + n, factor == n, isPrime);
		}
	}

}
