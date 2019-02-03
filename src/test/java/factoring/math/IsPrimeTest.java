package factoring.math;

import static org.junit.Assert.assertEquals;

import java.math.BigInteger;

import org.junit.Test;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import factoring.FactorFinder;
import factoring.fermat.lehman.Lehman_FastOrig;
import factoring.trial.TrialMultiply;

public class IsPrimeTest {

	@Test
	public void testBelow1000000()
	{
		final int range = 1000000;
		final int begin = 364289;
		final FactorFinder trial = new  TrialMultiply(range);
		for (int n = begin; n < begin + range ; n+= 2)
		{
			final boolean isPrime = PrimeMath.isPrime(n);
			final long factor = trial.findFactors(n, null);
			assertEquals("" + n, factor == n, isPrime);
		}
	}

	@Test
	public void testBig()
	{
		final int range = 10000;
		for (long n = 341550071728327l; n< 341550071728321l + range ; n+= 2)
		{
			final boolean isPrime = PrimeMath.isPrime(n);
			final FactorAlgorithmBase trial = new Lehman_FastOrig(true);
			final BigInteger factor = trial.findSingleFactor(BigInteger.valueOf(n));
			assertEquals("" + n, factor.longValue() == n || factor.longValue() < 2, isPrime);
		}
	}

}
