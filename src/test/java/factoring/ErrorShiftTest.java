package factoring;

import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

public class ErrorShiftTest {

	@Test
	public void test10037_4339 ()
	{
		final int p = 10037;
		final int q = 4339;
		final Set<Long> factors = factors(p, q);
		final ErrorShiftFact fact = new ErrorShiftFact();
		final long factor = fact.findFactor(p, q);

		assertTrue(factors.contains(factor));
	}

	@Test
	public void test10037_431 ()
	{
		final int p = 10037;
		final int q = 60013;
		final Set<Long> factors = factors(p, q);
		final ErrorShiftFact fact = new ErrorShiftFact();
		final long factor = fact.findFactor(p, q);

		assertTrue(factors.contains(factor));
	}

	private Set<Long> factors(long p, long q) {
		final Set<Long> factors = new HashSet();
		factors.add(p);
		factors.add(q);
		return factors;
	}

}
