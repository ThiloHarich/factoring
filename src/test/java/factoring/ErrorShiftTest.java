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
		final int p = 11003;
		//		final int p = 10037;
		//		final int q = 4337;
		final int q = 20023; // 465 -> s=4 a=1
		//		final int q = 20029; // 11 -> s=1 a=.5,
		//		final int q = 24043; // 1500 -> s=67
		//		final int q = 24029; // 7 -> s=2 a=.5
		//		final int q = 30029; // 20 -> s=5 a=.5
		//		final int q = 32029; // 10 -> s=3
		//		final int q = 45007; // 50 -> s=4, a=1
		//		final int q = 60029; // 6  -> s=1;2
		//		final int q = 60017; // 6  -> s=1
		//		final int q = 60013; // 6  -> s=1
		//		final int q = 90023; // 3  -> s=2,5
		//		final int q = 90007; // 3  -> s=2,5
		final Set<Long> factors = factors(p, q);
		final ErrorShiftFact fact = new ErrorShiftFact();
		final long factor = fact.findFactor(p, q);

		assertTrue(factors.contains(factor));
	}

	private Set<Long> factors(long p, long q) {
		final Set<Long> factors = new HashSet<Long>();
		factors.add(p);
		factors.add(q);
		return factors;
	}

}
