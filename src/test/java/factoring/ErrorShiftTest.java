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
		final int p = 9721;
//		final int p = 9719;
//		final int p = 11003;
//				final int p = 10037;
//				final int q = 4337;
//		final int q = 20023; //  3 exp = 99 t=1
//		final int q = 20029; //  2 exp = 2 t=2
//		final int q = 24043; //  2 exp = 1 t=3
//		final int q = 24029; //  .7 -> s=2 t=1
//		final int q = 30029; // .5 ->
//		final int q = 32029; //  2 exp = .8 t=3,5
//		final int q = 45007; // 10 exp = 0.5  -> s=0, t=56
		final int q = 60029; //  4 exp = .4  -> s=1;2
//		final int q = 60017; // 10 exp = .6  -> s=2
//		final int q = 60013; //  2 exp = 2  -> t=1
//		final int q = 90023; //  8 exp = 1  -> s=2
//		final int q = 90007; // 15 exp = 1  -> s=3
//		int q = 140729;
		final Set<Long> factors = factors(p, q);
		final ErrorShift2DFact fact = new ErrorShift2DFact();
//		final ErrorYIncShiftFact fact = new ErrorYIncShiftFact();
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
