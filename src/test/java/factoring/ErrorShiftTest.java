package factoring;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.junit.Test;

public class ErrorShiftTest {

	@Test
	public void test10037_4339 ()
	{
		final int p = 10037;
		final int q = 4339;
		final Collection<Long> factors = factors(p, q);
		final ErrorShiftFact fact = new ErrorShiftFact();
		final long factor = fact.findFactor(p, q);

		assertTrue(factors.contains(factor));
	}

	List<Integer> qs ()
	{
		return Arrays.asList(20023, 20029, 24043, 24029, 30029, 45007, 60029, 60017, 60013, 90023, 90007);
	}


	private void checkFactors(int p) {
		int operations = 0;
		int fermatOperations = 0;
		for (int q : qs())
		{
			final Collection<Long> factors = factors(p, q);
			final ErrorShiftFact fact = new ErrorShiftFact();
//		final ErrorYIncShiftFact fact = new ErrorYIncShiftFact();
			final long factor = fact.findFactor(p, q);
			assertTrue(factors.contains(factor));
			operations += fact.operations;
			fermatOperations += fact.getOperationsFermat();
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

	private Collection<Long> factors(long p, long q) {
		return Arrays.asList(p,q);
	}

}
