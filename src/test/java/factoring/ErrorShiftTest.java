package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.*;

import factoring.shift.ErrorShiftFact;
import org.junit.Test;

import factoring.fermat.FermatFact;

public class ErrorShiftTest {

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



	private void checkFactors(int p) {
		final int operations = 0;
		final int fermatOperations = 0;
		for (final int q : qs())
		{
			final List<Long> factors = factors(p, q);
			//			final Factorizer fact = new ErrorShiftFact(true);
			final FermatFact fact = new FermatFact();
			//			final ErrorShift2DFact fact = new ErrorShift2DFact();
			//		final ErrorYIncShiftFact fact = new ErrorYIncShiftFact();
			fact.findFactors(p*q, factors);
			final long factor = factors.get(0);
			//			assertTrue(factors.contains(factor));
			if (!factors.contains(factor))
				System.err.println("Factor not found");
			//			operations += fact.operations;
			//			fermatOperations += fact.getOperationsFermat();
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


	private List<Long> factors(long p, long q) {
		return Arrays.asList(p,q);
	}

}
