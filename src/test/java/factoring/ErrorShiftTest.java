package factoring;

import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

public class ErrorShiftTest {

	@Test
	public void test10037_4339 ()
	{
		final Set<Long> factors = new HashSet();
		factors.add(10037l);
		factors.add(4339l);
		final ErrorShiftFact fact = new ErrorShiftFact();

		final long factor = fact.findFactor(10037, 4339);

		assertTrue(factors.contains(factor));

	}

}
