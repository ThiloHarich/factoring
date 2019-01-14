package factoring.math;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SqrtIntTest {


	@Test
	public void testSqrt(){
		for (int i = 1; i < (1 << 16); i++) {
			SquareRoot.sqrt(i);
			assertEquals("test failed for i :" + i,(int)Math.sqrt(i), SqrtInt.sqrt(i));

		}
	}
}
