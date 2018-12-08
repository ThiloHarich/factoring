package factoring.math;

import static org.junit.Assert.assertEquals;

import java.math.BigInteger;

import org.junit.Test;

public class BigDoubleTest {

	@Test
	public void testMod() {

		final long n = (1l << 51) + 389271;
		final double nInv = 1.0 / n;

		final long x = (1l << 49) - 1;
		final long y = (1l << 49) - 1;

		final long[] prod = BigDouble.multiply(x, y);
		final BigInteger prodBig = BigInteger.valueOf(x).multiply(BigInteger.valueOf(y));

		final BigInteger prodD = BigInteger.valueOf(prod[0]).shiftLeft(52).add(BigInteger.valueOf(prod[1]));
		System.out.println(prodBig);
		System.out.println(prodD);
		assertEquals(prodD, prodBig);

		final long prodMod = BigDouble.mod(prod, n, nInv);

		final BigInteger prodBigMod = prodD.mod(BigInteger.valueOf(n));
		final BigInteger prodBigDiv = prodD.divide(BigInteger.valueOf(n));
		System.out.println(prodBigDiv);
		System.out.println(prodMod);
		System.out.println(prodBigMod);
		assertEquals(prodBigMod.longValue(), prodMod);

	}
	@Test
	public void testSplit() {

		final long x = (1l << 49) - 1;

		final long[] split = BigDouble.split(x);

		final BigInteger shiftLeft = BigInteger.valueOf(split[0]).shiftLeft(26);
		final BigInteger lower = BigInteger.valueOf(split[1]);
		final BigInteger prodD = shiftLeft.add(lower);
		System.out.println(x);
		System.out.println(prodD);
		assertEquals(prodD.longValue(), x);
	}

}
