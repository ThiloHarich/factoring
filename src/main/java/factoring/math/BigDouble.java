package factoring.math;

import java.math.BigDecimal;
import java.math.BigInteger;

public class BigDouble {

	private static final double twoPow51 = 1l << 52;
	final static long maskLower  = (1l << 26) - 1;

	/**
	 * returns 	an array x with 1/n = x[0] * 2^26 + x[1]
	 * @param n
	 * @return
	 */
	@Deprecated
	public static double [] invert2Double (long n) {
		final BigInteger nBigInt = BigInteger.valueOf(n);
		final BigDecimal nBig = new BigDecimal (nBigInt);
		final int newScale = (52*3)/10;
		final BigDecimal nInvBig = BigDecimal.ONE.setScale(2*newScale).divide(nBig);
		final long twoPow52l = 1l << 52;
		final BigDecimal twoPow52 = BigDecimal.valueOf(twoPow52l);
		final BigDecimal nInvUpper = nInvBig.multiply(twoPow52);
		final double [] ns = new double [2];
		ns[0] = nInvUpper.longValue();
		ns[1] = nInvBig.subtract(nInvUpper).multiply(twoPow52).multiply(twoPow52).longValue();

		return ns;
	}

	//	/**
	//	 * calculates x*d and returns an array representing the value.
	//	 *
	//	 * x[0]*d[0] | x[1]*d[0]
	//	 *           | x[0]*d[1] | x[1]*d[1]
	//	 *   52 bit     53 bit      52 bit
	//	 *
	//	 * @param x - represent a value x[0] * 2^26 + x[1]
	//	 * @param d - represent a value d[0] * 2^26 + d[1]
	//	 * @return x*d = xd[0] * 2^52 + xd[1] * 2^26 + xd[2]
	//	 */
	//	public static double[] multiply2(long[] x, double[] d) {
	//
	//		final double[] xd = new double [3];
	//		xd[0] = x[0] * d[0];
	//		xd[1] = x[1] * d[0] + x[0]*d[1];
	//		xd[3] = x[1] * d[1];
	//
	//		return xd;
	//	}
	//
	//	/**
	//	 * calculates x*d and returns an array representing the value.
	//	 *
	//	 * @param x - a Long value
	//	 * @param d - represent a value d[0] * 2^26 + d[1]
	//	 * @return x*d = xd[0] * 2^52 + xd[1] * 2^26 + xd[2]
	//	 */
	//	public static double[] multiply2(long x, double[] d) {
	//		final long[] x26Bit = split(x);
	//
	//		return multiply2(x26Bit, d);
	//	}

	@Deprecated
	public static long[] multiply (long[] x, long yl) {

		final long [] y = split(yl);
		final long xH =     x[0] * y[0];
		final long xM = 2 * x[0] * y[1];
		final long xL =     x[1] * y[1];
		final long xL2 = xM >> 26;
		final long xH2 = xM & maskLower;
		final long[]result = {xH+xH2, xL + xL2};
		return result;
	}

	public static long[] split(long x) {
		final long[] xD = new long [2];

		xD[0] = (x >> 26);
		xD[1] = x & maskLower;
		return xD;
	}

	/**
	 * calculates x*y and returns an array representing the value.
	 * Needs 4 multiplications
	 *
	 *            51           25          0
	 * x[0]*d[0] |
	 *      |       x[1]*d[0] |
	 *      |       x[0]*d[1] |
	 *           |              x[1]*d[1] |
	 *       77
	 *
	 * @param x - a long value lower 2^52
	 * @param y - a long value lower 2^52
	 * @return an array xy with x*y = xy[0] * 2^52 + xy[1] * 2^26 + xy[2]
	 */
	public static long[] multiply(long x, long y) {
		final long[] x26 = split(x);
		final long[] y26 = split(y);

		final long xy0 = x26[0] * y26[0];
		final long xy1 = x26[1] * y26[0] + x26[0] * y26[1];
		final long xy2 = x26[1] * y26[1];

		final long[] xy = new long [2];
		xy[0] = xy0 + (xy1 >> 26);
		xy[1] =      ((xy1 & maskLower) << 26) + xy2;
		return xy;
	}

	/**
	 * calculates x*y and returns an array representing the value.
	 * Needs 4 multiplications
	 *
	 *            51           25          0
	 * x[0]*d[0] |
	 *      |       x[1]*d[0] |
	 *      |       x[0]*d[1] |
	 *           |              x[1]*d[1] |
	 *       77
	 *
	 * @param x - a long value lower 2^52
	 * @param y - a long value lower 2^52
	 * @return an array xy with x*y = xy[0] * 2^52 + xy[1] * 2^26 + xy[2]
	 */
	public static long[] square(long x) {
		final long[] x26 = split(x);

		final long xy0 = x26[0] * x26[0];
		final long xy1 = x26[1] * x26[0] << 1;
		final long xy2 = x26[1] * x26[1];

		final long[] xy = new long [2];
		xy[0] = xy0 + (xy1 >> 26);
		xy[1] =      ((xy1 & maskLower) << 26) + xy2;
		return xy;
	}

	/**
	 * calculates x*d and returns an array representing the value.
	 *
	 *
	 * @param x - representation of a 104 Bit long num as an long array of length 4, each holding 26 bits
	 * @param y - representation of a 52  Bit double num as an double array of length 2, each holding 26 bits
	 * @return an array xy with x*y = xy[0] * 2^52 + xy[1] * 2^26 + xy[2]
	 */
	@Deprecated
	public static long[] multiply2(long[] x, double [] y) {

		// 8 multiplications
		final long xy0 = (long) (x[0] * y[0]);
		final long xy1 = (long) (x[1] * y[0] + x[0]*y[1]);
		final long xy2 = (long) (x[2] * y[0] + x[1]*y[1]);
		final long xy3 = (long) (x[3] * y[0] + x[2]*y[1]);
		final long xy4 = (long) (x[3] * y[1]);

		final long[] xy = new long [5];
		xy[0] = xy0 >> 26;
		xy[1] = xy1 >> 26 + xy0 & maskLower;
		xy[2] = xy2 >> 26 + xy1 & maskLower;
		xy[3] = xy3 >> 26 + xy2 & maskLower;
		xy[4] = xy4 >> 26 + xy3 & maskLower;
		xy[4] =             xy4 & maskLower;

		return xy;
	}

	/**
	 * Calculates x % mod;
	 * It needs no division, but 7 multiplications. With single (double) precision
	 * we need 2.
	 * multiplication
	 * @param x an array of long representing x[0] * 2^52 + x[1] this represents a number of 104 bits.
	 * @param mod the modulus
	 * @param modInverse the inverse of mod = 1.0/mod
	 * @return
	 */
	public static long mod(long[] x, long mod, double modInverse) {

		final double xDivN0 = x[0] * modInverse;
		final double xDivN1 = x[1] * modInverse;

		final long xDivN = (long) (xDivN0 * (1l << 52) + xDivN1);
		//		final long xDivN = (long) ((x[0] * (1l << 52) + x[1]) * modInverse);
		final long [] xDiff = multiply(xDivN, mod);

		// TODO theoretically xL must be 0 but is -1
		final long xL = x[0] - xDiff[0];
		final long xH = x[1] - xDiff[1];

		return (xL << 52) +  xH;
	}

	public static long mod2(long[] x, long mod, double modInverse) {
		final long xDivN = (long) ((x[0] * twoPow51 + x[1]) * modInverse);
		final long [] xDiff = multiply(xDivN, mod);

		// TODO theoretically xL must be 0 but is -1
		final long xL = x[0] - xDiff[0];
		final long xH = x[1] - xDiff[1];

		return (xL << 52) +  xH;
	}


}
