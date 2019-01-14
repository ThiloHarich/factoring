package factoring.math;

/**
 * if x = a*2^k + b
 * sqrt(a*2^k + b) =
 * sqrt(a*2^k* (1 + b/a*2^-k)) =
 * sqrt(a*2^k)* sqrt(1 + b/a*2^-k)) =
 * sqrt(a)*2^(k/2)* sqrt(1 + b/a*2^-k)) ~
 * sqrt(a)*2^(k/2)* (1 + b/a*2^(-(k+1))
 *
 * we choose k = 2^i e.g. k=256=2^8 or k=1024
 * precalculate sqrt(a)*2^(k/2) and 1/a as doubles
 *
 * if x = a*2^(k+2i) +b
 * sqrt(a*2^(k+2i) + b) =
 * sqrt(a)*2^(k/2+i)* (1 + b/a*2^(-(k+1))
 *
 * @author thiloharich
 *
 */
public class SqrtInt {


	private static final int MAKE_EVEN = Integer.MAX_VALUE - 1;
	private static int mod = 256;    // 2^8
	private static int sqrtMod = 16; // 2^4 = sqrt(mod)
	private static int sqrtModBits = 4;
	//	private static int mod = 1024;    // 2^10
	//	private static int sqrtMod = 32; // 2^5 = sqrt(mod)
	//	private static int sqrtModBits = 5;
	static int [] sqrts = new int [mod];
	static double [] inverse = new double [mod];

	static {
		for (int i = 0; i < sqrts.length; i++) {
			sqrts[i] = (int) (Math.sqrt(i) * sqrtMod);
			inverse[i] = 1.0 / i;
			//			assertEquals(table[i], (int)sqrts[i]);
			//			System.out.println(table[i] + "\t" + sqrts[i]);
		}
	}

	public static boolean isSquare(int x) {
		final int bits = (25 - Integer.numberOfLeadingZeros(x)) & MAKE_EVEN;
		//		final int bits = 23 - Integer.numberOfLeadingZeros(x);
		// shift the number such that a are the first 8 bits
		final int shifts = x < mod ? 0 : bits;
		final int a = sqrts[x >> shifts];
		int sqrtA = (a >> (sqrtModBits - (shifts >> 1) ));
		if (sqrtA * sqrtA == x)
			return true;
		sqrtA++;
		if (sqrtA * sqrtA == x)
			return true;
		return false;
	}
	public static int sqrt(int x) {
		final int bits = 25 - Integer.numberOfLeadingZeros(x);
		//		final int bits = 23 - Integer.numberOfLeadingZeros(x);
		// shift the number such that a are the first 8 bits
		final int shifts = x < mod ? 0 : bits & MAKE_EVEN;
		final int index = x >> shifts;
		//		final int a = table[x >> shifts];
		final int a = sqrts[index];
		final int sqrtShifts = shifts >> 1;
		int sqrtA ;
		//		if (x < 65536) {
		sqrtA = (a >> (sqrtModBits - sqrtShifts )) +1;
		//		}
		//		else {
		//			// do a correction
		//			sqrtA = (a << (sqrtShifts - sqrtModBits)) +1;
		//			final int LOWER_BITS_MASK = (1 << bits) -1;
		//			final int b = x & LOWER_BITS_MASK;
		//			final double bDivA = b * inverse[index];
		//			final double pow = Math.pow(2,-(bits + 6 - shifts + 2));
		//			final double adjust = bDivA * pow;
		//			sqrtA = (int) (sqrtA * (1 + adjust));
		//	}
		//		return x < mod ? sqrtA - 1: sqrtA;
		return sqrtA * sqrtA > x ? sqrtA - 1: sqrtA;
		//		if (sqrtA * sqrtA > x)
		//			sqrtA -= 1;
		//		return sqrtA;

	}
}
