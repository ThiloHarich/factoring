package factoring;
import java.util.BitSet;

public class MyMath {

	private static final long MATH_SQRT_WORKS = 10354000L;
	//	private static final long MATH_SQRT_WORKS = 4503599627370496L;
	static boolean [] squaresMod;
	static int mod = 1024;
	static BitSet squares = new BitSet(mod);


	public static boolean isSquare(long n)
	{
		if (!isProbableSquare(n))
			return false;
		long sqrt = (long)Math.sqrt(n);
		return sqrt*sqrt == n;
	}

	public static long sqrt(long n) {
		long sqrt = (long)Math.sqrt(n);
		if (n < MATH_SQRT_WORKS)
		{
			return sqrt;
		}
		long sqrt2=sqrt, tmp;
		do{
			tmp = sqrt;
			sqrt = sqrt2;
			sqrt2 = tmp;

			final long nDivSqrt = (long) ((0.0 +n)/sqrt);
			sqrt2 = (sqrt + nDivSqrt)/2;
		}
		while (sqrt2 - sqrt > 0);
		return sqrt;
	}

	public static boolean isProbableSquare(long n) {
		final int nMod = (int) (n & (mod-1));
		if (!squaresMod[nMod])
			return false;
		return true;
	}
	public static boolean isSquareBitSet(long n)
	{
		final int nMod = (int) (n & (mod-1));
		if (!squares.get(nMod))
			return false;

		final long sqrt = (long)Math.sqrt(n);

		return sqrt*sqrt == n;
	}

	static {
		squaresMod = new boolean[mod];
		for (int i = 0; i < squaresMod.length; i++) {
			final int square = (i*i) % mod;
			squaresMod[square] = true;
			squares.set(square);
		}
	}

	private final static boolean isPerfectSquare(long n)
	{
		switch((int)(n & 0x3F))
		{
		case 0x00: case 0x01: case 0x04: case 0x09: case 0x10: case 0x11:
		case 0x19: case 0x21: case 0x24: case 0x29: case 0x31: case 0x39:
			long sqrt;
			//			if(n < 410881L)
			//			{
			//				//John Carmack hack, converted to Java.
			//				// See: http://www.codemaestro.com/reviews/9
			//				int i;
			//				float x2, y;
			//
			//				x2 = n * 0.5F;
			//				y  = n;
			//				i  = Float.floatToRawIntBits(y);
			//				i  = 0x5f3759df - ( i >> 1 );
			//				y  = Float.intBitsToFloat(i);
			//				y  = y * ( 1.5F - ( x2 * y * y ) );
			//
			//				sqrt = (long)(1.0F/y);
			//			}
			//			else
			{
				//Carmack hack gives incorrect answer for n >= 410881.
				sqrt = (long)Math.sqrt(n);
			}
			return sqrt*sqrt == n;

		default:
			return false;
		}
	}

	static long goodMask; // 0xC840C04048404040 computed below
	static {
		for (int i=0; i<64; ++i) goodMask |= Long.MIN_VALUE >>> (i*i);
	}

	public static boolean isSquareM(long x) {
		// This tests if the 6 least significant bits are right.
		// Moving the to be tested bit to the highest position saves us masking.
		if (goodMask << x >= 0) return false;
		final int numberOfTrailingZeros = Long.numberOfTrailingZeros(x);
		// Each square ends with an even number of zeros.
		if ((numberOfTrailingZeros & 1) != 0) return false;
		x >>= numberOfTrailingZeros;
		// Now x is either 0 or odd.
		// In binary each odd square ends with 001.
		// Postpone the sign test until now; handle zero in the branch.
		if ((x&7) != 1 | x <= 0) return x == 0;
		// Do it in the classical way.
		// The correctness is not trivial as the conversion from long to double is lossy!
		final long tst = (long) Math.sqrt(x);
		return tst * tst == x;
	}

	public static void main(String[] args) {
		long split1, split2;
		int count2 = 0;
		split1 = System.currentTimeMillis();
		final long maxValue = 2* MATH_SQRT_WORKS;
		for (long i = 0; i < maxValue; i+=9) {
			final boolean square = isSquare(i);
			if (square)
				count2++;
		}
		split2 = System.currentTimeMillis();
		System.out.println("time : "+ (Math.abs(0.0 + split1 - split2))/1000);
		System.out.println(count2);

		int count0 = 0;
		split1 = System.currentTimeMillis();
		for (long i = 0; i < maxValue; i+=9) {
			final boolean perfectSquare = isSquareM(i);
			if (perfectSquare)
				count0++;
		}
		split2 = System.currentTimeMillis();
		System.out.println("time : "+ (Math.abs(0.0 + split1 - split2))/1000);
		System.out.println(count0);

		int count = 0;
		split1 = System.currentTimeMillis();
		for (long i = 0; i < maxValue; i+=9) {
			final boolean perfectSquare = isPerfectSquare(i);
			if (perfectSquare)
				count++;
		}
		split2 = System.currentTimeMillis();
		System.out.println("time : "+ (Math.abs(0.0 + split1 - split2))/1000);
		System.out.println(count);


		int count3 = 0;
		split1 = System.currentTimeMillis();
		for (long i = 0; i < maxValue; i+=9) {
			final boolean square = isSquareBitSet(i);
			if (square)
				count3++;
		}
		split2 = System.currentTimeMillis();
		System.out.println("time : "+ (Math.abs(0.0 + split1 - split2))/1000);
		System.out.println(count3);
	}
}
