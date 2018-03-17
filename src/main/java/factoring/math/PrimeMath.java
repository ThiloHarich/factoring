package factoring.math;
import java.util.BitSet;

public class PrimeMath {

	private static final long MATH_SQRT_WORKS = 10354000L;
	//	private static final long MATH_SQRT_WORKS = 4503599627370496L;
	public static boolean [] squaresMod2Pow;
	//	static boolean [] squaresModSmallPrimes;
	//	static boolean [] squaresModSmallPrimes2;
	static long squaresMask;
	//	static int mod2Pow = 8;
	static int mod2Pow = 1024;
	static int modSmallPrimes = 3*3 * 5 * 7 * 11; // 3465
	static int modSmallPrimes2 = 13 * 17 * 19; // 4212
	static BitSet squares = new BitSet(mod2Pow);

	static {
		squaresMod2Pow = new boolean[mod2Pow];
		//		squaresModSmallPrimes = new boolean[modSmallPrimes];
		//		squaresModSmallPrimes2 = new boolean[modSmallPrimes2];
		//		for (int i = 0; i < squaresModSmallPrimes2.length; i++) {
		for (int i = 0; i < squaresMod2Pow.length; i++) {
			//			if (i < squaresMod2Pow.length)
			squaresMod2Pow[(i*i) % mod2Pow] = true;
			//			if (i < squaresModSmallPrimes.length)
			//				squaresModSmallPrimes[(i*i) % modSmallPrimes] = true;
			//			squaresModSmallPrimes2[(i*i) % modSmallPrimes2] = true;
			//			squares.set(square);
			//			squaresMask |= 1l << square;
		}
	}


	public static int log2(long value) {
		return Long.SIZE-Long.numberOfLeadingZeros(value);
	}

	public static int mod (long x, int mod)
	{
		int xMod = (int) (x % mod);
		if (xMod < 0)
			xMod += mod;
		return xMod;
	}

	public static boolean isSquareLong(long n)
	{
		final int nMod = (int) (n & (mod2Pow -1));
		if ((squaresMask & 1l << nMod) == 0)
			return false;
		final long sqrt = (long)Math.sqrt(n);
		return sqrt*sqrt == n;
	}

	/**
	 * This method check if the value n is a square y^2.
	 * All squares modulu mod2Pow are precomputed and stored.
	 * Thus it can be checked easily if this number might be a prime
	 * modulu mod2Pow, by looking up the square. Since mod2Pow is a power of 2
	 * the modulo (%) operation is fast.
	 *
	 * @param n
	 * @return
	 */
	public static boolean isSquare(long n)
	{
		//		if (n < 0)
		//			return false;
		if (!isProbableSquare(n))
			return false;
		final long sqrt = sqrt(n);
		return sqrt*sqrt == n;
	}

	/**
	 * calulates a^-1 mod2Pow modulus via a variant of the eucliedean algorithm.
	 * we have to determine the quotien q = dividend / divisor;
	 * One variant is to subtract divisor from the divident as long as this value
	 * is grater then zero. So we have q subtractions if q = dividend / divisor.
	 * We use a more balanced tree for the q. It looks like this:
	 *      .
	 *    /   \
	 *   1      .
	 *       /     \
	 *     .         .
	 *    / \       / \
	 *   2   .     .  >8
	 *      / \   / \
	 *     3   4 5   .
	 *              / \
	 *             6   .
	 *                / \
	 *               7  8
	 * @param a
	 * @param modulus
	 * @return
	 */
	public static long invert (long a, long modulus)
	{
		long ps1, ps2, parity, rem3, rem2, rem, q;

		if (a == 1)
			return 1;

		q        = modulus / a;
		rem      = modulus - a * q;
		rem3 = a;
		rem2  = rem;
		ps1      = q;
		ps2      = 1;
		parity   = ~0;

		while (rem2 > 1)
		{
			rem = rem3 - rem2;

			if (rem >= rem2)
			{
				rem -= rem2;         // q=1

				if (rem < rem2)
				{
					q   += ps1;
				}
				else
				{
					final long divisor4 = rem2 << 2;
					if (rem < divisor4)
					{
						rem -= rem2;        // q=2

						final long p2 = ps1<<1;
						if (rem < rem2)
						{
							q   += p2;
						} else
						{
							rem -= rem2;        // q=3

							if (rem < rem2)
							{
								q   += p2 + ps1;
							}else
							{
								q   += ps1<<2;     // q=4
								rem -= rem2;
							}
						}
					}else
					{
						final long divisor8 = rem2<<3;

						if (rem < divisor8)
						{
							rem -= divisor4;        // q=5

							if (rem < rem2)
							{
								q   += ps1 * 5;
							}
							else
							{
								rem -= rem2;        // q=6

								if (rem < rem2)
								{
									q   += ps1 * 6;
								} else
								{
									rem -= rem2;        // q=7

									if (rem < rem2)
									{
										q   += ps1 * 7;
									}else
									{
										q   += ps1<<3;     // q=8
										rem -= rem2;
									}
								}
							}
						}else
						{
							q = rem3 / rem2;        // q >8
							rem = rem3 - q * rem2;
							q *= ps1;
						}
					}
				}
			}
			rem3 = rem2;
			rem2 = rem;

			q += ps2;
			parity = ~parity;
			ps2 = ps1;
			ps1 = q;
		}

		if (parity == 0)
			return (ps1);
		else
			return (modulus - ps1);
	}

	/**
	 * returns the lowest integer s where s*s >= n
	 * @param n
	 * @return
	 */
	public static long sqrt(long n) {
		long sqrt = (long) Math.ceil(Math.sqrt(n));
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
		if (sqrt*sqrt < n)
			return sqrt + 1;
		return sqrt;
	}

	public static boolean isProbableSquare(long n) {
		final int nMod = (int) (n & 1023);
		return squaresMod2Pow[nMod];
		//		if (!squaresMod2Pow[nMod])
		//			return false;
		// % (division) is still much faster then calculating the square root
		//		nMod = (int) (n % modSmallPrimes);
		//		if (!squaresModSmallPrimes[nMod])
		//			return false;
		//		return true;
		//		nMod = (int) (n % modSmallPrimes2);
		//		return squaresModSmallPrimes2[nMod];

	}
	public static boolean isSquareBitSet(long n)
	{
		final int nMod = (int) (n & (mod2Pow -1));
		if (!squares.get(nMod))
			return false;

		final long sqrt = (long)Math.sqrt(n);

		return sqrt*sqrt == n;
	}
	/**
	 * calulates a^-1 mod modulus via a variant of the eucliedean algorithm.
	 * we have to determine the quotien q = dividend / divisor;
	 * One variant is to subtract divisor from the divident as long as this value
	 * is grater then zero. So we have q subtractions if q = dividend / divisor.
	 * We use a more balanced tree for the q. It looks like this:
	 *      .
	 *    /   \
	 *   1      .
	 *       /     \
	 *     .         .
	 *    / \       / \
	 *   2   .     .  >8
	 *      / \   / \
	 *     3   4 5   .
	 *              / \
	 *             6   .
	 *                / \
	 *               7  8
	 * @param a
	 * @param modulus
	 * @return
	 */
	public static long gcdKruppaTree10 (long a, long modulus)
	{
		long rem3, rem2, rem, q;

		if (a == 1)
			return 1;

		q        = modulus / a;
		rem      = modulus - a * q;
		rem3 = a;
		rem2  = rem;

		while (rem2 > 1)
		{
			rem = rem3 - rem2;

			if (rem >= rem2)
			{
				rem -= rem2;         // q=1

				if (rem < rem2)
				{
//					q   += ps1;
				}
				else
				{
					long divisor4 = rem2 << 2;
					if (rem < divisor4)
					{
						rem -= rem2;        // q=2

//						long p2 = ps1<<1;
						if (rem < rem2)
						{
//							q   += p2;
						} else
						{
							rem -= rem2;        // q=3

							if (rem < rem2)
							{
//								q   += p2 + ps1;
							}else
							{
//								q   += ps1<<2;     // q=4
								rem -= rem2;
							}
						}
					}else
					{
						long divisor8 = rem2<<3;

						if (rem < divisor8)
						{
							rem -= divisor4;        // q=5

							if (rem < rem2)
							{
//								q   += ps1 * 5;
							}
							else
							{
								rem -= rem2;        // q=6

								if (rem < rem2)
								{
//									q   += ps1 * 6;
								} else
								{
									rem -= rem2;        // q=7

									if (rem < rem2)
									{
//										q   += ps1 * 7;
									}else
									{
//										q   += ps1<<3;     // q=8
										rem -= rem2;
									}
								}
							}
						}else
						{
							q = rem3 / rem2;        // q >8
							rem = rem3 - q * rem2;
//							q *= ps1;
						}
					}
				}
			}
			rem3 = rem2;
			rem2 = rem;

//			q += ps2;
//			parity = ~parity;
//			ps2 = ps1;
//			ps1 = q;
		}
		return rem3;
	}

	public static long gcd10 (long a, long b)
	{
		long rem, rem3, rem2;

		if (b > a)
		{
			final long tmp = a;
			a = b;
			b = tmp;
		}
		if (b == 0)
			return a;

		rem3 = a;
		rem2 = b;
		rem  = a % b;


		while (rem > 0)
		{
			rem = rem3 - rem2;

			if (rem >= rem2)
			{
				rem -= rem2;// q=1

				if (rem >= rem2)
				{
					rem -= rem2;         // q=2

					if (rem >= rem2)
					{
						if (rem >> 2 < rem2)
						{
							rem = remLower4 (rem, rem2);	// q=3,4,5
						}
						else
						{
							rem -= rem2 << 2;        // q=6

							if (rem >> 2 < rem2)
							{

								rem = rem4 (rem, rem2);
								//								if (rem >= rem2)
								//								{
								//									rem = rem3 (rem, rem2);		// q=7,8,9
								//								}
							}
							else
							{
								//								rem += rem2 << 2;
								rem = rem3 % rem2 - 4;
							}
						}
					}
				}
			}
			rem3 = rem2;
			rem2 = rem;
		}
		return rem3;
	}
	protected static long rem4 (long rem, long rem2)
	{
		//		rem -= rem2;        // q=0

		if (rem >= rem2)
		{
			rem = remLower4 (rem, rem2);
		}
		return rem;
	}
	protected static long rem2 (long rem, long rem2)
	{
		rem -= rem2;        // q=0

		if (rem >= rem2)
		{
			rem -= rem2;    // q=1
		}
		return rem;
	}
	/**
	 * @param rem
	 * @param rem2
	 * @return
	 */
	protected static long remLower4 (long rem, long rem2)
	{
		rem -= rem2;        // q=0

		if (rem >= rem2)
		{
			rem = rem2 (rem, rem2);
		}
		return rem;
	}

	public static long gcd (long a, long b)
	{
		a = Math.abs (a);

		long rem, divident, divisor;

		if (b > a)
		{
			final long tmp = a;
			a = b;
			b = tmp;
		}
		if (b == 0)
			return a;

		divident = a;
		divisor = b;
		rem = a % b;

		if (rem == 0)
			return divisor;

		while (rem > 0)
		{
			rem = divident - divisor;

			if (rem >= divisor)
			{
				rem -= divisor; // rem/rem2 = 1

				if (rem >= divisor)
				{
					rem -= divisor; // rem/rem2 = 2

					if (rem >= divisor)
					{
						rem = divident % divisor;
					}
				}
			}

			divident = divisor;
			divisor = rem;
		}
		return divident;
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
		// Now xArray is either 0 or odd.
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
			//			final boolean square = isSquareBitSet(i);
			//			if (square)
			count3++;
		}
		split2 = System.currentTimeMillis();
		System.out.println("time : "+ (Math.abs(0.0 + split1 - split2))/1000);
		System.out.println(count3);
	}
}
