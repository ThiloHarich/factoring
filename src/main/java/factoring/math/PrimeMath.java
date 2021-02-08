package factoring.math;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.BitSet;

public class PrimeMath {

	private static final long MATH_SQRT_WORKS = 10354000L;
	private static final int GCD_MASK = 0xFFFF; //65535; //(2 << 16) -1;
	private static final int GCD_MASK_LENGHT = 16; //65535; //(2 << 16) -1;
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
	private static long[] gcdTable;

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
	/**
	 *
	 * @param n has to be odd
	 * @param a, 1 < a < n-1
	 * @return false if the number is composite. true if a^(n-1) == 1 mod n.
	 * According to fermats little theorem, all prime numbers (but also other numbers)
	 * are fulfilling this condition.
	 */
	public static boolean isProbablePrime (int n, int a) {
		// write n-1 = d*2^j, j>0 since n-1 is odd
		final long nMinus1 = n - 1;
		long d = nMinus1 >> 1;
		int j = 1;
		for(; (d & 1) == 0; j++, d >>= 1);
		long t = 1;
		long p = a;

		//square and multiply: t = a^d mod n
		while (d > 0) {
			if ((d & 1) != 0)
				t = t*p % n;
			p = p*p % n;
			d >>= 1;
		}
		// if a^d == +/- 1 mod n -> a^(d*j) == 1 -> n is probable prime
		if (t == 1 || t == nMinus1)
			return true;

		// if a^(d*k) == -1 mod n -> a^(d*j) == 1 -> n is probable prime
		for (int k=1 ; k<j && (t != nMinus1); k++) {
			t = t*t % n;
		}
		// if we found a^(d*k) == -1 -> n is probable prime
		return t == nMinus1;
	}

	public static boolean isProbablePrime (BigInteger n, long a) {
		final BigInteger nMinusOne = n.subtract(BigInteger.ONE);
		final int j = nMinusOne.getLowestSetBit();
		final BigInteger d = nMinusOne.shiftRight(j);
		final BigInteger two = BigInteger.valueOf(2);

		//		t = a^d mod n
		BigInteger t = BigInteger.valueOf(a).modPow(d, n);

		if (t.equals(BigInteger.ONE) || t.equals(nMinusOne))
			return true;

		for (int k=1 ; k<j && !t.equals(nMinusOne); k++) {
			t = t.modPow(two, n);
		}
		return t.equals(nMinusOne);
	}

	/**
	 * deterministic check if a long number is a prime.
	 * Is based on https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
	 * @param n
	 * @return
	 */
	public static boolean isPrime(long n) {
		if (n < Integer.MAX_VALUE) {
			final long[] primes61 = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
			return Arrays.binarySearch(primes61, n) > 0 || isPrime((int) n, new int [] {2, 7, 61});
		}
		if (n < 2152302898747l) {
			final long [] as = {2, 3, 5, 7, 11};
			return Arrays.binarySearch(as, n) > 0 || isPrime(n, as);
		}
		if (n < 341550071728321l){
			final long [] as = {2, 3, 5, 7, 11, 13, 17};
			return Arrays.binarySearch(as, n) > 0 || isPrime(n, as);
		}
		final long [] as = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
		return Arrays.binarySearch(as, n) > 0 || isPrime(n, as);
	}

	private static boolean isPrime(int n, final int[] as) {
		final long nMin1 = n - 1;
		for (int i = 0; i < as.length && as[i] < nMin1; i++) {
			final boolean probablePrime = isProbablePrime(n, as[i]);
			if (!probablePrime)
				return false;
		}
		return true;
	}

	public static boolean isPrime (long n, long [] as) {
		final long nMin1 = n - 1;
		for (int i = 0; i < as.length && as[i] < nMin1; i++) {
			final boolean probablePrime = isProbablePrime(BigInteger.valueOf(n), as[i]);
			if (!probablePrime)
				return false;
		}
		return true;
	}

	/**
	 * This a deterministic variant of the Rabin Miller test for up to 32 Bits.
	 * It is sufficient to call the Rabin miller test with a = 2, 7, 61
	 * @param n
	 * @return
	 */
	public static boolean isPrime32 (int n) {
		if (n <= 3 || n %2 == 0)
			return n <= 3;
		final int [] as = {2, 7, 61};
		return isPrime(n, as);
	}


	/**
	 * This a deterministic variant of the Rabin Miller test for up to 41 Bits.
	 * It is sufficient to call the rabin miller test with a = 2, 3, 5, 7, 11
	 * @param n
	 * @return
	 */
	public static boolean isPrime41Bit (long n) {
		if (n <= 3 || n %2 == 0)
			return n <= 3;
		final int [] as = {2, 3, 5, 7, 11};
		final long nMin1 = n - 1;
		for (int i = 0; i < as.length && as[i] < nMin1; i++) {
			final boolean probablePrime = isProbablePrime(BigInteger.valueOf(n), as[i]);
			//			final boolean probablePrime = isProbablePrime((int) n, as[i]);
			if (!probablePrime)
				return false;
		}
		return true;
	}
	//
	//	/**
	//	 * This a deterministic variant of the Rabin Miller test for up to 61 Bits.
	//	 * It is sufficient to call the rabin miller test with a = 2, 3, 5, 7, 11, 13, 17, 19, 23
	//	 * @param n
	//	 * @return
	//	 */
	//	public static boolean isPrime61Bit (long n) {
	//		final int [] as = {2, 3, 5, 7, 11, 13, 17, 19, 23};
	//		if (n <= 26) {
	//			for(final int a : as)
	//				if (n == a)
	//					return true;
	//			return false;
	//		}
	//		final long nMin1 = n - 1;
	//		for (int i = 0; i < as.length && as[i] < nMin1; i++) {
	//			final boolean probablePrime = isProbablePrime(BigInteger.valueOf(n), as[i]);
	//			if (!probablePrime)
	//				return false;
	//		}
	//		return true;
	//	}


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

	public static int mod (long x, long mod)
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
	 * Binary gcd implementation. From Till Neumann
	 * @param m
	 * @param n
	 * @return gcd(m, n)
	 *
	 */
	public static long gcdBinaryTill/*_binary*/(long m, long n) {
		m = Math.abs(m);
		n = Math.abs(n);
		final long mCmp1 = m-1;
		final long nCmp1 = n-1;

		if (mCmp1>0 && nCmp1>0)
		{
			final int m_lsb = Long.numberOfTrailingZeros(m);
			final int n_lsb = Long.numberOfTrailingZeros(n);
			final int shifts = Math.min(m_lsb, n_lsb);
			m = m>>m_lsb;
										n = n>>n_lsb;
								// now m and n are odd
								//LOG.debug("m=" + m + ", n=" + n + ", g=" + g);
								while (m > 0) {
									final long t = (m-n)>>1;
						if (t<0) {
							n = -t;
							n = n>>Long.numberOfTrailingZeros(n);
						} else {
							m = t>>Long.numberOfTrailingZeros(t);
						}
						//LOG.debug("m=" + m + ", n=" + n);
								}
								final long gcd = n<<shifts;
								//LOG.debug("gcd=" + gcd);
								return gcd;
		}
		if (mCmp1<0) return n;
		if (nCmp1<0) return m;
		// else one argument is 1
		return 1;
	}

	/**
	 * Stolen from java.math.MutableBigInteger#binaryGcd(int, int)
	 * @param a
	 * @param b
	 * @return
	 */
	static int gcdBinary(int a, int b) {
		if (b == 0)
			return a;
		if (a == 0)
			return b;

		// Right shift a & b till their last bits equal to 1.
		final int aZeros = Integer.numberOfTrailingZeros(a);
		final int bZeros = Integer.numberOfTrailingZeros(b);
		a >>>= aZeros;
		b >>>= bZeros;

						final int t = (aZeros < bZeros ? aZeros : bZeros);

						while (a != b) {
							if ((a+0x80000000) > (b+0x80000000)) {  // a > b as unsigned
								a -= b;
								a >>>= Integer.numberOfTrailingZeros(a);
							} else {
								b -= a;
								b >>>= Integer.numberOfTrailingZeros(b);
							}
						}
						return a<<t;
	}


	/**
	 * use the fastest implementation
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcd(long a, long b) {
		final Long e;
		return gcdBinary(a, b);
	}

	/**
	 * use the fastest implementation
	 *
	 * if we do this often we can reuse the gcd from previous checks.
	 * take a hash from a,b          (int)(a ^ b ^ (a >>> 32) ^(a >>> 32)) & MASK
	 * store all results for samml a,b -> shift b by MASK / 2
	 *
	 * a1 a2
	 * a3 a4
	 *    b1
	 * b2 b3
	 * b4
	 *
	 * if only a1, b1 are set they do not interfere
	 * Look it up.
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdCached(long a, long b) {
		if (gcdTable == null)
			gcdTable = new long [4*(GCD_MASK+1)];

		final int hash = (int)(a ^ b ^ (a >>> 32) ^(a >>> 32)) & GCD_MASK;

		if (gcdTable[4*hash] != 0 && gcdTable[4*hash] == a && gcdTable[4*hash+1] == b)
			return gcdTable[4*hash+2];

		final long gcd = gcdBinary(a, b);

		gcdTable[4*hash] = a;
		gcdTable[4*hash+1] = b;
		gcdTable[4*hash+2] = gcd;
		return gcd;
	}

	/**
	 * Stolen from java.math.MutableBigInteger#binaryGcd(int, int)
	 * and adapted to long
	 * @param a
	 * @param b
	 * @return
	 */
	static long gcdBinary(long a, long b) {
		if (b == 0)
			return a;
		if (a == 0)
			return b;

		// Right shift a & b till their last bits equal to 1.
		final int aZeros = Long.numberOfTrailingZeros(a);
		final int bZeros = Long.numberOfTrailingZeros(b);
		a >>>= aZeros;
		b >>>= bZeros;

		final int t = (aZeros < bZeros ? aZeros : bZeros);

		while (a != b) {
			if ((a+0x8000000000000000l) > (b+0x8000000000000000l)) {  // a > b as unsigned
				a -= b;
				a >>>= Long.numberOfTrailingZeros(a);
			} else {
				b -= a;
				b >>>= Long.numberOfTrailingZeros(b);
			}
		}
		return a<<t;
	}



	/**
	 * Old Euclidean algorithm based on substrations from school book.
	 * Slow -> do not use it. Just to compare other variants against.
	 * @param a
	 * @param b
	 * @return
	 */
	@Deprecated
	public static long gcdBySubstractionOnly (long a, long b) {
		if (a == 0)
			return b;
		while (b > 0)
			if (a>b)
				a = a - b;
			else
				b = b - a;
		return a;
	}

	/**
	 * New Euclidean algorithm based on modulus from school book.
	 * Slow ->  do not use it. Just to compare other variants against.
	 * @param a
	 * @param b
	 * @return
	 */
	@Deprecated
	public static long gcdByMod (long a, long b) {
		while (b > 0) {
			final long tmp = a % b;
			a = b;
			b = tmp;
		}
		return a;
	}

	/**
	 * New Euclidean algorithm based on modulus from school book.
	 * Slow ->  do not use it. Just to compare other variants against.
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdByMod (double a, double b) {
		while (b > 0) {
			final double tmp = a % b;
			a = b;
			b = tmp;
		}
		return (long) a;
	}
	/**
	 * A variant of the Euclidean algorithm.
	 * We try to not use the expensive modulus operation "%".
	 * Instead we try to use
	 * 1) subtract the lower number from the bigger one and if the bigger number is still bigger then the lower
	 * 2) we reduce the bigger number by shifted values of the lower number
	 * the cost of calculating big / small by shifts is limited by
	 * by shifts is log(big / small). Since the probability that big / small > i is less then 1/i,
	 * the expected cost = Sum log(i)/i = 1/2 log(x)^2.
	 * Instead of reducing the bigger number by the modulus of the lower number in step 2,
	 *
	 * It is around 2 times faster then the Algorithms which are just based on
	 * subtractions or mod. It is ~ 10 % slower then the binary euclidean algorithm.
	 *
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdBySubstraction2 (long a, long b)
	//	public static long gcdBySubstractionAndShifts (long a, long b)
	{
		if (a == 1 || b == 1)
			return 1;

		long big = Math.max(a, b);
		long small = Math.min(a, b);

		// usually we check here for small > 0; but handling small == 1 separately saves time
		while (small > 1)
		{
			// do substractions for some time. The number of iterated substractions has a low impact on the performance
			big -= small;         // big/small =1
			if (big > small) {
				big -= small;         // big/small =2
				if (big > small) {
					big -= small;         // big/small =3
					if (big > small) {
						big -= small;         // big/small =4
						if (big > small) {
							big -= small;         // big/small =5
							if (big > small) {
								big -= small;         // big/small =6
								if (big > small) {
									big -= small;         // big/small =7
									if (big > small) {
										big -= small;         // big/small =8
										if (big > small) {
											big -= small;         // big/small =9
											if (big > small) {
												big -= small;         // big/small =10
												if (big > small) {
													big -= small;         // big/small =11
													if (big > small) {
														big -= small;         // big/small =12
														if (big > small) {
															big -= small;         // big/small =13
															if (big > small) {
																big -= small;         // big/small =14
																if (big > small) {
																	big -= small;         // big/small =15
																	if (big > small) {
																		big -= small;         // big/small =16
																		if (big > small) {
																			// we will calculate  big = big  % small in a more efficient way.
																			// since big and small are usually not completely different in size, we
																			// approximate the quotient big/small by just looking at the length of the
																			// numbers. Then we subtract the small number shifted by the (log) of the approximation of the quotient from the
																			// bigger one.
																			while (big > small) {
																				final int quotientApprox = Long.numberOfLeadingZeros(small) - Long.numberOfLeadingZeros(big);
																				final int shifts = quotientApprox < 1 ? 0 : quotientApprox -1;
																				big -= small << shifts;
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			final long tmp = small;
			small = big;
			big = tmp;
		}
		if (small == 1)
			return 1;

		return big;
	}
	/**
	 * A variant of the Euclidean algorithm.
	 * We try to not use the expensive modulus operation "%".
	 * Instead we try to use
	 * 1) subtract the lower number from the bigger one and if the bigger number is still bigger then the lower
	 * 2) we reduce the bigger number by shifted values of the lower number
	 * the cost of calculating big / small by shifts is limited by
	 * by shifts is log(big / small). Since the probability that big / small > i is less then 1/i,
	 * the expected cost = Sum log(i)/i = 1/2 log(x)^2.
	 * Instead of reducing the bigger number by the modulus of the lower number in step 2,
	 *
	 * It is around 2 times faster then the Algorithms which are just based on
	 * subtractions or mod. It is ~ 10 % slower then the binary euclidean algorithm.
	 *
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdBySubtraction (long a, long b)
	//	public static long gcdBySubstractionAndShifts (long a, long b)
	{
		if (a == 1 || b == 1)
			return 1;

		long big = Math.max(a, b);
		long small = Math.min(a, b);

		// usually we check here for small > 0; but handling small == 1 separately saves time
		while (small > 1)
		{
			// do substractions for some time. The number of iterated substractions has a low impact on the performance
			big -= small;         // big/small =1
			if (big > small) {
				big -= small;         // big/small =2
				if (big > small) {
					big -= small;         // big/small =3
					if (big > small) {
						big -= small;         // big/small =4
						if (big > small) {
							big -= small;         // big/small =5
							if (big > small) {
								big -= small;         // big/small =6
								if (big > small) {
									big -= small;         // big/small =7
									if (big > small) {
										big -= small;         // big/small =8
										if (big > small) {
											big -= small;         // big/small =9
											if (big > small) {
												big -= small;         // big/small =10
												if (big > small) {
													big -= small;         // big/small =11
													if (big > small) {
														big -= small;         // big/small =12
														if (big > small) {
															big -= small;         // big/small =13
															if (big > small) {
																big -= small;         // big/small =14
																if (big > small) {
																	big -= small;         // big/small =15
																	if (big > small) {
																		big -= small;         // big/small =16
																		if (big > small) {
																			// we will calculate  big = big  % small in a different way.
																			// since big and small are usually not completely different in size, we
																			// approximate the quotient big/small by just looking at the length of the
																			// numbers. Then we subtract the small number shifted by the (log) of the approximation of the quotient from the
																			// bigger one.
																			big = modByShiftAndSub(big, small);
																			//																			big = big % small;
																			//																																																									}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			final long tmp = small;
			small = big;
			big = tmp;
		}
		if (small == 1)
			return 1;

		return big;
	}
	/**
	 * A variant of the Euclidean algorithm.
	 * We try to not use the expensive modulus operation "%".
	 * Instead we try to use
	 * 1) subtract the lower number from the bigger one and if the bigger number is still bigger then the lower
	 * 2) we reduce the bigger number by shifted values of the lower number
	 * the cost of calculating big / small by shifts is limited by
	 * by shifts is log(big / small). Since the probability that big / small > i is less then 1/i,
	 * the expected cost = Sum log(i)/i = 1/2 log(x)^2.
	 * Instead of reducing the bigger number by the modulus of the lower number in step 2,
	 *
	 * It is around 2 times faster then the Algorithms which are just based on
	 * subtractions or mod. It is ~ 10 % slower then the binary euclidean algorithm.
	 *
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdDouble (double a, double b)
	{
		if (a == 1 || b == 1)
			return 1;

		double big = Math.max(a, b);
		double small = Math.min(a, b);
		final Double s;

		// usually we check here for small > 0; but handling small == 1 separately saves time
		while (small > 1)
		{
			// do substractions for some time. The number of iterated substractions has a low impact on the performance
			big -= small;         // big/small =1
			if (big > small) {
				big -= small;         // big/small =2
				if (big > small) {
					big -= small;         // big/small =3
					if (big > small) {
						big -= small;         // big/small =4
						if (big > small) {
							big -= small;         // big/small =5
							if (big > small) {
								big -= small;         // big/small =6
								if (big > small) {
									// we will calculate  big = big  % small in a different way.
									// since big and small are usually not completely different in size, we
									// approximate the quotient big/small by just looking at the length of the
									// numbers. Then we subtract the small number shifted by the (log) of the approximation of the quotient from the
									// bigger one.
									//																			big = modByShiftAndSub(big, small);
									big = big % small;
								}
							}
						}
					}
				}
			}
			final double tmp = small;
			small = big;
			big = tmp;
		}
		if (small == 1)
			return 1;

		return (long) big;
	}

	/**
	 * A variant of the Euclidean algorithm.
	 * We try to not use the expensive modulus operation "%".
	 * Instead we try to use
	 * - subtract the lower number from the bigger one and if the bigger number is still bigger then the lower
	 * - we reduce the bigger number by shifted values of the lower number
	 * We do not need to reduce the bigger number by the modulus of the lower number.
	 * It is around 2 times faster then the Algorithms which are just based on
	 * subtractions or mod. It is ~ 10 % slower then the binary euclidean algorithm.
	 *
	 * @param a
	 * @param b
	 * @return
	 */
	public static long gcdBySubtraction2Var (long a, long b)
	{
		if (a == 1 || b == 1)
			return 1;

		long big = Math.max(a, b);
		long small = Math.min(a, b);

		// usually we check here for small > 0; but handling small == 1 separately saves time
		while (small > 1)
		{
			// introduce a new variable, but looks like it is not reducing the hole time
			long newSmall = big - small;         // q=1
			if (newSmall > small) {
				newSmall -= small;         // q=2
				if (newSmall > small) {
					newSmall -= small;         // q=2
					if (newSmall > small) {
						newSmall -= small;         // q=2
						if (newSmall > small) {
							newSmall -= small;         // q=2
							if (newSmall > small) {
								newSmall -= small;         // q=2
								if (newSmall > small) {
									newSmall -= small;         // q=2
									if (newSmall > small) {
										newSmall -= small;         // q=2
										if (newSmall > small) {
											newSmall -= small;         // q=2
											if (newSmall > small) {
												newSmall -= small;         // q=2
												if (newSmall > small) {
													newSmall -= small;         // q=2
													if (newSmall > small) {
														newSmall -= small;         // q=2
														if (newSmall > small) {
															newSmall -= small;         // q=2
															if (newSmall > small) {
																newSmall = modByShiftAndSub(newSmall, small);
																//																while (newSmall > small) {
																//																	final int quotientApprox = Long.numberOfLeadingZeros(small) - Long.numberOfLeadingZeros(newSmall);
																//																	final int shifts = quotientApprox < 1 ? 0 : quotientApprox -1;
																//																	newSmall -= small << shifts;
																//																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			big = small;
			small = newSmall;
		}
		if (small == 1)
			return 1;

		return big;
	}


	/**
	 * calculates number % modulus by reducing number by a power of modulus.
	 *
	 * @param number
	 * @param modulus
	 * @return
	 */
	public static long modByShift2(long number, long modulus) {
		while (number > modulus) {
			final int quotientApprox = Long.numberOfLeadingZeros(modulus) - Long.numberOfLeadingZeros(number);
			final int shifts = quotientApprox < 1 ? 0 : quotientApprox -1;
			number -= modulus << shifts;
		}
		return number;
	}

	/**
	 * calculates number % modulus by reducing number by a power of modulus.
	 *
	 * @param number
	 * @param modulus
	 * @return
	 */
	public static long modByShift(long number, long modulus) {
		//		if (modulus == 0)
		//			return 0;
		while (number > modulus) {
			final int quotientApprox = Long.numberOfLeadingZeros(modulus) - Long.numberOfLeadingZeros(number);
			final int shifts = quotientApprox < 1 ? 0 : quotientApprox -1;
			number -= modulus << shifts;
		}
		if (number == modulus)
			return 0;
		return number;
	}
	/**
	 * calculates number % modulus by reducing number by a power of modulus.
	 *
	 * @param number
	 * @param modulus
	 * @return
	 */
	public static long modByShiftAndSub(long number, long modulus) {
		//		if (modulus == 0)
		//			return 0;
		while (number > modulus) {
			final int quotientApprox = Long.numberOfLeadingZeros(modulus) - Long.numberOfLeadingZeros(number);
			// just use subtractions for small numbers
			if (quotientApprox < 4) {
				while (number > modulus)
					number -= modulus;
				return number == modulus ? 0 : number;
			}

			final int shifts = quotientApprox < 1 ? 0 : quotientApprox -1;
			number -= modulus << shifts;
		}
		return number == modulus ? 0 : number;
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

	/**
	 * Calculates number % modulus, when modulusInv = 1.0 / modulus.
	 * This can be around 4 times faster then number % modulus, but will only work
	 * for number < 2^52, since double uses only 52 bits for the fraction.
	 * If the modulus is precalculated and the method is called at least twice
	 * with the same modulusInv this is overall faster then using the % on long numbers.
	 * @param number
	 * @param modulus
	 * @param modulusInv
	 * @return
	 */
	public static long mod(long number, long modulus, final double modulusInv) {
		final double DISCRIMINATOR = 1.0/(1<<10);
		final double numberDivMod = number * modulusInv;
		final long divRounded = (long)numberDivMod;
		final long remainder = number - modulus * divRounded;
		return remainder;
	}

	public static double mod(double number, double modulus, final double modulusInv) {
		final double numberDivMod = number * modulusInv;
		final double divRounded = Math.floor(numberDivMod);
		final double remainder = number - modulus * divRounded;
		return remainder;
	}

	/**
	 * Calculates number % modulus,
	 * @param number
	 * @param modulus
	 * @param modulusInv
	 * @return
	 */
	public static long mod(long number, long modulus, int multiplier, int shifts) {
		final long divRounded = (number * multiplier) >> shifts;
		long remainder = number - modulus * divRounded;
		if (remainder > modulus)
			remainder -= modulus;
		return remainder;
	}

	public static int [] barrettParams(long mod) {
		final int shifts = 64 - Long.numberOfLeadingZeros(mod -1);
		final int multiplier = (int) ((1l << shifts) / mod);
		final int [] params = {multiplier, shifts};

		return params;
	}


	//	public static long mod(int number, double i, double iInv) {
	//		// TODO Auto-generated method stub
	//		return 0;
	//	}


}
