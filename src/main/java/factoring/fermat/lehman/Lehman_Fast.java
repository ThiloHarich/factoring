/*
 * java-math-library is a Java library focused on number theory, but not necessarily limited to it. It is based on the PSIQS 4.0 factoring project.
 * Copyright (C) 2018 Tilman Neumann (www.tilman-neumann.de)
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program;
 * if not, see <http://www.gnu.org/licenses/>.
 */
package factoring.fermat.lehman;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * Fast implementation of Lehman's factor algorithm.
 * Works flawlessly for N up to 60 bit.<br><br>
 *
 * It is quite surprising that the exact sqrt test of <code>test = a^2 - 4kN</code> works for N >= 45 bit.
 * At that size, both a^2 and 4kN start to overflow Long.MAX_VALUE.
 * But the error - comparing correct results vs. long results - is just the same for both a^2 and 4kN
 * (and a multiple of 2^64).
 *  Thus <code>test</code> is correct and <code>b</code> is correct, too. <code>a</code> is correct anyway.
 *
 * @authors Tilman Neumann + Thilo Harich
 */
public class Lehman_Fast extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Lehman_FastOrig.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final TDiv63Inverse tdiv = new TDiv63Inverse(1<<21);

	private static double[] sqrt, sqrtInv;

	static {
		// Precompute sqrts for all possible k. 2^21 entries are enough for N~2^63.
		final int kMax = 1<<21;
		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0/sqrtI;
		}
	}

	private long N;
	private long fourN;
	private double sqrt4N;
	private final boolean doTDivFirst;
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Lehman_Fast(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
	}

	@Override
	public String getName() {
		return "Lehman_Fast(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	/**
	 * By looking at the upper limit of the inner 'a' loop we can directly see that below
	 * k limit / 64 the inner loop for 'a' can just iterate over two possible sucessing 'a'
	 * values. k limit here is n ^1/3.
	 * Additionally  we know from a mod 2 argument that only one out of the two values is possible.
	 * The loop above kLimit / 64 shrinks down to a simple loop over k - not over a.
	 * Additionally we ensure that we hit at least one 'a' value in the loop over the hole k range.
	 * In practice it is better to let the loop for the lower values of k run to n^1/3 / 8 instead of n^1/3 / 64.
	 * So we have a low k range k < n^1/3 / 8
	 * Mid range n^1/3 / 8 <= k <= n^1/3
	 * And a high range k >= n^1/3. The high range is not needed for the correctness and is outside of
	 * the Lehman numbers, but it is used for finding factors quickly
	 *
	 * We do the following steps
	 * <li> for number with small factors we use a trial division to find the small factors
	 * <li> iterate over k = 0 mod 6 and k = 3 mod 6 in the mid range, since they have a high chance to hit factors
	 * <li> iterate over k = 0 mod 6 and k = 3 mod 6 in the high range, since they have a high chance to hit factors
	 * <li> iterate in the low range
	 * <li> iterate over k != 0 and k != 3 mod 6 in the mid range theoretically needed, but it seems like we can skip it
	 * <li> use a trial division to find the small factors, if there are still any
	 * <li> find factors for numbers were the rounding performance improvement has failed
	 *
	 * The complicated numbers are the once in the upper middle range.
	 *
	 *
	 * @param N
	 * @return
	 */
	public long findSingleFactor(long N) {
		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		// do trial division before Lehman loop ?
		long factor;
		tdiv.setTestLimit(cbrt);
		if (doTDivFirst && (factor = tdiv.findSingleFactor(N))>1) return factor;

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);

		// kLimit must be 0 mod 6, since we also want to search above of it
		final int kLimit = ((cbrt + 6) / 6) * 6;
		// For kTwoA = kLimit / 64 the range for a is at most 2. We make it 0 mod 6, too.
		// interestingly it works a faster when we switch a the point with only one a possible
		final int kTwoA = (((cbrt >> 6) + 6) / 6) * 6;

		// We are investigating solutions of a^2 - sqrt(k*n) = y^2 in three k-ranges:
		// * The "small range" is 1 <= k < kTwoA, where we may have more than two 'a'-solutions per k.
		//   Thus, an inner 'a'-loop is required.
		// * The "middle range" is kTwoA <= k < kLimit, where we have at most two possible 'a' values per k.
		// * The "high range" is kLimit <= k < 2*kLimit. This range is not required for the correctness
		//   of the algorithm, but investigating it for some k==0 (mod 6) improves performance.

		// We start with the middle range cases k == 0 (mod 6) and k == 3 (mod 6),
		// which have the highest chance to find a factor. for small N gcd might return N
		final int kBegin = 3;
		if ((    factor = lehmanOdd315  (kBegin,        kLimit >> 1))      > 1 && factor < N ||
				(factor = lehmanEven    (kTwoA,         kLimit))      > 1 && factor < N ||
				(factor = lehmanOdd     (kTwoA + 3,     kLimit))      > 1 && factor < N )
			return factor;


		// Now investigate the small range
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability
		for (int k=1; k < kTwoA ; k++) {
			final double sqrt4kN = sqrt4N * sqrt[k];
			// only use long values
			final long aStart = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			long aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInv[k]);
			long aStep;
			if ((k & 1) == 0) {
				// k even -> make sure aLimit is odd
				aLimit |= 1L;
				aStep = 2;
			} else {
				aStep = ((k + N) & 3) == 0 ? 8 : 4;
				aLimit = adjustMod8(aLimit, k + N);
			}

			// we have increase the upper limit aLimit -> we can be sure that the a loop will be executed at least once
			// thus the time in determining the limits above is not wasted
			for (long a=aLimit; a >= aStart; a-=aStep) {
				final long test = a*a - k * fourN;
				// Here test<0 is possible because of double to long cast errors in the 'a'-computation.
				// But then b = Math.sqrt(test) gives 0 (sic!) => 0*0 != test => no errors.
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		if ((factor = lehmanEven (kLimit,        kLimit << 1)) > 1 && factor < N)
			return factor;

		// Complete middle range, according to lehman this is needed
		if ((    factor = lehmanOdd (kTwoA + 1, kLimit)) > 1 ||
				(factor = lehmanEven(kTwoA + 2, kLimit)) > 1 ||
				(factor = lehmanEven(kTwoA + 4, kLimit)) > 1 ||
				(factor = lehmanOdd (kTwoA + 5, kLimit)) > 1)
			return factor;



		// do trial division after Lehman loop ?
		if (!doTDivFirst && (factor = tdiv.findSingleFactor(N))>1)
			return factor;

		// If sqrt(4kN) is very near to an exact integer then the fast ceil() in the 'aStart'-computation
		// may have failed. Then we need a "correction loop":
		for (int k=kTwoA + 1; k <= kLimit; k++) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - k*fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}

		return 0; // fail
	}


	/**
	 * Find gcd(a+b) for a^2 - k* N = b^2, where k is odd.
	 * In this case a has to be even and k + n == a mod 4.
	 * @param kBegin
	 * @param kEnd
	 * @return
	 */
	private long lehmanOdd(int kBegin, final int kEnd) {
		for (int k = kBegin; k <= kEnd; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else {
				final long adjust1 = (kPlusN - a) & 15;
				final long adjust2 = (-kPlusN - a) & 15;
				a += adjust1<adjust2 ? adjust1 : adjust2;
			}
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}
	private long lehmanOdd315(int kBegin, final int kEnd) {
		//		final int multiplier = 3 * 3 * 5 * 7 * 11 * 13;
		//		final int multiplier = 9 * 5 * 7;
		final int multiplier = 9 * 25;
		final double sqrt1260N = Math.sqrt(multiplier) * sqrt4N;
		final long N1260 = 4 * multiplier * N;
		for (int k = kBegin; k <= kEnd; k += 6) {
			long a = (long) (sqrt1260N * sqrt[k] + ROUND_UP_DOUBLE);
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else {
				final long adjust1 = (kPlusN - a) & 15;
				final long adjust2 = (-kPlusN - a) & 15;
				a += adjust1<adjust2 ? adjust1 : adjust2;
			}
			final long test = a*a - k * N1260;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}

	/**
	 * Find gcd(a+b) for a^2 - k* N = b^2, where k is even.
	 * In this case a has to be odd.
	 * @param kBegin
	 * @param kEnd
	 * @return
	 */
	private long lehmanEven315(int kBegin, final int kEnd) {
		//		final int multiplier = 3 * 3 * 5 * 7 * 11 * 13;
		final int multiplier = 5;
		final double sqrt1260N = Math.sqrt(multiplier) * sqrt4N;
		final long N1260 = 4 * multiplier * N;
		for (int k = kBegin; k <= kEnd; k += 6) {
			// k even -> a must be odd
			final long a = (long) (sqrt1260N  * sqrt[k] + ROUND_UP_DOUBLE) | 1L;
			final long test = a*a - k * N1260 ;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}
	/**
	 * Find gcd(a+b) for a^2 - k* N = b^2, where k is even.
	 * In this case a has to be odd.
	 * @param kBegin
	 * @param kEnd
	 * @return
	 */
	private long lehmanEven(int kBegin, final int kEnd) {
		for (int k = kBegin; k <= kEnd; k += 6) {
			// k even -> a must be odd
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1L;
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}



	/**
	 * Here we make sure that for k even a fulfills
	 * a^2 - k*N = y^2 mod 8.
	 * k + n == a mod 8 if k+n == 0 mod 4 or else
	 * k + n == a mod 4
	 * @param a
	 * @param kPlusN k + N
	 * @return
	 */
	private static long adjustMod8(long a, long kPlusN) {
		if ((kPlusN & 3) == 0) {
			return a + ((kPlusN - a) & 7);
		} else {
			final long adjust1 = (kPlusN - a) & 15;
			final long adjust2 = (-kPlusN - a) & 15;
			return a + adjust1<adjust2 ? adjust1 : adjust2;
		}

	}
}
