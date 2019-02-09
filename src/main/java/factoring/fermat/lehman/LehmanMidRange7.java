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
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.trial.TrialMultiplyCorrection;

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
public class LehmanMidRange7 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Lehman_FastOrig.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final TrialMultiplyCorrection trialDivision = new TrialMultiplyCorrection(1<<21);

	private static double[] sqrt, sqrtInv;

	public static int [] foundInStep = new int[20];
	// 5 * 7 = 35
	private static final double sqrt35    = Math.sqrt(5 * 7);
	// 5 * 7 * 11 = 385
	private static final double sqrt385   = Math.sqrt(5 * 7 * 11);
	// 5 * 7 * 11 * 13 = 5005
	private static final double sqrt5005  = Math.sqrt(5 * 7 * 11 * 13);
	private static final double sqrt85085  = Math.sqrt(85085);


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
	private final int trialDivisionWhen;
	private final int numMultipliers;
	private final Gcd63 gcdEngine = new Gcd63();


	/**
	 * @param doTrialDivisionFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < n^1/3 frequently.
	 * If you have numbers which are semiprimes where the 2 or 3 factors are
	 * all above n^1/3 set it to false. In this case the trial division is done after the
	 * Leman loop.
	 * @param factorSize depending on the size of the second largest factor the algorithm can be steered
	 * to show best performance. If there might be a factor lower then n^1/3 then 0 is the best value.
	 * Then the trial division for numbers up to n^1/3 is done right at the beginning.
	 * If the numbers is a semiprime and one factor is known to be slightly higher then n^1/3.
	 * choose a value 1. Trial division is then applied as a last step.
	 * If the second highest factor is much greater then n^1/3 choose 2,
	 * and we apply a loop with an additional multiplier. If factors are known to be around n^1/2
	 * choose value 3.
	 */
	public LehmanMidRange7(int trialDivisionWhen, int numMultipliers) {
		this.trialDivisionWhen = trialDivisionWhen;
		this.numMultipliers = numMultipliers;
	}

	@Override
	public String getName() {
		return "LehmanMidRange(" + trialDivisionWhen + ", " + numMultipliers + ")";
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
	 *
	 *
	 * We do the following steps
	 * <li> for number with small factors we use a trial division to find the small factors
	 * <li> iterate over k = 0 mod 6 until k <= cbrt(n)/3 -> main loop for finding numbers
	 * <li> check  k = 0 mod 6 * product of small primes in parallel -> works good for high factors
	 * <li> do trial division for medium factors
	 * <li> iterate over k = 3 mod 6 until k <= cbrt(n) -> investigate in even a
	 * <li> iterate in the low range k < cbrt(n) / 64
	 * <li> iterate over k != 0 and k != 3 mod 6 in the mid range theoretically needed, but needed very rarely
	 * <li> use a trial division to find the small factors, if there are still any
	 * <li> find factors for numbers were the rounding performance improvement has failed
	 *
	 * @param N
	 * @return
	 */
	public long findSingleFactor(long N) {
		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		// do trial division before Lehman loop ?
		long factor;
		trialDivision.setMaxFactor(cbrt);
		if (trialDivisionWhen == 0 && (factor = trialDivision.findFactor(N))>1) {
			foundInStep[0]++;
			return factor;
		}

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);

		final int multiplier = 6;
		// kLimit must be 0 mod 6, since we also want to search above of it
		int kLimit = ((cbrt + multiplier) / multiplier) * multiplier;
		final int twoA = cbrt / 64;
		final int kTwoA = ((twoA + multiplier) /multiplier) * multiplier;
		kLimit = Math.max(kLimit, kTwoA);
		final int kMod6Split = (kLimit /18) * 6;

		// We start with the middle range cases k == 0 mod 6 and for bigger factors mod 6*5*7 or mod 6*5*7*11
		if ((factor = lehmanKEven (kTwoA, kMod6Split, numMultipliers)) > 1 && factor < N)
			return factor;


		// Now investigate in k = 3 mod 6 which have the next highest chance for big factors
		// handles ~ 2% of the cases (45 bit)
		if ((factor = lehmanKOdd  (kTwoA + 3, kLimit)) > 1 && factor < N ) {
			foundInStep[6]++;
			return factor;
		}

		if (trialDivisionWhen == 1 && (factor = trialDivision.findFactor(N))>1) {
			foundInStep[7]++;
			return factor;
		}

		// Now investigate the small range
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability
		for (int k=1; k < kTwoA ; k++) {
			final double sqrt4kN = sqrt4N * sqrt[k];
			// only use long values
			long a = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			final long aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInv[k]);
			long aStep;
			if ((k & 1) == 0) {
				//				 k even -> make sure aLimit is odd
				a |= 1L;
				aStep = 2;
			} else {
				aStep = ((k + N) & 3) == 0 ? 8 : 4;
				a = adjustMod8(a, k + N);
			}

			// make sure we have that the a loop will be executed at least once
			// thus the time in determining the limits above is not wasted
			do {
				final long test = a*a - k * fourN;
				// Here test<0 is possible because of double to long cast errors in the 'a'-computation.
				// But then b = Math.sqrt(test) gives 0 (sic!) => 0*0 != test => no errors.
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					foundInStep[8]++;
					return gcdEngine.gcd(a+b, N);
				}
				a += aStep;
			} while (a <= aLimit);

		}
		if (trialDivisionWhen == 2 && (factor = trialDivision.findFactor(N))>1) {
			foundInStep[9]++;
			return factor;
		}

		// also for factors above n^1/3 we look in the higher range, but after the regular Lehman step
		if ((factor = lehmanKEven (kMod6Split, kLimit << 1, 0)) > 1 && factor < N) {
			foundInStep[10]++;
			return factor;
		}

		// now handle very seldom happening cases
		// Complete middle range, according to lehman this is needed.
		if (    (factor = lehmanKOdd  (kTwoA + 1, kLimit)) > 1 ||
				(factor = lehmanKEven (kTwoA + 2, kLimit, 0)) > 1 ||
				(factor = lehmanKEven (kTwoA + 4, kLimit, 0)) > 1 ||
				(factor = lehmanKOdd  (kTwoA + 5, kLimit)) > 1) {
			foundInStep[11]++;
			return factor;
		}

		// If sqrt(4kN) is very near to an exact integer then the fast ceil() in the 'aStart'-computation
		// may have failed. Then we need a "correction loop":
		for (int k=kTwoA + 1; k <= kLimit; k++) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - k*fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				foundInStep[12]++;
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
	private long lehmanKOdd(int kBegin, final int kEnd) {
		for (int k = kBegin; k <= kEnd; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			a = adjustMod8(a, k + N);
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}

	/**
	 * Main lehman based loop for finding factors.
	 * Find gcd(a+b) for a^2 - k* N = b^2, where k is a multiple of 6.
	 * In this case a has to be odd.
	 * For numbers with two factors of the same size we additionally check for multiples of 6 like 6 * (5*7),
	 * 6 * (5*7*11) and 6 * (5*7*11*13).
	 * This brings  40% performance (at 45 bits) for semiprimes with factors of the same size (in bits).
	 * @param kBegin
	 * @param kEnd
	 * @param factorSize TODO
	 * @return
	 */
	private long lehmanKEven(int kBegin, final int kEnd, int numMultipliers) {
		int k = kBegin;
		long a,b, test;
		for (; k <= kEnd; k += 6) {
			final double sqrt4kn = sqrt4N * sqrt[k];
			final long fourkn = k * fourN;

			// if the factor is in the n^1/3 range we focus on numbers in the range
			a = (long) (sqrt4kn + ROUND_UP_DOUBLE) | 1L;
			test = a*a - fourkn;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				foundInStep[1]++;
				return gcdEngine.gcd(a+b, N);
			}

			// if we have a factor bigger then n^1/3 we can find more factors when looking
			// into multiples of 6, even more
			if (numMultipliers >= 2) {
				//				a = (long) (sqrt5005* sqrt4kn + ROUND_UP_DOUBLE) | 1L;
				//				test = a*a - 5005 * fourkn;
				a = (long) (sqrt385 * sqrt4kn + ROUND_UP_DOUBLE) | 1L;
				test = a*a - 385 * fourkn;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					foundInStep[3]++;
					return gcdEngine.gcd(a+b, N);
				}
			}
			// if we have a factor bigger then n^1/3 we can find more factors when looking
			// into multiples of 6, even more
			if (numMultipliers >= 3) {
				//				a = (long) (sqrt385 * sqrt4kn + ROUND_UP_DOUBLE) | 1L;
				//				test = a*a - 385 * fourkn;
				a = (long) (sqrt5005* sqrt4kn + ROUND_UP_DOUBLE) | 1L;
				test = a*a - 5005 * fourkn;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					foundInStep[4]++;
					return gcdEngine.gcd(a+b, N);
				}
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
		if ((kPlusN & 3) == 0)
			return a + ((kPlusN - a) & 7);
		return a + ((kPlusN - a) & 3);
	}

}
