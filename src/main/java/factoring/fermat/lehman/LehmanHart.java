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

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.ConfigUtil;
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
public class LehmanHart extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_FastOrig.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final TrialMultiplyCorrection trialDivision = new TrialMultiplyCorrection(1<<21);

	private static double[] sqrt, sqrtInv;

	// 5 * 7 = 35
	private static final double sqrt35    = Math.sqrt(5 * 7);
	// 5 * 7 * 11 = 385
	private static final double sqrt385   = Math.sqrt(5 * 7 * 11);
	// 5 * 7 * 11 * 13 = 5005
	private static final double sqrt5005  = Math.sqrt(5 * 7 * 11 * 13);
	// 5 * 7 * 11 * 13 * 17 = 5005
	private static final double sqrt85085 = Math.sqrt(5 * 7 * 11 * 13 * 17);

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
	private final int trialPhase;
	private final Gcd63 gcdEngine = new Gcd63();


	/**
	 * Full constructor.
	 * @param doTrialDivisionFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < n^1/3 frequently.
	 * If you have numbers which are semiprimes where the 2 or 3 factors are
	 * all above n^1/3 set it to false. In this case the trial division is done after the
	 * Leman loop.
	 * @param trialPhase depending on the size of the second largest factor the algorithm can be steered
	 * to show best performance. It all factors are known to be lower then n^1/3 then 0 is the best value.
	 * Then the trial division for numbers up to n^1/3 is done right at the beginning.
	 * This usually is also the best value for factoring ascending or completely random numbers.
	 * If the numbers is a semiprime and the lowest factor is known to be higher then n^1/3.
	 * choose a value 2 or higher. Trial division is then applied as a last step.
	 * If numbers have one high factor > n^1/3, but the second largest factor might be lower then n^1/3
	 * choose 1. Then trial division is done in the middle. This is usually the best choice if you do not know
	 * anything about the numbers.
	 */
	public LehmanHart(int trialPhase) {
		this.trialPhase = trialPhase;
	}

	@Override
	public String getName() {
		return "LehmanMidRange(" + trialPhase + ")";
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
	 * So we have a low k range k < n^1/3 / 64
	 * Mid range n^1/3 / 64 <= k <= n^1/3
	 * And a high range k >= n^1/3. The high range is not needed for the correctness and is outside of
	 * the Lehman numbers, but it is used for finding factors quickly.
	 * If hard numbers have a small factor near n^1/3 we want find these factors by trial division.
	 * We increase the point when switching from the Lehman phase to the trial division phase by
	 * a small multiplier (like 1.4). This decrease the limit of k (the multiplier of N) by the
	 * square of this multiplier, but increases the inner loop over 'a' by the square of this multiplier.
	 * Before switching to find factors in the lower range, we investigate in good numbers in the high range.
	 *
	 *
	 * We do the following steps
	 * <li> for number with small factors we use a trial division to find the small factors
	 * <li> iterate over k = 0 mod 6 in the mid range and k = 0 mod 30 in mid and high range in parallel
	 * <li> iterate over k = 0 mod 6 and k = 3 mod 6 in the mid range, since they have a high chance to hit factors
	 * <li> iterate over k = 0 mod 6 and k = 3 mod 6 in the high range, since they have a high chance to hit factors
	 * <li> iterate in the low range
	 * <li> iterate over k != 0 and k != 3 mod 6 in the mid range theoretically needed, but it seems like we can skip it
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
		if (trialPhase == 0 && (factor = trialDivision.findFactor(N))>1)
			return factor;

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);


		// We start with the middle range cases k == 0 mod 6 and mod 6*5*7 and mod 6*5*7*11 ,
		// which have the highest chance to find a factor. For factors around n^1/2 ~ 90% of the cases.
		// Then we investigate in k = 3 mod 6. Then there are only very few cases open for big factors.
		//		for small N gcd might return N
		if ((factor = lehman (1)) > 1 && factor < N)
			return factor;


		return 0; // fail
	}



	private long lehman(int kBegin) {
		long factor = 1;
		for (int k = kBegin; ;) {
			if (factor == 1 || factor == N)
				factor = findFactorEven(6*k, factor);
			if (factor == 1 || factor == N)
				factor = findFactorOdd(k, factor);
			k++;
			if (factor == 1 || factor == N)
				factor = findFactorEven(6*k, factor);
			if (factor == 1 || factor == N)
				factor = findFactorEven(k, factor);

			if (factor > 1 && factor < N)
				return factor;
			k++;
		}
	}

	private long findFactorOdd(int k, long factor) {
		long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
		a = adjustMod8(a, k + N);
		final long test = a*a - (k) * fourN;
		final long b = (long) Math.sqrt(test);
		if (b*b == test) {
			factor = gcdEngine.gcd(a+b, N);
		}
		return factor;
	}

	private long findFactorEven(int k, long factor) {
		final double sqrt4kn = sqrt4N * sqrt[k];
		final long fourkn = (k) * fourN;
		long a = (long) (sqrt4kn + ROUND_UP_DOUBLE) | 1L;
		long test = a*a - fourkn;
		long b = (long) Math.sqrt(test);
		if (b*b == test) {
			factor = gcdEngine.gcd(a+b, N);
		}
		a = (long) (sqrt5005 * sqrt4kn + ROUND_UP_DOUBLE) | 1L;
		test = a*a - 5005 * fourkn;
		b = (long) Math.sqrt(test);
		if (b*b == test && (factor == 1 || factor == N)) {
			factor = gcdEngine.gcd(a+b, N);
		}
		// and k'' = 11 * k' = 385 * k  = 11 * 7 * 5 * k
		a = (long) (sqrt385 * sqrt4kn + ROUND_UP_DOUBLE) | 1L;
		test = a*a - 385 * fourkn;
		b = (long) Math.sqrt(test);
		if (b*b == test && (factor == 1 || factor == N)) {
			factor = gcdEngine.gcd(a+b, N);
		}
		return factor;
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

	/**
	 * Test.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		ConfigUtil.initProject();

		// These test number were too hard for previous versions:
		final long[] testNumbers = new long[] {
				5640012124823L,
				7336014366011L,
				19699548984827L,
				52199161732031L,
				73891306919159L,
				112454098638991L,

				32427229648727L,
				87008511088033L,
				92295512906873L,
				338719143795073L,
				346425669865991L,
				1058244082458461L,
				1773019201473077L,
				6150742154616377L,

				44843649362329L,
				67954151927287L,
				134170056884573L,
				198589283218993L,
				737091621253457L,
				1112268234497993L,
				2986396307326613L,

				26275638086419L,
				62246008190941L,
				209195243701823L,
				290236682491211L,
				485069046631849L,
				1239671094365611L,
				2815471543494793L,
				5682546780292609L,
		};

		final LehmanHart lehman = new LehmanHart(1);
		for (final long N : testNumbers) {
			final long factor = lehman.findSingleFactor(N);
			LOG.info("N=" + N + " has factor " + factor);
		}
	}
}
