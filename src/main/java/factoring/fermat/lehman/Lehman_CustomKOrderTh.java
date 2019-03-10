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
import de.tilman_neumann.util.ConfigUtil;

/**
 * A variant of Lehman's algorithm that allows to arrange the k's in arrays of different priorities.
 * Testing multiples of 15 first, followed by multiples of 3, then the rest works not so bad...
 *
 * @authors Tilman Neumann + Thilo Harich
 */
public class Lehman_CustomKOrderTh extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Lehman_CustomKOrderTh.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final int K_MAX = 1<<20;
	final int maxLevel = 5;


	private final TDiv63Inverse tdiv = new TDiv63Inverse(K_MAX);

	private final double[] sqrts;
	private final double[] sqrtInvs;
	private final int[] ks;
	private final double[] sqrts1;
	private final double[] sqrtInvs1;
	private final int[] ks1;

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
	public Lehman_CustomKOrderTh(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// arrange k in different arrays
		sqrts = new double[K_MAX+1];
		sqrtInvs = new double[K_MAX+1];
		ks = new int [K_MAX+1];
		sqrts1 = new double[K_MAX+1];
		sqrtInvs1 = new double[K_MAX+1];
		ks1 = new int [K_MAX+1];

		for (int pos = 0; pos < K_MAX; pos ++) {
			ks[pos] = 105 * (pos+1);
			sqrts[pos] = Math.sqrt(105 * (pos+1));
			sqrtInvs[pos] = 1 / sqrts[pos];
			ks1[pos] = pos+1;
			sqrts1[pos] = Math.sqrt(pos+1);
			sqrtInvs1[pos] = 1 / sqrts1[pos];
		}
	}

	@Override
	public String getName() {
		return "Lehman_CustomKOrder(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		// N==9 would require to check if the gcd is 1 < gcd < N before returning it as a factor
		if (N==9) return 3;

		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		// do trial division before Lehman loop ?
		long factor;
		tdiv.setTestLimit(cbrt);
		if (doTDivFirst && (factor = tdiv.findSingleFactor(N))>1) return factor;

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);

		final int kLimit = cbrt;
		// For kTwoA = kLimit / 64 the range for a is at most 2. kLimit / 128 seems to work as well...
		final int kTwoA = (cbrt + 127) >> 7;

		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability
		if ((factor = test(kTwoA, 8 *kLimit, sixthRootTerm, ks, sqrts, sqrtInvs)) > 1) return factor;

		// do trial division now?
		if (!doTDivFirst && (factor = tdiv.findSingleFactor(N))>1) return factor;

		// finish Lehman loops
		if ((factor = test(kTwoA, kLimit, sixthRootTerm, ks1, sqrts1, sqrtInvs1)) > 1) return factor;

		// If sqrt(4kN) is very near to an exact integer then the fast ceil() in the 'aStart'-computation
		// may have failed. Then we need a "correction loop":
		if ((factor = correctionLoop(kLimit)) > 1) return factor;

		return 1; // fail
	}

	private long test(int kTwoA, int kLimit, double sixthRootTerm, int[] ksx, double[] sqrtsx, double[] sqrtInvsx) {
		long aLimit, aStart;
		int i, k;
		long gcd;
		for (i=0; (k = ksx[i]) < kTwoA; i++) {
			final double sqrt4kN = sqrt4N * sqrtsx[i];
			aStart = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInvsx[i]);
			aLimit = adjustA(aLimit, k);

			// processing the a-loop top-down is faster than bottom-up
			final long fourkN = k * fourN;
			for (long a=aLimit; a >= aStart; a-=2) {
				final long test = a*a - fourkN;
				// Here test<0 is possible because of double to long cast errors in the 'a'-computation.
				// But then b = Math.sqrt(test) gives 0 (sic!) => 0*0 != test => no errors.
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N)
						return gcd;
				}
			}
		}

		for ( ; i <kLimit; i++) {
			k = ksx[i];
			long a = (long) (sqrt4N * sqrtsx[i] + ROUND_UP_DOUBLE);
			a = adjustA(a, k);
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N)
					return gcd;
			}
		}
		return 1;
	}

	private long adjustA(long a, int k) {
		if ((k & 1) == 0) {
			// k even -> make sure aLimit is odd
			a |= 1L;
		} else {
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else {
				final long adjust1 = (kPlusN - a) & 15;
				final long adjust2 = (-kPlusN - a) & 15;
				a += adjust1<adjust2 ? adjust1 : adjust2;
			}
		}
		return a;
	}

	private long correctionLoop(int kLimit) {
		int i=0, k;
		for (; (k = ks1[i])<kLimit; i++) {
			final long a = (long) (sqrt4N * sqrts1[i] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - k*fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return 1;
	}

	/**
	 * Test.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		ConfigUtil.initProject();

		// These test number were too hard for previous versions:
		final long[] testNumbers = new long[] {
				// odd semiprimes
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

				// special case
				9,
		};

		final Lehman_CustomKOrderTh lehman = new Lehman_CustomKOrderTh(false);
		for (final long N : testNumbers) {
			final long factor = lehman.findSingleFactor(N);
			LOG.info("N=" + N + " has factor " + factor);
		}
	}
}
