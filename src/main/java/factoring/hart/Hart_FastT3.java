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
package factoring.hart;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.ConfigUtil;

/**
 * Pretty simple yet fast variant of Hart's one line factorizer.
 *
 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
 * But the smaller possible factors get, it will become slower and slower.
 *
 * For any kind of test numbers except very hard semiprimes, Hart_TDiv_Race will be faster.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class Hart_FastT3 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_FastT3.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 1; // 315

	/** Size of arrays */
	private static final int I_MAX = 1<<20;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Hart_FastT3(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt[i] = Math.sqrt(i*K_MULT);
		}
	}

	@Override
	public String getName() {
		return "Hart_FastT(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	/**
	 * Find a factor of long N.
	 * @param N
	 * @return factor of N
	 */
	public long findSingleFactor(long N) {
		// do trial division before the Hart loop ?
		long factor;
		if (doTDivFirst) {
			tdiv.setTestLimit((int) Math.cbrt(N));
			if ((factor = tdiv.findSingleFactor(N))>1) return factor;
		} else {
			// at least do trial division by multiplier factors
			if ((N%3)==0) return 3;
			if ((N%5)==0) return 5;
			if ((N%7)==0) return 7;
		}

		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		try {
			for (int kFactorN=1; ;kFactorN++) {
				// 1 , 3
				for (int kFactor1=1; kFactor1 <= 3 ;kFactor1 += 3) {
					// 5
					//					for (int kFactor2=5; kFactor2 <= 11 ;kFactor2 += 2) {
					final int kFactor2=5;
					// 5, 15
					//						for (int kFactor3=5; kFactor3 <= 9 ;kFactor3 += 2) {
					final int k = kFactor1 * kFactor2  * kFactorN;
					// odd k -> adjust a mod 8
					a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
					a = adjustA(N, a, k);
					test = a*a - k * fourN;
					b = (long) Math.sqrt(test);
					if (b*b == test) {
						if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {

							return gcd;
						}
					}
					//						}
					//					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// this may happen if this implementation is tested with doTDivFirst==false and N having
			// very small factors, or if N is too big
			return 0;
		}
	}

	private long adjustA(long N, long a, int k) {
		if ((k & 1) == 0)
			a |= 1;
		else {
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

				// test numbers that required large arrays
				135902052523483L,
				1454149122259871L,
				5963992216323061L,
				26071073737844227L,
				8296707175249091L,
				35688516583284121L,
				//35245060305489557L, // too big for I_MAX
				//107563481071570333L, // too big for I_MAX
				//107326406641253893L, // too big for I_MAX
				//120459770277978457L, // too big for I_MAX

				// failures with random odd composites
				949443, // = 3 * 11 * 28771
				996433, // = 31 * 32143
				1340465, // = 5 * 7 * 38299
				1979435, // = 5 * 395887
				2514615, // = 3 * 5 * 167641
				5226867, // =  3^2 * 580763
				10518047, // = 61 * 172427
				30783267, // = 3^3 * 1140121
				62230739, // = 67 * 928817
				84836647, // = 7 * 17 * 712913
				94602505,
				258555555,
				436396385,
				612066705,
				2017001503,
				3084734169L,
				6700794123L,
				16032993843L, // = 3 * 5344331281, 34 bit number
				26036808587L,
				41703657595L, // = 5 * 8340731519, 36 bit number
				68889614021L,
				197397887859L, // = 3^2 * 21933098651, 38 bit number

				2157195374713L,
				8370014680591L,
				22568765132167L,
				63088136564083L,

				// more test numbers with small factors
				// 30 bit
				712869263, // = 89 * 8009767
				386575807, // = 73 * 5295559
				569172749, // = 83 * 6857503
				// 40 bit
				624800360363L, // = 233 * 2681546611
				883246601513L, // = 251 * 3518910763

				// problems found by Thilo
				35184372094495L,
				893, // works
				35, // works
				9 // works
		};

		final Hart_FastT3 holf = new Hart_FastT3(false);
		for (final long N : testNumbers) {
			final long factor = holf.findSingleFactor(N);
			LOG.info("N=" + N + " has factor " + factor);
		}
	}
}
