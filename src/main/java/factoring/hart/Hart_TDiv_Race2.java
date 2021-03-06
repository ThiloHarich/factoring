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
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.jml.primes.exact.AutoExpandingPrimesArray;
import de.tilman_neumann.util.ConfigUtil;

/**
 * A factoring algorithm racing Hart's one line factorizer against trial division.
 *
 * This is the fastest algorithm for test numbers <= 44 bit when their nature is unknown.
 * Hart_TDiv_Race_Unsafe is faster for moderate semiprimes N >= 45 bit.
 * Hart_Fast and Hard_Fast_Unsafe are faster for hard semiprimes.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class Hart_TDiv_Race2 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_TDiv_Race2.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT =3;

	/** Size of arrays, sufficient to factor all numbers <= 52 bit. */
	private static final int I_MAX = 1<<20;

	private static final int DISCRIMINATOR_BITS = 10; // experimental result
	private static final double DISCRIMINATOR = 1.0/(1<<DISCRIMINATOR_BITS);

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final double[] sqrt;
	private final int[] primes;
	private final double[] reciprocals;

	private final int [] divisors = new int []{3,3,7};
	private final boolean[] smooth =  new boolean [3*3*7]; // 63, 2431, 46189

	private final AutoExpandingPrimesArray SMALL_PRIMES = AutoExpandingPrimesArray.get();
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 */
	public Hart_TDiv_Race2() {
		sqrt = new double[I_MAX];
		primes = new int[I_MAX];
		reciprocals = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt[i] = Math.sqrt(i*K_MULT);
			final int p = SMALL_PRIMES.getPrime(i);
			primes[i] = p;
			reciprocals[i] = 1.0/p;
		}
		for (int i=1; i<smooth.length; i++) {
			for (final int divisor : divisors){
				if ( i% divisor == 0)
					smooth[i] = true;
			}
		}
	}

	@Override
	public String getName() {
		return "Hart_Thilo";
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
		final int pMinBits = 21 - Long.numberOfLeadingZeros(N);

		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		int k = K_MULT;
		try {
			int i=1;
			if (pMinBits>0) {
				// for the smallest primes we must do standard trial division
				final int pMin = 1<<pMinBits;
				for ( ; primes[i]<pMin; i++, k += K_MULT) {
					// tdiv step
					if (N%primes[i]==0) {
						return primes[i];
					}

					a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
					a = adjustA(N, a, k, i);
					test = a*a - k * fourN;
					b = (long) Math.sqrt(test);
					if (b*b == test) {
						if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) return gcd;
					}
				}
			}

			// continue with Hart and fast inverse trial division
			for (; ; i++, k += K_MULT) {
				// tdiv step
				//LOG.debug("test p[" + i + "] = " + primes[i]);
				final long nDivPrime = (long) (N*reciprocals[i] + DISCRIMINATOR);
				if (nDivPrime * primes[i] == N) {
					// nDivPrime is very near to an integer
					if (N%primes[i]==0) {
						//LOG.debug("Found factor " + primes[i]);
						return primes[i];
					}
				}

				// odd k -> adjust a mod 8
				a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				a = adjustA(N, a, k, i);
				test = a*a - k * fourN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) return gcd;
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			LOG.error("Hart_TDiv_Race: Failed to factor N=" + N + " because the arrays are too small.");
			return 0;
		}
	}

	private long adjustA(long N, long a, int k, int i) {
		if ((i & 1) == 0)
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
				16032993843L, // fine here
				26036808587L,
				41703657595L, // fine here
				68889614021L,
				197397887859L, // fine here

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
				9
		};

		final Hart_TDiv_Race2 holf = new Hart_TDiv_Race2();
		for (final long N : testNumbers) {
			final long factor = holf.findSingleFactor(N);
			LOG.info("N=" + N + " has factor " + factor);
		}
	}
}
