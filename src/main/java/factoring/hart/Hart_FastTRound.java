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

/**
 * Pretty simple yet fast variant of Hart's one line factorizer.
 * Compared with the regular Lehman algorithm, the Hart algorithm does not
 * need a second loop to iterate over the numbers 'a' for a given 'k' in the equation a^2 - 4kn.
 * So the upper bound for this does not has to be calculated. For each k the value sqrt(k) in order
 * to determine a = ceil(sqrt(4kn))
 * will be calculated only once and then stored in an array. This speeds up the sieving buy
 * a big factor since calculating the sqrt is expensive.
 * We choose k to be a multiple of 315 = 3*3*5*7 this causes that
 * a^2 - 4kn = b^2 mod 3*3*5*7 which increases the chance to find a solution a^2 - 4kn = b^2 pretty much.
 * After selecting 'a' we ensure that a^2 - 4kn = b^2 mod 64 by increasing 'a' at most by 8.
 *
 *
 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
 * But the smaller possible factors get, it will become slower and slower.
 *
 * For any kind of test numbers except very hard semiprimes, Hart_TDiv_Race will be faster.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class Hart_FastTRound extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_FastTRound.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 1; // 315

	/** Size of arrays this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	//	static int mod = 11 * 13 * 17 * 2;
	static int mod = 11*13*17*19;
	//	static int mod2 = mod * 2;
	//	static int mod4 = mod * 4;
	static boolean[] squares = new boolean [mod];

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();


	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Hart_FastTRound(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt[i] = Math.sqrt(i*K_MULT);
		}
		for (int i=0; i<mod; i++) {
			squares[(i*i) % mod] = true;
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
		if (doTDivFirst) {
			// do trial division before the Hart loop
			tdiv.setTestLimit((int) Math.cbrt(N));
			final long factor = tdiv.findSingleFactor(N);
			if (factor > 1) return factor;
		}


		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		int k = K_MULT;
		try {
			for (int i=1; ;i++, k += K_MULT) {
				// calculating the sqrt here is 5 times slower then storing it
				a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				if ((i & 1) == 0) {
					a |= 1;
					// find the next odd a
					//					test = a*a - k * fourN;
					//					int bMod = (int) (test % mod);
					//					if (!squares[bMod]) {
					//						int aMod = (int)(a % mod);
					//						do {
					//							// a = a + 2
					//							// (a+2)^2 = a^2 + 4a + 4 = a^2 + 4*(a+1)
					//							aMod = aMod == mod-1 ? 0 : aMod +1;
					//							a += 2;
					//
					//							bMod += (aMod) << 2;
					//							bMod = bMod >= mod4 ? bMod - mod4 : bMod;
					//							bMod = bMod >= mod2 ? bMod - mod2 : bMod;
					//							bMod = bMod >= mod  ? bMod - mod  : bMod;
					//						}while (!squares[bMod]);
					//					}
				}
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
				test = a*a - k * fourN;
				// For even i with a odd we try to find a good a modulo mod
				a = ((i & 1) == 0) && !squares[(int) (test % mod)] ? a+2 : a;
				//				if ((i & 1) == 0) {
				//					final int testMod = (int) (test % mod);
				//					//					while (!squares[testMod]) {
				//					if (!squares[testMod]) {
				//						a += 2;
				//						test = a*a - k * fourN;
				//						//						testMod = (int) (test % mod);
				//					}
				//				}
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
						return gcd;
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			//			e.printStackTrace();
			LOG.error("Hart_Fast: Failed to factor N=" + N + ". Either it has factors < cbrt(N) needing trial division, or the arrays are too small.");
			// TODO Or N is square
			return 1;
		}
	}
}
