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
 * need a second loop to iterate over the numbers a for a given k in the equation a^2 - 4kn.
 * So the upper bound for this does not has to be calculated. The value for a = ceil(sqrt(4kn))
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
public class Hart_FastFMA extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_FastFMA.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	private static final int K4_MULT = 4*3*3*5*7; // 315
	//	private static final int K_MULT = 1; // 315

	/** Size of arrays this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	public long fma=0;
	public long adjust = 0;
	public long sqrtTime = 0;
	public long test = 0;


	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Hart_FastFMA(boolean doTDivFirst) {
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
		if (doTDivFirst) {
			// do trial division before the Hart loop
			tdiv.setTestLimit((int) Math.cbrt(N));
			final long factor = tdiv.findSingleFactor(N);
			if (factor > 1) return factor;
		}


		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long b, gcd;
		long a0, a1, a2, a3, test0, test1, test2, test3;
		long k0 =   K_MULT;
		long k1 = 2*K_MULT;
		long k2 = 3*K_MULT;
		long k3 = 4*K_MULT;
		try {
			for (int i=1; ; ) {
				long split1 = System.currentTimeMillis();
				// calculating the sqrt here is 5 times slower then storing it
				//				a0 = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE);
				//				a1 = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE);
				//				a2 = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE);
				//				a3 = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE);
				a0 = (long) Math.fma(sqrt4N, sqrt[i++], ROUND_UP_DOUBLE);
				a1 = (long) Math.fma(sqrt4N, sqrt[i++], ROUND_UP_DOUBLE);
				a2 = (long) Math.fma(sqrt4N, sqrt[i++], ROUND_UP_DOUBLE);
				a3 = (long) Math.fma(sqrt4N, sqrt[i++], ROUND_UP_DOUBLE);

				//				long split2 = System.nanoTime();
				long split2 = System.currentTimeMillis();
				fma += split2 - split1;

				a0 = adjustA(fourN, a0, (int) k0);
				a1 = adjustA(fourN, a1, (int) k1);
				a2 = adjustA(fourN, a2, (int) k2);
				a3 = adjustA(fourN, a3, (int) k3);

				split1 = System.currentTimeMillis();
				adjust += split1 - split2;

				// we need to do this here since the increased k is needed
				test0 = a0 * a0 - k0 * fourN;
				test1 = a1 * a1 - k1 * fourN;
				test2 = a2 * a2 - k2 * fourN;
				test3 = a3 * a3 - k3 * fourN;

				split2 = System.currentTimeMillis();
				test += split2 - split1;

				b = (long) Math.sqrt(test0);
				if (b*b == test0) {
					if ((gcd = gcdEngine.gcd((a0)+b, N))>1 && gcd<N) {
						return gcd;
					}
				}

				b = (long) Math.sqrt(test1);
				if (b*b == test1) {
					if ((gcd = gcdEngine.gcd((a1)+b, N))>1 && gcd<N) {
						return gcd;
					}
				}

				b = (long) Math.sqrt(test2);
				if (b*b == test2) {
					if ((gcd = gcdEngine.gcd((a2)+b, N))>1 && gcd<N) {
						return gcd;
					}
				}

				b = (long) Math.sqrt(test3);
				if (b*b == test3) {
					if ((gcd = gcdEngine.gcd((a3)+b, N))>1 && gcd<N) {
						return gcd;
					}
				}
				split1 = System.currentTimeMillis();
				sqrtTime += split1 - split2;

				k0 += K4_MULT;
				k1 += K4_MULT;
				k2 += K4_MULT;
				k3 += K4_MULT;
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			LOG.error("Hart_Fast: Failed to factor N=" + N + ". Either it has factors < cbrt(N) needing trial division, or the arrays are too small.");
			// TODO Or N is square
			return 1;
		}
	}
	private long adjustA(long N, long x, int k) {
		if ((k&1)==0)
			return x | 1;
		final long kNp1 = k*N+1;
		if ((kNp1 & 3) == 0)
		{
			return x + ((kNp1 - x) & 7);
		}
		else if ((kNp1 & 7) == 2)
		{
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			final long diff = adjust1<adjust2 ? adjust1 : adjust2;
			return x + diff;
		}
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1<adjust2 ? adjust1 : adjust2);
	}

}
