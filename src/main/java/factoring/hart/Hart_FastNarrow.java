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
 * Here we try to not only use the smallest a' > a which fulfills a mod 32 Argument. We also
 * use the next a'' > a'
 * @authors Thilo Harich & Tilman Neumann
 */
public class Hart_FastNarrow extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_FastNarrow.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 1155;
	//	private static final int K_MULT = 1; // 315

	/** Size of arrays this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final long MAX_STEP = 2;

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Hart_FastNarrow(boolean doTDivFirst) {
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
		long a, b, test, gcd;
		int k = K_MULT;
		try {
			for(int range = 1; ; range *= 2)
			{
				for (int i=range; i < range * 2;i++, k+= K_MULT) {
					//					k = i*K_MULT;
					// calculating the sqrt here is 5 times slower then storing it
					a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
					//					a = adjustA(N, a, k, true);
					a = adjustA(N, a, k);
					if (a > 0) {
						test = a*a - k * fourN;
						b = (long) Math.sqrt(test);
						if (b*b == test) {
							if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
								return gcd;
							}
						}
					}
				}
				//				for (int i=range; i < range * 2;i+=2) {
				//					k = i*K_MULT;
				//					// calculating the sqrt here is 5 times slower then storing it
				//					a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				//					a = adjustA(N, a, k, false);
				//					if (a > 0) {
				//						test = a*a - k * fourN;
				//						b = (long) Math.sqrt(test);
				//						if (b*b == test) {
				//							if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
				//								return gcd;
				//							}
				//						}
				//					}
				//				}
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
		else if ((kNp1 & 7) == 2) {
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			final long diff = adjust1<adjust2 ? adjust1 : adjust2;
			//					final int xMod = (int) ((x) % 16);
			//					final int resMod = (int) ((x+ diff) % 16);
			return x + diff;
		}
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1<adjust2 ? adjust1 : adjust2);
	}

	private long adjustA(long N, long x, int k, boolean narrow) {
		if ((k&1)==0)
			return narrow  ? x | 1 : (x | 1) +2;
		final long kNp1 = k*N+1;
		if ((kNp1 & 3) == 0)
		{
			final long step = (kNp1 - x) & 7;
			return narrow ? x + step : x + step + 8;
		}
		else if ((kNp1 & 7) == 2) {
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			if (adjust1<adjust2)
				return narrow ? x + adjust1 : x + adjust2;
			else
				return narrow ? x + adjust2 : x + adjust1;
		}
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		if (adjust1<adjust2)
			return narrow ? x + adjust1 : x + adjust2;
		else
			return narrow ? x + adjust2 : x + adjust1;
	}

}
