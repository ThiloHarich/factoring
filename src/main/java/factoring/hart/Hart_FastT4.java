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
 *
 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
 * But the smaller possible factors get, it will become slower and slower.
 *
 * For any kind of test numbers except very hard semiprimes, Hart_TDiv_Race will be faster.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class Hart_FastT4 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Hart_FastT4.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 1; // 315

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
	public Hart_FastT4(boolean doTDivFirst) {
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
		}
		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		int k = K_MULT;
		try {
			for (int i=1; ;i++, k += K_MULT) {
				// odd k -> adjust a mod 8
				final long aOrig = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				a = adjustA(N, aOrig, k);
				//				a = adjustA(N, aOrig, k, 27);
				final long fourKN = k * fourN;
				test = a*a - fourKN;
				//				final long testMod = test % 27;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
						//						final long knMod8 = (k*N) % 8;
						//						if (knMod8 == 1) {
						final long mod = 16;
						//							final int step = (int) ((a-aOrig) % mod );
						final int aMod = (int) (a % mod);
						final int knMod = (int) ((k*N) % 8);
						//							final int diff = (int) ((aMod - knMod) % mod);
						//							System.out.print(aMod + ",");
						//							//						System.out.print((k*N) % 8);
						//						}
						return gcd;
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// this may happen if this implementation is tested with doTDivFirst==false and N having
			// very small factors, or if N is too big
			return 0;
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

	private long adjustA(long N, long x, int k, int mod) {
		final long kNp1 = k*N+1;
		final long adjust1 = Math.floorMod( kNp1 - x, mod);
		final long adjust2 = Math.floorMod(-kNp1 - x, mod);
		return x + (adjust1<adjust2 ? adjust1 : adjust2);
	}
}


