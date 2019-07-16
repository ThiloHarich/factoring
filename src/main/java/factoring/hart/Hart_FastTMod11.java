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
public class Hart_FastTMod11 extends FactorAlgorithm {
	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 1; // 315

	/** Size of arrays */
	private static final int I_MAX = 1<<18;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();
	private boolean[] sqares11;

	private final double inv11 = 1/ 11.0;

	public static int [] ks = new int[I_MAX];

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public Hart_FastTMod11(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		// sqrt ((1 - 1 * 2^-12)) /13)
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
		else {
			// at least do trial division by multiplier factors
			if ((N%3)==0) return 3;
			if ((N%5)==0) return 5;
			if ((N%7)==0) return 7;
		}

		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		int k = K_MULT;
		try {
			for (int i=1; ;i++, k += K_MULT) {
				// odd k -> adjust a mod 8
				// calculating the sqrt here is 5 times slower then storing it
				a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
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
				test = a*a - k * fourN;
				final int testMod11 = (int) (test - (long) (test*inv11)*11);
				if (sqares11[testMod11]) {
					b = (long) Math.sqrt(test);
					if (b*b == test) {
						if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
							//						System.out.print("," + i );
							//						System.out.print("," + i % 9 + ";" + a % 9);
							return gcd;
						}
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// this may happen if this implementation is tested with doTDivFirst==false and N having
			// very small factors, or if N is too big
			return 0;
		}
	}


}
