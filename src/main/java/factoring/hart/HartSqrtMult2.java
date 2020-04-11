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

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.SortedMultiset;
import factoring.smooth.SmoothNumbers;
import org.apache.log4j.Logger;

import java.math.BigInteger;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Instead of looking for a solution x^2 - kn = y^2 we look for
 * x^2 - kn = s * y^2, where s is smooth (and free of squares).
 * If we have much more values of this kind as real squares, we might combine
 * such relations (or the different s) to a a square.
 * But unfortunately (i expected more such numbers) it happens to seldom to get
 * a benefit out of it.
 *
 * We had to calculate the s, which takes more time then calculating just the sqrt.
 *
 * Here s is chosen as a product of the first primes which fit in a long value,
 * such that we can determine by some squaring mod s and a final gcd.
 *
 *
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class HartSqrtMult2 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartSqrtMult2.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
//	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 3*7*17; // 315
	//	private static final int K_MULT = 5*11*13; // 315
	//	private static final int K_MULT = 45;
		private static final int K_MULT = 1; // 315

	/** Size of arrays this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt;
//	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final TDiv63Inverse tdiv = new TDiv63Inverse(31);
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public HartSqrtMult2(boolean doTDivFirst) {
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
//			tdiv.setTestLimit((int) Math.cbrt(N));
			tdiv.setTestLimit(31);
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
//				a = adjustA(N, a, k);
				test = a*a - k * fourN;
				int shifts = Long.numberOfTrailingZeros(test);
				long smooth = SmoothNumbers.remainder(test, SmoothNumbers.PRIMORIAL);
				SortedMultiset<BigInteger> factors = tdiv.factor(BigInteger.valueOf(test));
				long nonSquare = 1;
				long test2 = test;
				for (Map.Entry<BigInteger, Integer> entry: factors.entrySet()){
					if (entry.getValue() % 2 == 1 && entry.getKey().longValue() <= 31){
						nonSquare *= entry.getKey().longValue();
						test2 /= entry.getKey().longValue();
					}
				}
//				assertEquals(test/nonSquare, test2);
				if (test/nonSquare != test2){
					System.out.println();
				}
//				long test2 = test >> shifts;
				b = (long) Math.sqrt(test);
				long b2 = (long) Math.sqrt(test2);
				if (b2*b2 == test2) {
					System.out.print(k+ ",");
				}
				if (b*b == test) {
					System.out.println(k + "XXX");
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
						return gcd;
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			LOG.error("Hart_Fast: Failed to factor N=" + N + ". Either it has factors < cbrt(N) needing trial division, or the arrays are too small.");
			// TODO Or N is square
			return 1;
		}
	}
	/**
	 * Increases x to return the next possible solution for x for x^2 - 4kn = b^2.
	 * Due to performance reasons we give back solutions for this equations modulo a
	 * power of 2, since we can determine the solutions just by additions and binary
	 * operations.
	 *
	 * if k is even x must be odd.
	 * if k*n == 3 mod 4 -> x = k*n+1 mod 8
	 * if k*n == 1 mod 8 -> x = k*n+1 mod 16 or -k*n+1 mod 16
	 * if k*n == 5 mod 8 -> x = k*n+1 mod 32 or -k*n+1 mod 32
	 *
	 * @param N
	 * @param x
	 * @param k
	 * @return
	 */
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
			return x + diff;
		}
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1<adjust2 ? adjust1 : adjust2);
	}

}
