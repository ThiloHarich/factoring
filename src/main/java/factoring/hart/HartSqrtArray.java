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
 * This implementations introduces some improvements that make it the fastest factoring algorithm
 * for numbers with more then 20? and less then 50 bit.
 * It avoids the complexity of calculating the square root when factoring multiple numbers,
 * by storing the square roots of k in the main loop.
 * It uses an optimized trial division algorithm to factorize small numbers.
 * It uses a well chosen multiplier m = 3*3*5*7 which is odd.
 * After calculating a number 'a' above sqrt(4*m*k) a will be adjusted to satisfy
 * some modulus a power of 2 argument.
 * It reuses the idea of rounding up by adding a well choosen constant (Warren D. Smith)
 *
 * It tires to find solutions for a^2 - 4*m*i*n = b^2 from fermat we then know that
 * gcd(a+b, n) and gcd(a-b, n) are divisors of n.
 *
 * This is done by one simple loop over k were we generate numbers a = sqrt(4*m*k*n).
 * By storing sqrt(k) in an array this can be calculated fast.
 *
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
public class HartSqrtArray extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartSqrtArray.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315

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

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public HartSqrtArray(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt[i] = Math.sqrt(i * (double)K_MULT);
		}
	}

	@Override
	public String getName() {
		return "Hart_FastT(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger n) {
		return BigInteger.valueOf(findSingleFactor(n.longValue()));
	}

	/**
	 * Find a factor of long N.
	 * @param n
	 * @return factor of N
	 */
	public long findSingleFactor(long n) {
		if (doTDivFirst) {
			// do trial division before the Hart loop
			tdiv.setTestLimit((int) Math.cbrt(n));
			final long factor = tdiv.findSingleFactor(n);
			if (factor > 1) return factor;
		}

		final long fourN = n<<2;
		final double sqrt4N = Math.sqrt(fourN);
		int k = K_MULT;
		for (int i=1; ;i++, k += K_MULT) {
			// calculating the sqrt here is 5 times slower then storing it
			long a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
			a = adjustA(n, a, k);
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			long gcd;
			if (b*b == test && (gcd = gcdEngine.gcd(a+b, n))>1 && gcd<n) {
				return gcd;
			}
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
	 * @param n
	 * @param x
	 * @param k
	 * @return
	 */
	private long adjustA(long n, long x, int k) {
		if ((k&1)==0)
			return x | 1;
		final long kNp1 = k*n+1;
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
