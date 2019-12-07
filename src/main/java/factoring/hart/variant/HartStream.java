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
package factoring.hart.variant;

import java.math.BigInteger;
import java.util.OptionalInt;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * A hart lehman factoring algorithm which uses a Lambda approach.
 * It has half the speed of a iterative implementation.
 * Since we precalculate the sqrt of the k we can use a stream of the
 * sqrt's.
 *
 * @authors Thilo Harich
 */
public class HartStream extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartStream.class);

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

	private final boolean doTDivFirst;
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public HartStream(boolean doTDivFirst) {
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

	public boolean hasSquare (double sqrt, double sqrt4n, long n, int k) {
		long a = (long) (sqrt4n * sqrt + ROUND_UP_DOUBLE);
		a = adjustA(n, a, k);
		final long fourN = n<<2;
		final long test = a*a - k * fourN;
		final long b = (long) Math.sqrt(test);
		if (b*b == test) {
			long gcd;
			if ((gcd = gcdEngine.gcd(a+b, n))>1 && gcd<n) {
				return true;
			}
		}
		return false;
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


		final OptionalInt goodIndex = IntStream.range(0, sqrt.length) // .parallel().unordered() is ~ 20 times slower then the sequential stream
				.filter(i -> hasSquare(sqrt[i], sqrt4N, N, i*K_MULT)) // .mapToLong(i -> getFactor(sqrt[i], sqrt4N, N, i*K_MULT)) is 50% slower
				.findAny(); // .findFirst() makes no difference
		if (goodIndex.isPresent()) {
			final int i = goodIndex.getAsInt();
			long a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
			final int k = i*K_MULT;
			a = adjustA(N, a, k);
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			final long gcd = gcdEngine.gcd(a+b, N);
			return gcd;
		}

		//		final OptionalLong factor = IntStream.range(0, sqrt.length)
		//				.filter(i -> hasSquare(sqrt[i], sqrt4N, N, i*K_MULT)).mapToLong(i -> getFactor(sqrt[i], sqrt4N, N, i*K_MULT)).findAny();
		//
		//		if (factor.isPresent())
		//			return factor.getAsLong();
		return 1;
	}
	private long getFactor(double sqrt, double sqrt4n, long n, int k) {
		long a = (long) (sqrt4n * sqrt + ROUND_UP_DOUBLE);
		a = adjustA(n, a, k);
		final long fourN = n<<2;
		final long test = a*a - k * fourN;
		final long b = (long) Math.sqrt(test);
		final long gcd = gcdEngine.gcd(a+b, n);
		return gcd;
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
