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
package factoring.hart.sieve;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.math.PrimeMath;
import factoring.smooth.SmoothNumbers;

/**
 * A factoring algorithm for numbers up to 2^63 numbers, which generates numbers
 * x^2 - kn. From there numbers the primefactors lower then m are calculated (with gcd).
 * For smooth numbers over m we try to find a product out of these numbers which is a square.
 * This is generate numbers with a modified Hart Algorithm and then combine the solutions 
 * in a Quadratic Sieve manner without sieving.
 * 
 * @author Thilo Harich
 */
public class HartSmooth extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartSmooth.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
//	private static final int K_MULT1 = 315;  // 315
	private static final int K_MULT1 = 1;  // 315
	// product of 11 primes <= 37 
	private static final long PRIMORIAL = 3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l;  
	static int [] primes = {2,3,5,7,11,13,17,19,23,29,31,37};
	
	long[] primeIndex = new long [38];
	
	/** 
	 * Size of arrays: this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52.
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt1;
//	private final double[] sqrt2;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
	 * But the smaller possible factors get, it will become slower and slower.
	 */
	public HartSmooth(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt1 = new double[I_MAX];
//		sqrt2 = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt1[i] = Math.sqrt(i*K_MULT1);
//			if (i%7 != 0) {
//				sqrt2[i] = Math.sqrt(i*K_MULT2);
//			}
		}
		// for all the primes we have an index to a long value which stores the exponent of this
		// prime in the prime factorization. As long as they do not interfere we are fine
		for (int i = 0; i < primes.length; i++) {
			primeIndex[primes[i]] = 1l << (5*i);			
		}
	}

	@Override
	public String getName() {
		return "Hart_Fast2Mult(" + doTDivFirst + ")";
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
		
		// test for exact squares
		final double sqrtN = Math.sqrt(N);
		final long floorSqrtN = (long) sqrtN;
		if (floorSqrtN*floorSqrtN == N) return floorSqrtN;

		final long fourN = N<<2;
		SmoothNumbers smooth = new SmoothNumbers(2,3,5,7,11,13,17,19,23,29,31,37);
		int smoothCount = 0;
		final double sqrt4N = sqrtN*2;
		long a;
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		int k1 = K_MULT1;
		double sum = 0;
		try {
			for (int i=1; smoothCount < 12; i++, k1 += K_MULT1) {
				a = (long) (sqrt4N * sqrt1[i] + ROUND_UP_DOUBLE);
//				long test = a*a - i*fourN;
				long test = a*a - k1 * fourN;
//				double xInv = 1.0 / test;
//				sum += xInv;
				
				final boolean isSmooth = smooth.isSmoothBy3Inv(test, PRIMORIAL);
				smoothCount += isSmooth ? 1 : 0;
//				if (isSmooth) {
//					int[] factors = smooth.factorizeIfSmooth1(test, 37);
//				}
//				final double sqrt = Math.sqrt(test);
//				double rem = sqrt - Math.floor(sqrt);
//				for (int p : primes) {
//					if (rem * p == Math.round(rem * p))
//						System.out.println(p);
//				}
//				long b = (long) sqrt;
//				if (b*b == test) {
//					return b;
//				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			LOG.error("Hart_Fast2Mult: Failed to factor N=" + N + ". Either it has factors < cbrt(N) needing trial division, or the arrays are too small.");
			return 1;
		}
		return PRIMORIAL;
	}
	

	private String factors(long a, long[] primeFacts, int[] composites) {
		String primefactors = "";
		for (int j = 0; j < 12; j++) {
			int exponent = (int) ((primeFacts[0] >> 5*j) & 31);
			if (exponent > 0)
				primefactors += primes[j] + "^" + exponent + " * ";
		}
		final String compString = Arrays.stream(composites).filter(i -> i != 0).mapToObj(Integer::toString).collect(Collectors.joining(" * "));
		return "" + a + " = " +primefactors + compString;
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
		if ((k&1)==0) return x | 1;
		
		final long kNp1 = k*N+1;
		if ((kNp1 & 3) == 0) return x + ((kNp1 - x) & 7);
		
		if ((kNp1 & 7) == 2) {
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			return x + (adjust1 < adjust2 ? adjust1 : adjust2);
		}
		
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1 < adjust2 ? adjust1 : adjust2);
	}

}
