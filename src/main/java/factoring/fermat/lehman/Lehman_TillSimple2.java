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
package factoring.fermat.lehman;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * Faster implementation of Lehmans factor algorithm following https://programmingpraxis.com/2017/08/22/lehmans-factoring-algorithm/.
 * Many improvements inspired by Thilo Harich (https://github.com/ThiloHarich/factoring.git).
 * Works for N <= 45 bit.
 *
 * This version does trial division after the main loop.
 *
 * @author Tilman Neumann
 */
public class Lehman_TillSimple2 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple2.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final boolean[] isSquareMod1024 = isSquareMod1024();

	private static boolean[] isSquareMod1024() {
		final boolean[] isSquareMod_1024 = new boolean[1024];
		for (int i = 0; i < 1024; i++) {
			isSquareMod_1024[(i * i) & 1023] = true;
		}
		return isSquareMod_1024;
	}

	/** A multiplicative constant to adjust the trial division limit. */
	private final float tDivLimitMultiplier;

	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt, sqrtInv;

	public Lehman_TillSimple2(float tDivLimitMultiplier) {
		this.tDivLimitMultiplier = tDivLimitMultiplier;
		SMALL_PRIMES.ensurePrimeCount(10000); // for tDivLimitMultiplier ~ 2 we need more than 4793 primes
		initSqrts();
	}

	private void initSqrts() {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = (int) (tDivLimitMultiplier*Math.cbrt(1L<<45) + 1);
		//LOG.debug("kMax = " + kMax);

		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0/sqrtI;
		}
		LOG.info("Lehman: Built sqrt tables for multiplier " + tDivLimitMultiplier + " with " + sqrt.length + " entries");
	}

	@Override
	public String getName() {
		return "Lehman_TDivLast(" + tDivLimitMultiplier + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		final double cbrt = Math.ceil(Math.cbrt(N));

		// 1. Main loop for small k, where we can have more than four a-value

		final int cbrtInt = (int) (cbrt);
		final int kLimit = (cbrtInt >> 2) + 1;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		// make it odd
		final int fourA = (kLimit >> 4 - 1) | 1;
		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		final double sixthRootTerm = Math.pow(N, 1/6.0); // double precision is required for stability
		int k=1;
		for (; k <= fourA; k++) {
			final double sqrt4kN = sqrt4N * sqrt[k];
			// only use long values
			long aStart = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			final long aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInv[k]);
			long aStep;
			if ((k & 1) == 0) {
				// k even -> make sure aLimit is odd
				//				aLimit |= 1l;
				aStart |= 1l;
				aStep = 2;
			} else {
				final long kn = k*N;
				// this extra case gives ~ 5 %
				if ((kn & 3) == 3) {
					aStep = 8;
					//					aLimit += ((7 - kn - aLimit) & 7);
					aStart += ((7 - kn - aLimit) & 7);
				} else
				{
					aStep = 4;
					//					aLimit += ((k + N - aLimit) & 3);
					aStart += ((k + N - aLimit) & 3);
				}
			}

			// processing the a-loop top-down is faster than bottom-up
			//			for (long a=aLimit; a >= aStart; a-=aStep) {
			for (long a=aStart; a <= aLimit; a += aStep) {
				//				final long test = a*a - (kn << 2);
				final long test = a*a - k * fourN;
				if (isSquareMod1024[(int) (test & 1023)]) {
					final long b = (long) Math.sqrt(test);
					if (b*b == test) {
						final long gcd = gcdEngine.gcd(a+b, N);
						return gcd;
					}
				}
			}
		}
		// 2. continue main loop for larger even k, where we can have only 4 a values per k
		for (int k1 = k; k1 <= kLimit; k1+=2) {
			final long a = (long) (sqrt4N * sqrt[k1] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k1 * fourN;
			if (isSquareMod1024[(int) (test & 1023)]) {
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		// 3. continue main loop for larger odd k
		int k2 = k + 1;
		for ( ; k2 <= kLimit; k2 += 2) {
			long a = (long) (sqrt4N * sqrt[k2] + ROUND_UP_DOUBLE);
			a += (k2 + N - a) & 3;
			final long test = a*a - k2 * fourN;
			if (isSquareMod1024[(int) (test & 1023)]) {
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		// 4. Check via trial division whether N has a nontrivial divisor d <= cbrt(N), and if so, return d.
		final int tDivLimit = (int) (tDivLimitMultiplier*cbrt);
		int i=0, p;
		while ((p = SMALL_PRIMES.getPrime(i++)) <= tDivLimit) {
			if (N%p==0) return p;
		}

		// Nothing found. Either N is prime or the algorithm didn't work because N > 45 bit.
		return 0;
	}
}
