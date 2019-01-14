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

import java.util.Collection;

import org.apache.log4j.Logger;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 * Faster implementation of Lehman's factor algorithm.
 * Works flawlessly for N <= 56 bit.
 *
 * This version will do trial division before or after the main loop depending on factorSemiprimes.
 *
 * @authors Tilman Neumann + Thilo Harich
 */
public class LehmanMod30 implements FactorizationOfLongs, FactorFinder {
	private static final Logger LOG = Logger.getLogger(Lehman_Fast30.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static double[] sqrt, sqrtInv, sqrt6, sqrt30;

	private static final int kMax = 1<<22;


	static {
		// Precompute sqrts for all possible k. 2^22 entries are enough for N~2^66.
		final long start = System.currentTimeMillis();
		sqrt = new double[kMax + 1];
		sqrt6 = new double[kMax/6 + 1];
		sqrt30 = new double[kMax/30 + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			if (i < sqrt6.length)
				sqrt6[i] = Math.sqrt(6*i);
			if (i < sqrt30.length)
				sqrt30[i] = Math.sqrt(30*i);
			sqrtInv[i] = 1.0/sqrtI;
		}
		final long end = System.currentTimeMillis();
		System.out.println("Time Init sqrt   : + " + (end - start));
	}

	long N;
	long fourN;
	double sqrt4N;
	boolean factorSemiprimes;

	FactorizationOfLongs smallFactorizer = new TrialInvFact(kMax);

	private long twentyfourN;
	private long N120;

	/**
	 * Only constructor.
	 * @param factorSemiprimes if the number to be factored might be a semiprimes were each factor is greater then n^1/3
	 * then set this to true. The algorithm will always work. But this might have a very positive effect on the performance.
	 */
	public LehmanMod30(boolean factorSemiprimes) {
		this.factorSemiprimes = factorSemiprimes;
	}


	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		N = n;
		int cbrt = (int) Math.cbrt(N);

		smallFactorizer.setMaxFactor(cbrt);
		// factor out all small factors if
		if (!factorSemiprimes) {
			final long nAfterTrial = smallFactorizer.findFactors(n, primeFactors);
			if (nAfterTrial == 1)
				return nAfterTrial;
			N = nAfterTrial;
		}
		// re-adjust the maximal factor we have to search for. If small factors were found by trial division,
		// which is quite often the case for arbitrary numbers, this cuts down the runtime dramatically.
		// n^1/3 is the upper limit of the outer loop of the original Lehman algorithm
		cbrt = (int) Math.cbrt(N);

		// In the first phase the inspect only k as multiples of 30
		final int bigStep = 30;
		final int bigStepDivSmallStep = 5;
		// limit divided by 30 -> limit for the big steps
		final int kLimit30 = (cbrt + bigStep) / bigStep;
		// limit for the small steps
		final int kLimit6 = kLimit30 * bigStepDivSmallStep;
		final int kLimit = kLimit30 * 30;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoA = (cbrt >> 6);
		//  twoA = 0 mod 30
		final int kTwoA30 = (kTwoA + bigStep)/ bigStep;
		final int kTwoA6 = kTwoA30*bigStepDivSmallStep;
		kTwoA = kTwoA30 * 30;
		fourN = N<<2;
		twentyfourN = N* 24;
		N120 = N * 120;
		sqrt4N = Math.sqrt(fourN);
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability

		// first investigate in solutions a^2 - sqrt(k*n) = y^2 were we only have two possible 'a' values per k
		// but do not go to far, since there we have a lower chance to find a factor
		// here we only inspect k = 30*i since they much more likely have a solution a^2 - sqrt(k*n) = y^2
		// this is we use a multiplier of 30 for odd 'a' we have much higher chances to find solutions
		long factor = lehman30(kTwoA30, kLimit30);
		if (factor>1 && factor != N)
			return factor;

		// inspect k = 6*i != 30*j
		factor = lehman6(kTwoA6 + 1, kLimit6);
		if (factor>1 && factor != N)
			return factor;

		// investigate in solution a^2 - sqrt(k*n) = y^2 were we might have more then two solutions 'a'
		for (int k=1 ; k < kTwoA30 * 30; k++) {
			final double sqrt4kN = sqrt4N * sqrt[k];
			// only use long values
			long aStart = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			final long aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInv[k]);
			long aStep;
			if ((k & 1) == 0) {
				// k even -> make sure aLimit is odd
				aStart |= 1l;
				aStep = 2;
			} else {
				final long kPlusN = k + N;
				if ((kPlusN & 3) == 0) {
					aStep = 8;
					aStart += (kPlusN - aStart) & 7;
				} else {
					aStep = 4;
					aStart += ((kPlusN - aStart) & 3);
				}
			}
			for (long a=aStart; a <= aLimit + (aStep >> 1); a+=aStep) {
				final long test = a*a - k * fourN;
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return PrimeMath.gcd(a+b, N);
				}
			}
		}
		// So far we have only odd solutions for 'a' now we try to get all the even solutions
		factor = lehmanOdd(kTwoA + 3, kLimit);
		if (factor > 1)
			return factor;

		// continue even k loop. Theoretically we might have missed solutions for k = 1,2,4,5 mod 6.
		// but when looking up to 6 * kLimit it seems to be we are not missing one of them.
		// to be sure we might just increase the limit even further

		factor = lehman30(kLimit30, 6*kLimit30);
		if (factor>1) return factor;

		factor = lehman6(kLimit6 + 1, kLimit);
		if (factor > 1)
			return factor;


		// Check via trial division whether N has a nontrivial divisor d <= cbrt(N).
		//LOG.debug("entering tdiv...");
		factor = smallFactorizer.findFactors(n, primeFactors);
		if (factor>1 && factor<N) return factor;

		for (int k=kTwoA6 + 1; k <= kLimit; k++) {
			final long fourKN = k*N<<2;
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - fourKN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				final long gcd = PrimeMath.gcd(a+b, N);
				if (gcd>1 && gcd<N) {
					return gcd;
				}
			}
		}

		return 0; // fail
	}

	private long lehmanOdd(int kBegin, final int kLimit) {
		for (int k = kBegin; k <= kLimit; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			// for k odd a must be even and k + n + a = 0 mod 4
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else
			{
				a += ((kPlusN - a) & 3);
			}
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return PrimeMath.gcd(a+b, N);
			}
		}

		return -1;
	}

	private long lehman6(int kBeginIndex, final int kEnd) {
		for (int k = kBeginIndex ; k <= kEnd;) {
			long a  = (long) (sqrt4N * sqrt6[k]   + ROUND_UP_DOUBLE) | 1;
			long k24N = k++ * twentyfourN;
			long test = a*a - k24N;
			long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return PrimeMath.gcd(a+b, N);
			}
			k24N += twentyfourN;
			a = (long) (sqrt4N * sqrt6[k++] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return PrimeMath.gcd(a+b, N);
			}
			k24N += twentyfourN;
			a = (long) (sqrt4N * sqrt6[k++] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return PrimeMath.gcd(a+b, N);
			}
			k24N += twentyfourN;
			a = (long) (sqrt4N * sqrt6[k] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return PrimeMath.gcd(a+b, N);
			}
			// jump over k = 0 mod 5 / 30
			k += 2;
		}
		return -1;
	}

	private long lehman30(int kBeginIndex, final int kEnd) {
		for (int k = kBeginIndex ; k <= kEnd; k++) {
			//			long a  = (long) (sqrt4N * sqrt30[k] + ROUND_UP_DOUBLE);
			final long a  = (long) (sqrt4N * sqrt30[k] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k * N120;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				//				System.out.print("|" + a % 30);
				return PrimeMath.gcd(a+b, N);
			}
		}
		return -1;
	}
}
