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
import factoring.trial.TrialInvFact;

/**
 * Faster implementation of Lehmans factor algorithm following https://programmingpraxis.com/2017/08/22/lehmans-factoring-algorithm/.
 * Many improvements inspired by Thilo Harich (https://github.com/ThiloHarich/factoring.git).
 * Works for N <= 45 bit.
 *
 * This version does trial division after the main loop.
 *
 * @author Tilman Neumann
 */
public class Lehman_TillSimple3 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple3.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	TrialInvFact smallFact;


	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt, sqrtInv;

	public Lehman_TillSimple3() {
		//		SMALL_PRIMES.ensurePrimeCount(10000); // for tDivLimitMultiplier ~ 2 we need more than 4793 primes
		smallFact = new TrialInvFact((int) (1L << 48/3));
		initSqrts();
	}

	private void initSqrts() {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = (int) (Math.cbrt(1L<<48) + 1);
		//LOG.debug("kMax = " + kMax);

		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0/sqrtI;
		}
	}

	@Override
	public String getName() {
		return "Lehman_TDivLast()";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		final double cbrt = Math.ceil(Math.cbrt(N));

		// 1. Main loop for small k, where we can have more than four a-value

		final int kLimit = ((int) cbrt  + 6) / 6 * 6;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoA = Math.max(((kLimit >> 6) - 1), 0) | 1;
		// twoA = 0 mod 6
		kTwoA = ((kTwoA + 6)/ 6) * 6;
		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability

		// investigate in solutions a^2 - sqrt(k*n) = y^2 were we only have two possible 'a' values per k
		// but do not go to far, since they might have a lower chance to produce a solution
		// here we only inspect k = 6*i since they much more likely have a solution x^2 - sqrt(k*n) = y^2
		for (int k = kTwoA ; k <= kLimit; k += 6) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				final long gcd = gcdEngine.gcd(a+b, N);
				if (gcd != N)
					return gcd;
			}
		}

		// investigate in solution a^2 - sqrt(k*n) = y^2 were we might have more then two solutions 'a'
		for (int k=1 ; k < kTwoA; k++) {
			final double sqrt4kN = sqrt4N * sqrt[k];
			// only use long values
			final long aStart = (long) (sqrt4kN + ROUND_UP_DOUBLE); // much faster than ceil() !
			long aLimit = (long) (sqrt4kN + sixthRootTerm * sqrtInv[k]);
			long aStep;
			if ((k & 1) == 0) {
				// k even -> make sure aLimit is odd
				aLimit |= 1l;
				aStep = 2;
			} else {
				//				final long kn = kSmall*N;
				// this extra case gives ~ 5 %
				final long kPlusN = k + N;
				if ((kPlusN & 3) == 0) {
					aStep = 8;
					aLimit += ((kPlusN - aLimit) & 7);
				} else
				{
					aStep = 4;
					aLimit += ((kPlusN - aLimit) & 3);
				}
			}

			// processing the a-loop top-down is faster than bottom-up
			for (long a=aLimit; a >= aStart; a-=aStep) {
				final long test = a*a - k * fourN;
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		// additional loop for k = 3 mod 6 in the middle range
		// this loop has fewer possible a's, and therefore gives less often factors
		for (int k = kTwoA + 3; k <= kLimit; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			a += (k + N - a) & 3;
			final long test = a*a - k * fourN;
			{
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		// switch back to the main loop for bigger k
		for (int k = kLimit ; k <= kLimit << 1; k += 6) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}

		//	continue additional loop for larger k = 3 mod 6.
		for (int k = kLimit + 3; k <= kLimit << 1; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			a += (k + N - a) & 3;
			final long test = a*a - k * fourN;
			{
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}


		// 4. Check via trial division whether N has a nontrivial divisor d <= cbrt(N), and if so, return d.
		return smallFact.findFactors(N, null);
		//		final int tDivLimit = (int) (tDivLimitMultiplier*cbrt);
		//		int i=0, p;
		//		while ((p = SMALL_PRIMES.getPrime(i++)) <= tDivLimit) {
		//			if (N%p==0) {
		//				return p;
		//			}
		//		}
	}
}
