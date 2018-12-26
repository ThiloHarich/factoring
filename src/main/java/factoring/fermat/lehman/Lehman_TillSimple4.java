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
 * Faster implementation of Lehmans factor algorithm .
 * Works for N <= 45 bit.
 *
 * This version also might do trial division before or after (if factorSemiprimes == true) the main loop.
 *
 * @author Tilman Neumann + Thilo Harich
 */
public class Lehman_TillSimple4 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple3.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	static TrialInvFact smallFact= new TrialInvFact((int) (1L << (48/3)));;

	boolean factorSemiprimes;


	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt, sqrtInv;
	long N;
	long fourN;
	double sqrt4N;

	/**
	 *
	 * @param factorSemiprimes if the number to be factored might be a semiprimes were each factor is greater then n^1/3
	 * then set this to true. The algorithm will always work. But this might have a very positive effect on the performance.
	 */
	public Lehman_TillSimple4(boolean factorSemiprimes) {
		this.factorSemiprimes = factorSemiprimes;
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
		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		long factor;
		smallFact.setMaxFactor(cbrt);
		if (!factorSemiprimes && (factor = smallFact.findFactors(N, null)) != N)
			return factor;

		// limit for must be 0 mod 6, since we also want to search above of it
		final int kLimit = (cbrt + 12) / 12 * 12;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoA = (cbrt >> 6);
		// twoA = 0 mod 6
		kTwoA = ((kTwoA + 12)/ 12) * 12;
		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability

		// first investigate in solutions a^2 - sqrt(k*n) = y^2 were we only have two possible 'a' values per k
		// but do not go to far, since there we have a lower chance to finda factor
		// here we only inspect k = 6*i since they much more likely have a solution x^2 - sqrt(k*n) = y^2
		// this is we use a multiplier of 6
		factor = lehmanAOdd(kTwoA + 0, kLimit);
		if (factor != N && factor > 1)
			return factor;



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



		// continue even loop, because we are looking at very high numbers this now done after the k = 3 mod 6 loop
		if ((factor = lehmanAOdd(kLimit, kLimit << 1)) > 1)
			return factor;

		// next best candidates
		if ((factor = lehmanAknMod4(kTwoA + 9, kLimit)) > 1 && factor != N)
			return factor;
		if ((factor = lehmanAknMod4(kTwoA + 3, kLimit)) > 1 && factor != N)
			return factor;


		// Check via trial division whether N has a nontrivial divisor d <= cbrt(N).
		return smallFact.findFactors(N, null);
	}

	/**
	 * This method tries to choose an a such that n + k = a mod 4.
	 * This can be applied for odd kBegin.
	 * If a*a - 4 * k * N is a square it returns true.
	 * @param kBegin
	 * @param kLimit
	 * @return
	 */
	long lehmanAknMod4(int kBegin, final int kLimit) {
		for (int k = kBegin; k <= kLimit; k += 12) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			// for k = 0 mod 6 a must be even and k + n + a = 0 mod 4
			a += (k + N - a) & 3;
			//			a = (a & 1) == 1 ? a+1 : a;
			final long test = a*a - k * fourN;
			{
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}
		return -1;
	}

	/**
	 * makes a odd can be applied to even kBegin.
	 * @param kBegin
	 * @param kEnd
	 * @return
	 */
	long lehmanAOdd(int kBegin, final int kEnd) {
		for (int k = kBegin ; k <= kEnd; k += 6) {
			// for k = 0 mod 6 a must be odd
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		return -1;
	}
}