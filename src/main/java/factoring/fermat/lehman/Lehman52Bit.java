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
public class Lehman52Bit extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman52Bit.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final BigInteger big24 = BigInteger.valueOf(24);

	static TrialInvFact smallFact= new TrialInvFact((int) (1L << (48/3)));;

	boolean factorSemiprimes;


	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt, sqrtInv;
	long N;
	long fourN;
	BigInteger fourNBig = BigInteger.valueOf(fourN);
	double sqrt4N;

	/**
	 *
	 * @param factorSemiprimes if the number to be factored might be a semiprimes were each factor is greater then n^1/3
	 * then set this to true. The algorithm will always work. But this might have a very positive effect on the performance.
	 */
	public Lehman52Bit(boolean factorSemiprimes) {
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
		final double cbrt = Math.ceil(Math.cbrt(N));

		long factor;
		smallFact.setMaxFactor((int) cbrt);
		if (!factorSemiprimes && (factor = smallFact.findFactors(N, null)) != N)
			return factor;

		final int kLimit = ((int) cbrt  + 6) / 6 * 6;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoA = Math.max(((kLimit >> 6) - 1), 0) | 1;
		// twoA = 0 mod 6
		kTwoA = ((kTwoA + 6)/ 6) * 6;
		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability

		// first investigate in solutions a^2 - sqrt(k*n) = y^2 were we only have two possible 'a' values per k
		// but do not go to far, since there we have a lower chance to finda factor
		// here we only inspect k = 6*i since they much more likely have a solution x^2 - sqrt(k*n) = y^2
		// this is we use a multiplier of 6
		factor = lehmanEven(kTwoA, kLimit);
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

		// additional loop for k = 3 mod 6 in the middle range
		// this loop has fewer possible a's, then k = 0 mod 6 and therefore gives less often factors
		if ((factor = lehmanOdd(kTwoA + 3, kLimit)) > 1)
			return factor;

		// continue even loop, because we are looking at very high numbers this now done after the k = 3 mod 6 loop
		if ((factor = lehmanEven(kLimit, kLimit << 1)) > 1)
			return factor;

		//	continue additional loop for larger k = 3 mod 6.
		if ((factor = lehmanOdd(kLimit + 3, kLimit << 1)) > 1)
			return factor;

		// the following code might not be reached
		if ((factor = lehmanEven(kTwoA + 2, kLimit)) > 1)
			return factor;
		if ((factor = lehmanEven(kTwoA + 4, kLimit)) > 1)
			return factor;

		// Check via trial division whether N has a nontrivial divisor d <= cbrt(N).
		return smallFact.findFactors(N, null);
	}

	long lehmanOdd(int kBegin, final int kLimit) {
		BigInteger fourKNBig = BigInteger.valueOf(kBegin).multiply(fourNBig);
		//		BigInteger twoPow26 = BigInteger.valueOf(1l << 26);
		//		BigInteger[] divAndRem = fourKNBig.divideAndRemainder(twoPow26);
		//		long[] fourKNs = {divAndRem[0].longValue(), divAndRem[1].longValue()};

		for (int k = kBegin; k <= kLimit; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			// for k = 0 mod 6 a must be even and k + n + a = 0 mod 4
			a += (k + N - a) & 3;
			//			long[] aSquare = BigDouble.square(a);
			//			aBig = BigInteger.valueOf(	)
			final long test = a*a - k * fourN;
			{
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
			fourKNBig = fourKNBig.add(fourKNBig.multiply(big24));
		}
		return -1;
	}

	long lehmanEven(int kBegin, final int kEnd) {
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
