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

import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 *
 * Works for N <= 45 bit.
 *
 * This version does trial division after the main loop.
 *
 * @author Tilman Neumann + Thilo Harich
 */
public class LehmanSimple2 implements FactorizationOfLongs, FactorFinder {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final Gcd63 gcdEngine = new Gcd63();
	final FactorFinder smallFactoriser;
	boolean hardNumbers;

	private double[] sqrt;

	public LehmanSimple2(boolean hardNumbers) {
		this.hardNumbers = hardNumbers;
		final int cubicRoot = (int) Math.cbrt(1L<<50);
		initSqrts(cubicRoot);
		smallFactoriser = new TrialInvFact(cubicRoot);
	}

	private void initSqrts(int cubicRoot) {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = 40 * cubicRoot;
		sqrt = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			sqrt[i] = Math.sqrt(i);
		}
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		int kLimit = (int)  Math.ceil(Math.cbrt(n));
		if (!hardNumbers) {
			smallFactoriser.setMaxFactor(kLimit);
			final long factor = smallFactoriser.findFactors(n, primeFactors);
			if (factor == 1) return factor;
			n = factor;
		}
		kLimit = (int)  Math.ceil(Math.cbrt(n));
		final long fourN = n<<2;
		final double sqrt4N = Math.sqrt(fourN);
		//		final long MN = fourN * 2*3*3*5*7;// = 630
		final long MN = fourN * 3*3*5*7;// = 315
		final double sqrtMN = Math.sqrt(MN);

		long aForK1 = (long) (sqrt4N  + ROUND_UP_DOUBLE);
		long aStep;
		if ((n & 3) == 3) {
			aStep = 8;
			aForK1 += ((7 - n - aForK1) & 7);
		} else
		{
			aStep = 4;
			aForK1 += ((1 + n - aForK1) & 3);
		}
		long gcd = 1;
		// we have to increase kLimit here, because we only take good but seldom k's
		// TODO what is the max of the loop
		for (int k=1; k < 40*kLimit; aForK1 += aStep) {
			for (int i=0; i < 8; i++) {
				// investigate in k = 0 mod multiplier only like the Hart variant does it
				long  a, test, b;
				// first check odd values of k
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE);
				final long kPlusN = k + n;
				if ((kPlusN*3 & 3) == 0) // 315 == 3 mod 8
				{
					a += ((kPlusN - a) & 7);
				} else
				{
					a += ((kPlusN - a) & 3);
				}
				// even values of k
				test = a*a - k++ * MN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = PrimeMath.gcd(a+b, n);
				}
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE) | 1;
				test = a*a - k++ * MN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = PrimeMath.gcd(a+b, n);
				}
			}
			// Here k is always 1 and we increase 'a' by aStep.
			// This phase should save us if the first phase does not find anything
			//			final long test = aForK1*aForK1 - fourN;
			//			final long b = (long) Math.sqrt(test);
			//			if (b*b == test) {
			//				gcd= gcdEngine.gcd(aForK1+b, n);
			//			}
			//			if (gcd > 1 && gcd < n)
			//				return gcd;
		}

		smallFactoriser.setMaxFactor(kLimit);
		final long factor = smallFactoriser.findFactors(n, primeFactors);
		if (factor == 1)
			return factor;

		// If sqrt(4kN) is very near to an exact integer then the fast ceil() in the 'aStart'-computation
		// may have failed. Then we need a "correction loop":
		for (int k=1; k <= kLimit; k++) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - k*fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, n);
			}
		}
		return -1;
	}
}
