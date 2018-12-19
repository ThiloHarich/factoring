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
public class Lehman_TillSuperSimple extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;


	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt;

	public Lehman_TillSuperSimple() {
		initSqrts();
	}

	private void initSqrts() {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = (int) (6*Math.cbrt(1L<<45) + 1);
		//LOG.debug("kMax = " + kMax);

		sqrt = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
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
		final int kLimit = (int)  Math.ceil(Math.cbrt(N));;

		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);

		long aForK1 = (long) (sqrt4N  + ROUND_UP_DOUBLE);
		long aStep;
		if ((N & 3) == 3) {
			aStep = 8;
			aForK1 += ((7 - N - aForK1) & 7);
		} else
		{
			aStep = 4;
			aForK1 += ((1 + N - aForK1) & 3);
		}


		// we have to increase kLimit here, because for each k we only take one a
		for (int k=2; k < 6*kLimit; k+=2) {
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			long test = a*a - k * fourN;
			long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
			// Here k is always 1 and we increase a by aStep.
			test = aForK1*aForK1 - fourN;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(aForK1+b, N);
			}
			aForK1 += aStep;
		}
		//		 4. Check via trial division whether N has a nontrivial divisor d <= cbrt(N), and if so, return d.
		final int tDivLimit = (int) (Math.cbrt(N));
		int i=0, p;
		while ((p = SMALL_PRIMES.getPrime(i++)) <= tDivLimit) {
			if (N%p==0) return p;
		}
		return N;
	}
}
