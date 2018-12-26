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
public class Lehman_TillSuperSimple2 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;


	private final Gcd63 gcdEngine = new Gcd63();

	private double[] sqrt;

	private double sqrt4N;

	private long fourN;

	long N;

	public Lehman_TillSuperSimple2() {
		initSqrts();
	}

	private void initSqrts() {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = (int) (6*Math.cbrt(1L<<50) + 1);
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
		this.N = N;
		final int kLimit = (int)  Math.ceil(Math.cbrt(N));;

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);

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
		for (int k=1; k < 8*kLimit; k++) {
			long factor;

			//			 0 mod 6 will be check double as often as the rest
			if ((factor = lehmanEven(12 * k)) > 1 && factor < N)
				return factor;

			if ((factor = lehmanEven(12 * k + 6)) > 1 && factor < N)
				return factor;

			if ((factor = lehmanOdd(6 * k + 3)) > 1 && factor < N)
				return factor;

			//			if ((factor = lehmanEven(6 * k + 4)) > 1 && factor < N)
			//				return factor;
			//
			//			if ((factor = lehmanEven(6 * k + 2)) > 1 && factor < N)
			//				return factor;


			// Here k is always 1 and we increase 'a' by aStep.
			final long test = aForK1*aForK1 - fourN;
			final long b = (long) Math.sqrt(test);
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

	/**
	 * This method tries to choose an a such that n + k = a mod 4.
	 * This can be applied for odd kBegin.
	 * If a*a - 4 * k * N is a square it returns true.
	 * @param kBegin
	 * @param kLimit
	 * @return
	 */
	long lehmanOdd(int k) {
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
		return -1;
	}

	/**
	 * makes a odd can be applied to even kBegin.
	 * @param kBegin
	 * @param kEnd
	 * @return
	 */
	long lehmanEven(int k) {
		// for k = 0 mod 6 a must be odd
		final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
		final long test = a*a - k * fourN;
		final long b = (long) Math.sqrt(test);
		if (b*b == test) {
			return gcdEngine.gcd(a+b, N);
		}
		return -1;
	}

}
