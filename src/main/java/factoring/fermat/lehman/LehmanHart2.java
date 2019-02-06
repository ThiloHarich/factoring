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

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.trial.TrialMultiplyCorrection;

public class LehmanHart2 extends FactorAlgorithmBase {

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static double[] sqrt, sqrtInv;


	static {
		// Precompute sqrts for all possible k. 2^21 entries are enough for N~2^63.
		final int kMax = 1<<21;
		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0/sqrtI;
		}
	}

	private static final TrialMultiplyCorrection trialDivision = new TrialMultiplyCorrection(1<<21);

	private long N;
	private long fourN;
	private double sqrt4N;
	private final Gcd63 gcdEngine = new Gcd63();

	private final double sqrt35 = Math.sqrt(35);


	public LehmanHart2() {
	}

	@Override
	public String getName() {
		return "LehmanHart()";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	/**
	 *
	 * @param N
	 * @return
	 */
	public long findSingleFactor(long N) {
		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		// do trial division before Lehman loop ?
		long factor;

		fourN = N<<2;
		sqrt4N = Math.sqrt(fourN);

		trialDivision.setMaxFactor(cbrt);
		if ((factor=trialDivision.findFactor(N))>1)
			return factor;

		long a,b,test;
		int k = 6;
		for (; ;) {
			a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) | 1L;
			test = a*a - k * fourN;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				factor = gcdEngine.gcd(a+b, N);
				if (factor < N)
					return factor;
			}
			//			a = (long) (sqrt4N * sqrt35  * sqrt[k] + ROUND_UP_DOUBLE) | 1L;
			//			test = a*a - k * 35 * fourN;
			//			b = (long) Math.sqrt(test);
			//			if (b*b == test) {
			//				factor = gcdEngine.gcd(a+b, N);
			//				if (factor < N)
			//					return factor;
			//			}
			k += 3;
			a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			a = adjustMod8(a, k + N);
			test = a*a - (k) * fourN;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				factor = gcdEngine.gcd(a+b, N);
				if (factor < N)
					return factor;
			}
			k += 3;
		}
	}


	/**
	 * Here we make sure that for k even a fulfills
	 * a^2 - k*N = y^2 mod 8.
	 * k + n == a mod 8 if k+n == 0 mod 4 or else
	 * k + n == a mod 4
	 * @param a
	 * @param kPlusN k + N
	 * @return
	 */
	private static long adjustMod8(long a, long kPlusN) {
		if ((kPlusN & 3) == 0)
			return a + ((kPlusN - a) & 7);
		return a + ((kPlusN - a) & 3);
	}

}
