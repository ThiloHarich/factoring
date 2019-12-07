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
package factoring;

import static de.tilman_neumann.jml.base.BigIntConstants.I_1;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.TestsetGenerator;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.ConfigUtil;

/**
 * Simple implementation of Lehmans factor algorithm following https://programmingpraxis.com/2017/08/22/lehmans-factoring-algorithm/,
 * useful for illustrating the basic algorithm.
 * Works for N <= 45 bit.
 *
 * @author Tilman Neumann
 */
public class Lehman_Analyzer3 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(Lehman_Analyzer3.class);

	// algorithm options
	/** number of test numbers */
	private static final int N_COUNT = 100000;
	/** the bit size of N to start with */
	private static final int START_BITS = 40;
	/** the increment in bit size from test set to test set */
	private static final int INCR_BITS = 1;
	/** maximum number of bits to test (no maximum if null) */
	private static final Integer MAX_BITS = 63;

	private final Gcd63 gcdEngine = new Gcd63();

	private final Set<Integer>[][] aValues;

	private final int[] numSol;
	private final int[][] numSol2;

	// 0, 6, 12, 18, 24
	private static final int MOD = 30;

	public Lehman_Analyzer3() {
		aValues = new SortedSet[MOD][MOD];
		for (int i=0; i<MOD; i++) {
			aValues[i] = new SortedSet[MOD];
			for (int j=0; j<MOD; j++) {
				aValues[i][j] = new TreeSet<>();
			}
		}
		numSol = new int[MOD];
		numSol2 = new int[MOD][MOD];
	}

	@Override
	public String getName() {
		return "Lehman_Analyzer3";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		final int cbrt = (int) Math.ceil(Math.cbrt(N));
		final double sixthRoot = Math.pow(N, 1/6.0); // double precision is required for stability

		final int scm = MOD;
		//		final int kLimit = (cbrt + scm) / scm * scm;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoA = (cbrt >> 6);
		// twoA = 0 mod 6
		kTwoA = ((kTwoA + scm)/ scm) * scm;

		for (int k=1; k < kTwoA; k+= 1) {
			final long fourKN = k*N<<2;
			final double fourSqrtK = Math.sqrt(k<<4);
			final double sqrt4kN = Math.sqrt(fourKN);
			final int sqrt4kNCeil = (int) Math.ceil(sqrt4kN); // ceil() is required for stability
			final int limit = (int) (sqrt4kN + sixthRoot / fourSqrtK);
			for (int a = sqrt4kNCeil; a <= limit; a++) {
				final long test = a*(long)a - fourKN;
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		//		for (int k = kTwoA; k <= 2*cbrt; k += 6) {
		//			final long fourKN = k*N<<2;
		//			final long a = ((long) Math.ceil(Math.sqrt(fourKN))) | 1; // ceil() is required for stability
		//			// for k = 0 mod 6 a must be odd
		//			final long test = a*a - fourKN;
		//			final long b = (long) Math.sqrt(test);
		//			if (b*b == test) {
		//				//				System.out.print(k/6 & 3);
		//				return gcdEngine.gcd(a+b, N);
		//			}
		//		}

		//		for (int k = kTwoA; k <= cbrt << 2; k += 12) {
		//			final long fourKN = k*N<<2;
		//			final long a = ((long) Math.ceil(Math.sqrt(fourKN))) | 1; // ceil() is required for stability
		//			// for k = 0 mod 6 a must be odd
		//			final long test = a*a - fourKN;
		//			final long b = (long) Math.sqrt(test);
		//			if (b*b == test) {
		//				//				System.out.print(k/6 & 3);
		//				return gcdEngine.gcd(a+b, N);
		//			}
		//		}

		//		for (int k = kTwoA + 3; k <= cbrt; k += 6) {
		//			final long fourKN = k*N<<2;
		//			long a = ((long) Math.ceil(Math.sqrt(fourKN))); // ceil() is required for stability
		//			final long kPlusN = k + N;
		//			if ((kPlusN & 3) == 0) {
		//				a += ((kPlusN - a) & 7);
		//			} else
		//			{
		//				a += ((kPlusN - a) & 3);
		//			}
		//			final long test = a*a - fourKN;
		//			final long b = (long) Math.sqrt(test);
		//			if (b*b == test) {
		//				//				System.out.print(k/6 & 3);
		//				return gcdEngine.gcd(a+b, N);
		//			}
		//		}

		for (int k=kTwoA; k <= cbrt; k+= 1) {
			final long fourKN = k*N<<2;
			final double fourSqrtK = Math.sqrt(k<<4);
			final double sqrt4kN = Math.sqrt(fourKN);
			final int sqrt4kNCeil = (int) Math.ceil(sqrt4kN); // ceil() is required for stability
			final double range = sixthRoot / fourSqrtK;
			final int limit = (int) (sqrt4kN + range);
			for (int a = sqrt4kNCeil; a <= limit; a++) {
				final long test = a*(long)a - fourKN;
				final long b = (long) Math.sqrt(test);
				if (b*b == test) {
					aValues[k%MOD][(int) ((N+ k)%MOD)].add(a%MOD);
					//					aValues[(k)%MOD][(int) ((k*N)%MOD)].add(a%MOD);
					//					aValues[(k)%MOD][0].add(a%MOD);
					numSol[(k)%MOD]++;
					//					numSol2[k%MOD][(int) ((4*N*k)%MOD)]++;
					numSol2[k%MOD][0]++;
					//					numSol2[k%MOD][(int) ((k*N)%MOD)]++;
					return gcdEngine.gcd(a+b, N);
				}
			}
		}

		// Nothing found. Either N is prime or the implementation is buggy. For N > 45 bit it won't work.
		return 0;
	}

	private void testRange(int bits) {
		final BigInteger N_min = I_1.shiftLeft(bits-1);
		// find N-set for square tests
		//ArrayList NSet = TestsetGenerator.generate(bits, N_COUNT);
		final ArrayList<BigInteger> NSet = TestsetGenerator.generate(bits, N_COUNT);
		LOG.info("Test N with " + bits + " bits, i.e. N >= " + N_min);

		for (final BigInteger N : NSet) {
			final BigInteger factor = this.findSingleFactor(N);
		}

		for (int i=0; i<MOD; i++) {
			boolean logged = false;
			for (int j=0; j<MOD; j++) {
				if (aValues[i][j].size() > 0) {
					//					LOG.info("Success a-values %" + MOD + " for k %" + MOD + "==" + i + ", (4nk)%" + MOD + "==" + j + " num all : " + numSol[i] + " num : " + numSol2[i][j] +  ", aval : "+ aValues[i][j]);
					LOG.info("Success a-values % " + MOD + " for k %" + MOD + " =="  + i + ",  k*n %" + MOD + " == " + j + " num all : " + numSol[i] + " num : " + numSol2[i][j] +  ", aval : "+ aValues[i][j]);
					logged = true;
				}
			}
			//			if (logged) LOG.info("");
		}

	}

	public static void main(String[] args) {
		ConfigUtil.initProject();
		int bits = START_BITS;
		while (true) {
			// test N with the given number of bits, i.e. 2^(bits-1) <= N <= (2^bits)-1
			final Lehman_Analyzer3 testEngine = new Lehman_Analyzer3();
			testEngine.testRange(bits);
			bits += INCR_BITS;
			if (MAX_BITS!=null && bits > MAX_BITS) break;
		}
	}
}
