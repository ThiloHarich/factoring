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
package factoring.hart;

import java.math.BigInteger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.primes.Primes;

/**
 * Pretty simple yet fast variant of Hart's one line factorizer.
 *
 * When called with doTDivFirst=false, this variant is marginally slower than Hart_Fast_HardSemiprimes
 * for hard semiprimes, but much better on random composites.
 *
 * If test numbers are known to be random composites, then doTDivFirst=true will improve performance significantly.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class HartSimple extends FactorAlgorithm {
	// Size of number is ~ 2^52
	private static final double DISCRIMINATOR = 1.0/(1<<10); // experimental result

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 * The multiplier ensures that the generated test values are always a square mod K_MULT
	 * Since the K_MULT consists out of 4 primes these numbers have a 2^4 = 16 times
	 * higher chance of being a square then random numbers. This is very helpful
	 */
	private static final int K_MULT = 3 * 3 * 5* 7;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static double[] sqrt;

	private static int maxFactor = 1 << 19;

	static {
		// Precompute sqrts for all k required for N <= MAX_N and multiplier K_MULT
		sqrt = new double[maxFactor+1];
		for (int i = 1; i <= maxFactor; i++) {
			sqrt[i] = Math.sqrt(i*K_MULT);
		}
		System.out.println("Hart_Fast: Initialized sqrt array with " + maxFactor + " entries");
	}

	private final Gcd63 gcdEngine = new Gcd63();
	double[] primesInv;
	int[] primes;

	/**
	 *
	 * @param isHardSemiprime Set this to false for most of the numbers.
	 * You should set this to true only if it is a semiprime with both factors near n^1/2.
	 */
	public HartSimple() {
		final Primes primesGen = Primes.initPrimesEratosthenes(maxFactor );
		primes = primesGen.primes;
		primesInv = primesGen.primesInv;
	}

	@Override
	public String getName() {
		return "HartMod9";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	/**
	 * Find a factor of long N.
	 * @param N
	 * @return factor of N
	 */
	public long findSingleFactor(long N) {
		long a,b,test, gcd;
		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		final int multiplicationWorks = Math.min (1<<(21 - Long.numberOfLeadingZeros(N)), 0);

		int primeIndex = 0;
		//for the smallest primes we must do standard trial division
		for (; primes[primeIndex] < multiplicationWorks; primeIndex++) {
			if (N%primes[primeIndex]==0) {
				return primes[primeIndex];
			}
		}

		for (int i = 1, k = i * K_MULT; ;k += K_MULT, primeIndex++, i++) {
			// do trial division
			if ((long) (N * primesInv[primeIndex] + DISCRIMINATOR) * primes[primeIndex] == N) {
				return primes[primeIndex];
			}
			a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
			// adjust a
			if ((i & 1) == 0)
				a |= 1;
			else {
				final long kPlusN = k + N;
				a += (kPlusN & 3) == 0 ? ((kPlusN - a) & 7) : ((kPlusN - a) & 3);
			}
			test = a*a - k * fourN;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd < N)
					return gcd;
			}
		}
	}
}
