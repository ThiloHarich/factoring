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
 *
 * a^2 - 4m * k*N = b^2 mod 9
 * a^2 - b^2  = 4m * k*N  mod 9
 * a^2 = 0,1,4,7  mod 9
 * 4m = 0 mod 9 -> we have highest chances -> m=0 mod 9
 * it is sufficient to ensure
 * a^2 - 4m * k*N = 0,1,4,7 mod 9
 * a^2 = 4m * k*N + 0,1,4 mod 9
 * 0 -> 0 (0, 3, 6)
 * 1 -> 1 (1, 8)
 * 2 -> 4 (2, 5)
 * 3 -> 0
 * 4 -> 7 (4)
 * 5 -> 4
 * 6 -> 0
 * 7 -> 4
 * 8 -> 1
 *
 * mod 8:
 * 0 -> 0 (0,4)
 * 1 -> 1 (1,3,5,7) <- preferred solution
 * 2 -> 4 (2,6)
 * 3 -> 1
 * 4 -> 0
 * 5 -> 1
 * 6 -> 4
 * 7 -> 1
 *
 * -> we want a^2 - 4m * k*N = 1 mod 8
 * -> we want a^2  = 1 + 4m * k*N mod 8 , m,N
 * -> we want a^2  = 1 + 4k mod 8
 * 1) k even
 * a^2  = 1 mod 8  -> a odd, every second b fulfills the equation
 * 2) k odd has no solution
 * here b^2 must be 0,2 ->  2 out of 8 values for b will work

 * mod 16:
 * 0 ->  0 (0,4,8,12)
 * 1 ->  1 (1,15,7,9)
 * 2 ->  4 (2,6,10,14)
 * 3 ->  9 (3,5,11,13)
 * 4 ->  0
 * 5 ->  9
 * 6 ->  4
 * 7 ->  1
 * 8 ->  0
 * 9 ->  1
 *10 ->  4
 *11 ->  9
 *12 ->  0
 *13 ->  9
 *14 ->  4
 *15 ->  1
 * mod 32:
 * 0 ->  0 (0,8,16,24)
 * 1 ->  1 (1,15,17,31)
 * 2 ->  4 (2,6,10,14,18,22,26,30) <- preferred sol.
 * 3 ->  9 4*
 * 4 -> 16 (4,12,20,28)
 * 5 -> 25 4*
 * 6 ->  4
 * 7 -> 17 4*
 * 8 ->  0
 * 9 -> 17
 *10 ->  4
 *11 -> 25
 *12 -> 16
 *13 ->  9
 *14 ->  4
 *15 ->  1
 *16 ->  0
 * k = even
 * a^2 - 4m *k* N = 4, mod 32 , only possible for m*k*N = -1,0,3 mod 8
 * a^2  = 4 + 4m *k* N mod 32
 * a^2  = 4 * (1 + m *k* N) mod 32  ->  a^2 = 0,4,16 <-> a*a & 3 == 0
 * ->  1+ m*k * N = 0,1,4 mod 8
 * ->  m*k * N = 7,0,3 mod 8
 * m = 315 = 27 mod 32
 * -> a^2 - b^2 != 2,6
 * 4m * k*N != 2,6 mod 9
 *  k != (4*m*N) ^-1 * i mod 9 , i=2,6
 *  sei m= 5*7 = 35, N != 0 mod 3
 *  k != (4*35*N) ^-1 * i mod 9 , i=2,6
 *  k != 5^-1 * N ^-1 * i mod 9 , i=2,6
 *  k != 2 * N ^-1 * i mod 9 , i=2,6
 *  k != 4 * N ^-1  mod 9 and   <-> k*N != 4 mod 9
 *  k != 3 * N ^-1  mod 9 and   <-> k*N != 3 mod 9
 *  -> check it
 *  a^2 - 4* 35 * 4 = b^2 mod 9
 *  a^2 - 5 * 4 = b^2 mod 9
 *  a^2 - 5 * 3 = b^2 mod 9

 *  is this also working for other mods?
 */
public class HartSimple3 extends FactorAlgorithm {
	// Size of number is ~ 2^52
	private static final double DISCRIMINATOR = 1.0/(1<<10); // experimental result

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 * The multiplier ensures that the generated test values are always a square mod K_MULT
	 * Since the K_MULT consists out of 4 primes these numbers have a 2^4 = 16 times
	 * higher chance of being a square then random numbers. This is very helpful
	 */
	private static final int K_MULT = 3 * 3 * 5 * 7 /*11 * 13*/;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static double[] sqrt;

	// works for numbers up to 49 Bits
	private static int maxFactor = 1 << 21;

	//	squares mod 64 : 0,1,4,6,16,17,25,33,36,41,49,57
	//	a^2 - kN = squares
	//	a^2 = squares - kN
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
	public HartSimple3() {
		final Primes primesGen = Primes.initPrimesEratosthenes(maxFactor );
		primes = primesGen.primes;
		primesInv = primesGen.primesInv;
	}

	@Override
	public String getName() {
		return "HartSimple";
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
		final int multiplicationWorks =  Long.numberOfLeadingZeros(N) > 21 ? 0 : 1<<(21 - Long.numberOfLeadingZeros(N));

		int primeIndex = 0;
		//for the smallest primes we must do standard trial division
		for (; primes[primeIndex] < multiplicationWorks; primeIndex++) {
			if (N%primes[primeIndex]==0) {
				return primes[primeIndex];
			}
		}

		for (int sqrtIndex = 1, k = sqrtIndex * K_MULT; ;k += K_MULT, primeIndex++, sqrtIndex++) {
			// do trial division
			if ((long) (N * primesInv[primeIndex] + DISCRIMINATOR) * primes[primeIndex] == N)
				return primes[primeIndex];

			//			final long kN = k * N;
			//			if ((sqrtIndex & 1) ==  1 && (kN & 7) > 2 && (kN & 7) != 4)
			//				continue;

			final double aDouble = sqrt4N * sqrt[sqrtIndex];
			a = (long) (aDouble + ROUND_UP_DOUBLE);

			// a^2 - 4m * k*n = b^2   mod 9*64 ; 4m = 4*5*7 = 140
			// a^2 - 140 * k*n = b^2   mod 9*64 ; 4m = 4*5*7 = 140  , we can
			// adjust a
			if ((sqrtIndex & 1) == 0)
				a |= 1;
			else {
				final long kPlusN = k + N;
				if ((kPlusN & 3) == 0) {
					a += ((kPlusN - a) & 7);
				} else {
					final long adjust1 = (kPlusN - a) & 15;
					final long adjust2 = (-kPlusN - a) & 15;
					a += adjust1<adjust2 ? adjust1 : adjust2;
				}
			}

			test = a*a - k * fourN;
			b = (long) Math.sqrt(test);
			if (test == b*b && (gcd = gcdEngine.gcd(a+b, N))>1 && gcd < N) {
				//					System.out.print("," + test % 32);
				if ((b*b) % 32 == 4)
					System.out.print("," + a % 32);
				return gcd;
			}

		}
	}
}
