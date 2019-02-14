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
import factoring.trial.TrialMultiplyUnrol;

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
public class HartFlexible extends FactorAlgorithm {

	// Size of number is ~ 2^52
	private static final int DISCRIMINATOR_BITS = 10; // experimental result
	private static final double DISCRIMINATOR = 1.0/(1<<DISCRIMINATOR_BITS);

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

	private final int trialDivisonWhen;

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
	private final TrialMultiplyUnrol trialDivision = new TrialMultiplyUnrol(maxFactor);
	double[] primesInv;
	int[] primes;

	/**
	 *
	 * @param isHardSemiprime Set this to false for most of the numbers.
	 * You should set this to true only if it is a semiprime with both factors near n^1/2.
	 */
	public HartFlexible(int trialDivisonWhen) {
		final Primes primesGen = Primes.initPrimesEratosthenes(3 * maxFactor );
		primes = primesGen.primes;
		primesInv = primesGen.primesInv;
		this.trialDivisonWhen = trialDivisonWhen;
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
		long nDivPrime;
		long multiplicationWorks = 0;
		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);

		if (trialDivisonWhen <= 0) {
			final int testLimit = (int) Math.cbrt(N);
			trialDivision.setTestLimit(testLimit);
			final long factor = trialDivision.findSingleFactor(N);
			if (factor < N && factor > 1)
				return factor;
			else {
				final int NBits = 64-Long.numberOfLeadingZeros(N);
				final int multiplicationWorksBits = NBits - 53 + DISCRIMINATOR_BITS;
				multiplicationWorks = Math.min(1<<multiplicationWorksBits, 0);
			}
		}
		int primeIndex = -1;
		//for the smallest primes we must do standard trial division
		for (; primes[++primeIndex] < multiplicationWorks;) {
			if (N%primes[primeIndex]==0) {
				return primes[primeIndex];
			}
		}

		int i = (((K_MULT + N) & 3) != 0) ? 1 : 3;
		try {
			for (int k = i * K_MULT; ;) {
				// do trial division
				for (int adjustAMod = 3; adjustAMod <= 7; adjustAMod += 4) {
					// for most of the numbers it is beneficial to check for small primes multiple times
					// a loop is much slower then unrolling it for unkown reasons.
					nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
					if (nDivPrime * primes[primeIndex] == N) {
						return primes[primeIndex];
					}
					//						nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
					//						if (nDivPrime * primes[primeIndex] == N) {
					//							return primes[primeIndex];
					//						}
					if (trialDivisonWhen != 2) {
						nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
						if (nDivPrime * primes[primeIndex] == N) {
							return primes[primeIndex];
						}
						nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
						if (nDivPrime * primes[primeIndex] == N) {
							return primes[primeIndex];
						}
						nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
						if (nDivPrime * primes[primeIndex] == N) {
							return primes[primeIndex];
						}
						//						nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
						//						if (nDivPrime * primes[primeIndex] == N) {
						//							return primes[primeIndex];
						//						}
					}


					a = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE);
					// k odd -> a + k + N = 0 mod 4, k + N = 0 mod 8 if k+N = 0 mod 4, we adjust a accordingly
					a += (((k + N) - a) & adjustAMod);
					test = a*a - k * fourN;
					b = (long) Math.sqrt(test);
					if (b*b == test && (gcd = gcdEngine.gcd(a+b, N)) > 1 && gcd < N) {
						return gcd;
					}
					k += K_MULT;

					// even k -> a must be odd
					a = (long) (sqrt4N * sqrt[i++] + ROUND_UP_DOUBLE) | 1L;
					test = a*a - k * fourN;
					b = (long) Math.sqrt(test);
					if (b*b == test && (gcd = gcdEngine.gcd(a+b, N)) > 1 && gcd < N) {
						return gcd;
					}
					k += K_MULT;

					nDivPrime = (long) (N * primesInv[++primeIndex] + DISCRIMINATOR);
					if (nDivPrime * primes[primeIndex] == N) {
						return primes[primeIndex];
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// should never happen in this implementation; if it does then N > MAX_N
			System.out.println(this.getClass().getSimpleName() + " failed to factor N=" + N + ". Cause: " + e);
			return 0;
		}
	}
}
