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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.ConfigUtil;

/**
 * Faster implementation of Lehman's factor algorithm.
 * Works flawlessly for N <= 56 bit.
 *
 * This version will do trial division before or after the main loop depending on factorSemiprimes.
 *
 * @authors Tilman Neumann + Thilo Harich
 */
public class Lehman_Fast33 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_Fast33.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private static final int DISCRIMINATOR_BITS = 10; // experimental result
	private static final double DISCRIMINATOR = 1.0/(1<<DISCRIMINATOR_BITS);

	private static double[] sqrt, sqrtInv;

	private static double[] primesInv;
	private static int[] primes;


	static {
		// Precompute sqrts for all possible k. 2^22 entries are enough for N~2^66.
		final long start = System.currentTimeMillis();
		final int kMax = 1<<25;
		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0/sqrtI;
		}
		final long end = System.currentTimeMillis();
		System.out.println("Time Init sqrt   : + " + (end - start));
		// finds the prime factors up to kMax by the sieve of eratosthenes.
		// Not optimized, since this is only called once when initializing.
		// needs n * log(log(n)) time
		//		kMax = 1 << 19;
		final double logMaxFactor = Math.log(kMax);
		final int maxPrimeIndex = (int) ((kMax) / (logMaxFactor - 1.1)) + 4;
		primesInv = new double [maxPrimeIndex];
		primes = new int [maxPrimeIndex];
		int primeIndex = 0;
		final boolean [] noPrimes = new boolean [kMax];
		for (int i = 2; i <= Math.sqrt(kMax); i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
			for (int j = i * i; j < kMax; j += i) {
				noPrimes[j] = true;
			}
		}
		for (int i = (int) (Math.sqrt(kMax)+1); i < kMax; i++) {
			if (!noPrimes[i]) {
				primes[primeIndex] = i;
				primesInv[primeIndex++] = 1.0 / i;
			}
		}
		for (int i=primeIndex; i < primes.length; i++) {
			primes[i] = Integer.MAX_VALUE;
		}
		final long end2 = System.currentTimeMillis();
		System.out.println("Time Init primes : + " + (end2 - end));
	}

	long N;
	long fourN;
	double sqrt4N;
	double sqrt24N;
	double sqrt120N;
	boolean factorSemiprimes;
	private final Gcd63 gcdEngine = new Gcd63();

	private long n24;
	private long N120;

	public Multimap<Integer, Integer> remainders = ArrayListMultimap.create();

	/**
	 * Only constructor.
	 * @param factorSemiprimes if the number to be factored might be a semiprimes were each factor is greater then n^1/3
	 * then set this to true. The algorithm will always work. But this might have a very positive effect on the performance.
	 */
	public Lehman_Fast33(boolean factorSemiprimes) {
		this.factorSemiprimes = factorSemiprimes;
	}

	@Override
	public String getName() {
		return "Lehman_Fast";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		this.N = N;
		final int cbrt = (int) Math.cbrt(N);

		long factor;
		if (!factorSemiprimes && (factor = findSmallFactors(N, cbrt)) != N)
			return factor;

		// limit for must be 0 mod 6, since we also want to search above of it
		final int bigStep = 30;
		final int kLimitDiv30 = (cbrt + bigStep) / bigStep;
		final int kLimitDiv6 = kLimitDiv30 * 5;
		// For kLimit / 64 the range for a is at most 2, this is what we can ensure
		int kTwoADiv6 = (cbrt >> 6);
		// twoA = 0 mod 5 -> twoA*6 = 0 mod 30
		final int kTwoADiv30 = (kTwoADiv6 + bigStep)/ bigStep;
		kTwoADiv6 = kTwoADiv30 * 5;
		fourN = N<<2;
		n24 = N* 24;
		N120 = N * 120;
		sqrt4N = Math.sqrt(fourN);
		sqrt24N = Math.sqrt(n24);
		sqrt120N = Math.sqrt(N120);
		final double sixthRootTerm = 0.25 * Math.pow(N, 1/6.0); // double precision is required for stability

		// first investigate in solutions a^2 - sqrt(k*n) = y^2 were we only have two possible 'a' values per k
		// but do not go to far, since there we have a lower chance to find a factor
		// here we only inspect k = 30*i since they much more likely have a solution a^2 - sqrt(k*n) = y^2
		// this is we use a multiplier of 30 for odd 'a' we have much higher chances to find solutions

		factor = lehman30(kTwoADiv30, kLimitDiv30);
		if (factor>1 && factor != N)
			return factor;

		// inspect k = 6*i != 30*j
		factor = lehman6(kTwoADiv6+1, kLimitDiv6);
		if (factor>1 && factor != N)
			return factor;

		final int kTwoA = kTwoADiv6*6;
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
				} else {
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

		// So far we have only odd solutions for 'a' now we try to get all the even solutions
		factor = lehmanOdd(kTwoA + 3, 6 * kLimitDiv6);
		if (factor > 1)
			return factor;

		// continue even k loop. Theoretically we might have missed solutions for k = 1,2,4,5 mod 6.
		// but when looking up to 6 * kLimit it seems to be we are not missing one of them.
		// to be sure we might just increase the limit even further

		factor = lehman6(kLimitDiv6 + 1, 6 * kLimitDiv6);
		if (factor > 1)
			return factor;

		// seems like we do not need this
		factor = lehman30(kLimitDiv30, kLimitDiv30 * 6);
		if (factor>1 && factor != N)
			return factor;


		// Check via trial division whether N has a nontrivial divisor d <= cbrt(N).
		//LOG.debug("entering tdiv...");
		factor = findSmallFactors(N, cbrt);
		if (factor>1 && factor<N) return factor;

		for (int k=kTwoA + 1; k <= 6 *kLimitDiv6; k++) {
			final long fourKN = k*N<<2;
			final long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE) - 1;
			final long test = a*a - fourKN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				final long gcd = gcdEngine.gcd(a+b, N);
				if (gcd>1 && gcd<N) {
					return gcd;
				}
			}
		}

		return 0; // fail
	}

	private long lehmanOdd(int kBegin, final int kLimit) {
		for (int k = kBegin; k <= kLimit; k += 6) {
			long a = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			// for k odd a must be even and k + n + a = 0 mod 4
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else
			{
				a += ((kPlusN - a) & 3);
			}
			final long test = a*a - k * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}

		return -1;
	}

	private long lehman6(int kBeginIndex, final int kEnd) {
		for (int k = kBeginIndex ; k <= kEnd;) {
			long a  = (long) (sqrt24N * sqrt[k]   + ROUND_UP_DOUBLE) | 1;
			long k24N = k++ * n24;
			long test = a*a - k24N;
			long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
			k24N += n24;
			a = (long) (sqrt24N * sqrt[k++] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
			k24N += n24;
			a = (long) (sqrt24N * sqrt[k++] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
			k24N += n24;
			a = (long) (sqrt24N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			test = a*a - k24N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
			// jump over k = 0 mod 5
			k += 2;
		}
		return -1;
	}

	private long lehman30(int kBeginIndex, final int kEnd) {
		for (int k = kBeginIndex ; k <= kEnd; k++) {
			//			System.out.print("," + k);
			final long a  = (long) (sqrt120N * sqrt[k] + ROUND_UP_DOUBLE) | 1;
			final long test = a*a - k * N120;
			final long b = (long) Math.sqrt(test);
			if (b*b == test) {
				return gcdEngine.gcd(a+b, N);
			}
		}
		//		System.out.println();
		return -1;
	}


	/**
	 * Finds factors up to maxFactor by a modified trial division.
	 * Checks if n is divisible by p by multiplying n by 1/p.
	 * The double values of 1/p are pre-computed.
	 * @param n
	 * @param maxFactor
	 * @return
	 */
	public long findSmallFactors(long n, int maxFactor) {
		for (int primeIndex = 1; primes[primeIndex] <= maxFactor; primeIndex++) {
			final double nDivPrime = n * primesInv[primeIndex];
			if (((long)(nDivPrime+DISCRIMINATOR)) - ((long)(nDivPrime-DISCRIMINATOR)) == 1 &&
					n > 1 && n % primes[primeIndex] == 0) {
				return primes[primeIndex];
			}
		}
		return n;
	}


	/**
	 * Test.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		ConfigUtil.initProject();

		// These test number were too hard for previous versions:
		final long[] testNumbers = new long[] {
				646131439L,
				18019629551563L,
				// stuff Lehman_TillSimple3 did not like
				19699548984827L,
				5640012124823L,
				7336014366011L,
				52199161732031L,
				73891306919159L,
				112454098638991L,

				// stuff Lehman_TillSimple3_2 did not like
				32427229648727L,
				87008511088033L,
				92295512906873L,
				338719143795073L,
				346425669865991L,
				1058244082458461L,
				1773019201473077L,
				6150742154616377L,

				// stuff Lehman_Fast2, Lehman_Fast3 did not like
				44843649362329L,
				67954151927287L,
				134170056884573L,
				198589283218993L,
				737091621253457L,
				1112268234497993L,
				2986396307326613L,

				// stuff Lehman_TillSimple3_3 did not like
				26275638086419L,
				62246008190941L,
				209195243701823L,
				290236682491211L,
				485069046631849L,
				1239671094365611L,
				2815471543494793L,
				5682546780292609L,
		};

		final Lehman_Fast33 lehman = new Lehman_Fast33(true);
		for (final long N : testNumbers) {
			final long factor = lehman.findSingleFactor(N);
			LOG.info("N=" + N + " has factor " + factor);
		}
	}
}
