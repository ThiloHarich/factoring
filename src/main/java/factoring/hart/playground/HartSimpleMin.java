package factoring.hart.playground;

import java.math.BigInteger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.gcd.Gcd63;

public class HartSimpleMin extends FactorAlgorithm {

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
	// static double[] primes; with a table of primes the algorithm is much faster; here we just use odd numbers 3,5,7,9,..
	static double[] primesInv;

	// works for numbers up to 40 Bits
	private static int maxFactor = 1 << 19;

	static {
		// Precompute sqrts for all k required for N <= MAX_N and multiplier K_MULT
		sqrt = new double[maxFactor+1];
		primesInv = new double[maxFactor+1];
		for (int i = 1; i <= maxFactor; i++) {
			sqrt[i] = Math.sqrt(i*K_MULT);
			primesInv[i] = 1.0 / i;
		}
	}

	private final Gcd63 gcdEngine = new Gcd63();

	public HartSimpleMin() {
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

		for (int sqrtIndex = 1, k = sqrtIndex * K_MULT, prime = 3; ;k += K_MULT, prime += 2, sqrtIndex++) {
			// do trial division
			if ((long) (N * primesInv[prime] + DISCRIMINATOR) * prime == N)
				return prime;
			a = (long) (sqrt4N * sqrt[sqrtIndex] + ROUND_UP_DOUBLE);
			// adjust a mod 8
			if ((sqrtIndex & 1) == 0) {
				a |= 1;
			}
			else {
				final long kPlusN = k + N;
				a += (kPlusN & 3) == 0 ? ((kPlusN - a) & 7) : ((kPlusN - a) & 3);
			}
			test = a*a - k * fourN;
			b = (long) Math.sqrt(test);
			if (b*b == test && (gcd = gcdEngine.gcd(a+b, N))>1 && gcd < N) {
				return gcd;
			}
		}
	}
}
