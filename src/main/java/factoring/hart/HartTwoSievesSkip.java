package factoring.hart;

import static org.junit.Assert.assertTrue;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * Pretty simple yet fast variant of Hart's one line factorizer.
 * This implementations introduces some improvements that make it the fastest factoring algorithm
 * for numbers with more then ~20 and less then ~50 bit.
 * It avoids the complexity of calculating the square root when factoring multiple numbers,
 * by storing the square roots of k in the main loop. This has the highest performance impact.
 *
 * But there are some other improvements:
 *
 * It uses an optimized trial division algorithm to factorize small numbers.
 *
 * After calculating a number 'a' above sqrt(4*m*k) a will be adjusted to satisfy
 * some modulus a power of 2 argument.
 * It reuses the idea of rounding up by adding a well chosen constant (Warren D. Smith)
 *
 * We choose k to be a multiple of 315 = 3*3*5*7 and 45 = 3*3*5 this causes that
 * a^2 - 4kn = b^2 mod 3*3*5*7 or 3*3*5 which increases the chance to find a solution a^2 - 4kn = b^2 pretty much.
 * We iterate over k1 = 315 * i and k2 = 45 * i in parallel, but make sure that k2 != k1.
 *
 * General idea of this implementation:
 *
 * It tires to find solutions for a^2 - 4*m*i*n = b^2 from fermat we then know that
 * gcd(a+b, n) and gcd(a-b, n) are divisors of n.
 *
 * This is done by one simple loop over k were we generate numbers a = sqrt(4*m*k*n).
 * By storing sqrt(k) in an array this can be calculated fast.
 *
 * Compared with the regular Lehman algorithm, the Hart algorithm does not
 * need a second loop to iterate over the numbers 'a' for a given 'k' in the equation a^2 - 4kn.
 * This means that the upper bound for this loop - which would be a expensive sqrt call - does not has to be calculated.
 *
 * For each k the value sqrt(k) in order to determine a = ceil(sqrt(4kn))
 * the sqrt will be calculated only once and then stored in an array. This speeds up the sieving buy
 * a big factor since calculating the sqrt is expensive.
 *
 *
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class HartTwoSievesSkip extends FactorAlgorithm {

	private static final int SQUARES_MOD = 253; // 11*23
	//	private static final int SQUARES_MOD = 13; // 11*13*17*19
	public double hits = 0;
	public double total = 0;

	int mod = 11*13*17*19*23;

	private static final Logger LOG = Logger.getLogger(HartTwoSievesSkip.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT1 = 3*3*5*7;  // 315
	private static final int K_MULT2 = 3*3*5;    //  45

	/** Size of arrays this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	/**
	 * If the CPU supports https://en.wikipedia.org/wiki/Multiply–accumulate_operation#Fused_multiply–add set this to true
	 */
	private static final boolean USE_FUSED_MULTIPLY_ADD = true;

	private final boolean doTDivFirst;
	private final double[] sqrt1 = new double[I_MAX];
	private final double[] sqrt2 = new double[I_MAX];
	private final boolean [] sqrtMod = new boolean[mod];
	private final int [] k1s = new int[I_MAX];
	private final int [] k2s = new int[I_MAX];
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();


	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
	 * But the smaller possible factors get, it will become slower and slower.
	 */
	public HartTwoSievesSkip(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX

		for (int i=0, k =1, j =1; i<I_MAX-1; i++) {
			// for multiples of 3*3*5*7 we try to avoid multiples of 9
			//			if (i%9 != 0)
			{
				k1s[k] = i*K_MULT1;
				sqrt1[k] = Math.sqrt(i*K_MULT1);
				// a(k) = a(k-1) * (sqrt(1+1/(k-1)) - 1) + a(k-1)
				k++;
			}
			// for multiples of 3*3*5 we try to avoid number of the form 3*3*5*7 * i
			if (i%7 != 0) {
				k2s[j] = i * K_MULT2;
				sqrt2[j] = Math.sqrt(i * K_MULT2);
				// a(k) = a(k-1) * (sqrt(1+1/(k-1)) - 1) + a(k-1)
				//				sqrt1Add[j] = Math.sqrt(1.0 + 1.0/((i-1) * K_MULT2)) - 1;
				j++;
			}
		}
		for (int i = 0; i < sqrtMod.length; i++) {
			sqrtMod[(i*i) % mod] = true;
		}
	}

	@Override
	public String getName() {
		return "Hart_Fast2Mult(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger n) {
		return BigInteger.valueOf(findSingleFactor(n.longValue()));
	}

	/**
	 * Find a factor of long N.
	 * @param n
	 * @return factor of N
	 */
	public long findSingleFactor(long n) {
		if (doTDivFirst) {
			// do trial division before the Hart loop
			tdiv.setTestLimit((int) Math.cbrt(n));
			final long factor = tdiv.findSingleFactor(n);
			if (factor > 1) return factor;
		}

		final long fourN = n<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long gcd;
		long a = 0;
		long aOld = 0;
		int aMod = (int) a;
		// ensure k*n = 5 mod 8
		for (int i=1; ;i++) {
			// using the fusedMultiplyAdd operation defined in IEEE 754-2008 gives ~ 4-8 % speedup
			final int k1 = k1s[i];
			a = USE_FUSED_MULTIPLY_ADD ? (long)  Math.fma(sqrt4N, sqrt1[i], ROUND_UP_DOUBLE) : (long) (sqrt4N * sqrt1[i] + ROUND_UP_DOUBLE);
			a = adjustA(n, a, k1);
			final long aDiff = a - aOld;
			aMod += aDiff;
			aMod = aMod > mod ? aMod - mod : aMod;
			assertTrue(aMod < mod);
			final long test = a*a - k1 * fourN;
			final long b = (long) Math.sqrt(test);
			if (b*b == test && (gcd = gcdEngine.gcd(a+b, n))>1 && gcd<n) {
				assertTrue (sqrtMod[aMod]);
				return gcd;
			}
			aOld = a;
		}
	}

	/**
	 * Increases x to return the next possible solution for x for x^2 - 4kn = b^2.
	 * Due to performance reasons we give back solutions for this equations modulo a
	 * power of 2, since we can determine the solutions just by additions and binary
	 * operations.
	 * Based on mod arguments we can not prove that these are the only possible values
	 *
	 * if k is even x must be odd.
	 * if k*n == 3 mod 4 -> x = k*n+1 mod 8
	 * if k*n == 1 mod 8 -> x = k*n+1 mod 16 or -k*n+1 mod 16
	 * if k*n == 5 mod 8 -> x = k*n+1 mod 32 or -k*n+1 mod 32
	 *
	 * @param n
	 * @param x
	 * @param k
	 * @return
	 */
	private long adjustA(long n, long x, int k) {
		if ((k&1)==0)
			return x | 1;
		final long kNplus1 = k*n+1;
		if ((kNplus1 & 3) == 0)
		{
			return x + ((kNplus1 - x) & 7);
		}
		else if ((kNplus1 & 7) == 2) {
			final long adjust1 = ( kNplus1 - x) & 15;
			final long adjust2 = (-kNplus1 - x) & 15;
			return x + (adjust1 < adjust2 ? adjust1 : adjust2);
		}
		final long adjust1 = ( kNplus1 - x) & 31;
		final long adjust2 = (-kNplus1 - x) & 31;
		return x + (adjust1 < adjust2 ? adjust1 : adjust2);
	}
}
