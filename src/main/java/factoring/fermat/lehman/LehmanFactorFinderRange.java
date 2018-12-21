package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 * This is a version of the lehman factorization, which is a variant of the fermat
 * factorization.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster as the fastest (java) implementation I know, which is the yafu
 * lehman implementation found here https://github.com/DarkenCode/yafu/blob/master/factor/LehmanClean.c.
 * This project also contains the Java version of this here {@link LehmanYafuFact}.
 * This algorithm works for up to 41 bit. It is faster then other algorithms like SQUAFU for up to around 35 bits
 * for hard numbers with two or three (big) factors and random numbers. This is true for finding one factor and for
 * finding all factors.
 * It is well suited for numbers with small factors.
 *
 * The main features for the performance improvements over the yafu implementations is:
 *
 * It stores the square roots of the multiplier k. Thus the range for x (of x^2 - 4kn = y^2) can be calculated faster.
 * For k > n^1/3 / 16 the upper bound for x is less the
 * The lower bound plus 1. In this case at most one value of x has to be considered, but
 * the calculation of the lower and upper rage has to be done all the time.
 * Since calculating the ranges for x requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 * The calculation of a square root is much more time wasting then indexing the small array in a iterative way.
 *
 * It uses a version of trial division, where the multiple inverse of the primes are stored.
 * So instead of a division a multiplication is needed to find out if a number is dividable
 * by a prime. Since we already store the square roots for the lehman step, we are relaxed in using the space.
 *
 * It also uses an array of squares mod 1024 (a power of 2) and a product of small primes (different from 2)
 * to check if a number might be a square. So in most of the cases we get around calculating a squareroot of a number.
 *
 * This gives that the inner loop of the algorithm can be executed only by multiplications and index lookups in
 * small 2^14 = 16kb arrays.
 *
 *
 * The lehman equation x^2 - 4kn = y^2 is considered here for certain modulus.
 * Beside the arguments modulus 2 and 4 for calculating the step for x we added a new
 * When k*n % 4 == 3 we can increase x by 8. This comes out if you investigate
 * x^2 - 4kn = y^2 mod 32.
 *
 * If the switch from trial division to the lehman phase is done by n^1/3 ({@link #maxFactorMultiplier} == 1.
 * then the lehman phase (also) always finds prime factors.
 * So for maxFactorMultiplier == 1 this algorithm is well suited for finding all factors of numbers which are not
 * hard to factor, which is the case if the number is not having only prime factors above n^1/3.
 *
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * The main loop is build around the results mod 4.
 * It simply increases the value k, but knows for each kn mod 4 how the possible x look like.
 * It knows where the sequence x starts and how big the gap between two x are.
 *
 *
 * Open questions, possible improvements :
 * - can we get rid of storing the square roots? how can we calculate them efficiently?
 * - can we improve the range checks beside from the things we do in {@link LehmanFactorFinderMod12}?
 * - use a clever Pollard Rho variant for small factors
 *
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanFactorFinderRange implements FactorizationOfLongs, FactorFinder {

	// to be fast to decide if a number is a square we consider the
	// number modulo some values. First filer/modulus is a power of 2
	// since we can easily calculate the remainder by using bit and operation
	// here also low mod like 64 work fine as well
	//	static int mod_1024 = 1024;
	//
	//	private static boolean[] isSquareMod_1024;
	static int mod_1024 = 1024;

	private static boolean[] isSquareMod_1024;
	// for the other mod arguments we have to do a real expensive "%" operation
	// so we might have time to do a bit more expensive but space saving representation
	// of the squares

	float maxFactorMultiplier = 1.0f;
	float maxFactorMultiplierCube;
	boolean doTrialFirst;

	// This is a constant that is below 1 for rounding up double values to long
	protected static final double ROUND_UP_DOUBLE = 0.9999999665;

	double[] sqrt;
	double[] sqrtInv;
	// a fast way to do the trial division phase
	final FactorFinder smallFactoriser;
	int maxTrialFactor;

	static {
		isSquareMod_1024 = new boolean[mod_1024];
		for (int i = 0; i < mod_1024; i++) {
			isSquareMod_1024[(i * i) & 1023] = true;
		}
	}


	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 42 bits are not allowed since the precision of long does not allow more here
	 *
	 * @param bits                  the maximum number of bits the algorithm is going to work.
	 * @param maxFactorMultiplierIn Defines the size of the factors the algorithm is first looking for.
	 *                              The algorithm is working best for number where the lowest factor is
	 *                              maxFactorMultiplier * n^1/3 or higher. If maxFactorMultiplier is 1 or lower we first do trial division
	 *                              for numbers below n^1/3 then do the lehman factorizationByFactors above n^1/3. In this case the running time is
	 *                              bounded by max(maximal factor, c*n^1/3 / log(n)).
	 *                              For maxFactorMultiplier >= 1 we first inspect numbers above maxFactorMultiplier * n^1/3 with the lehman
	 *                              factorizationByFactors. After completing this phase we do
	 *                              trial division for numbers below maxFactorMultiplier * n^1/3<br>
	 *                              Use<br>
	 *                              maxFactorMultiplier <= 1<br>
	 *                              - if you do not know anything about the numbers to factorize<br>
	 *                              - if you know the maximal factors can not exceed maxFactorMultiplier * n^1/3 < n^1/3<br>
	 *                              here all found factors are prime factors and will be stored in the primeFactors Collection.
	 *                              maxFactorMultiplier = 3 (this is the optimal value) <br>
	 *                              - if you know for most of the numbers the maximal factors will exceed 3*n^1/3<br>
	 *                              In the last case {@link #findFactors(long, Collection)} might return a composite number
	 */
	public LehmanFactorFinderRange(int bits, float maxFactorMultiplierIn, boolean doTrialFirst) {
		if (bits > 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		maxFactorMultiplier = maxFactorMultiplierIn < 1 ? 1 : maxFactorMultiplierIn;
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.cbrt(1L << bits));
		maxFactorMultiplierCube = maxFactorMultiplier * maxFactorMultiplier * maxFactorMultiplier;
		smallFactoriser = new TrialInvFact(maxTrialFactor);
		this.doTrialFirst = doTrialFirst;
		initSquares();
	}

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
		final int kMax = (int) (maxTrialFactor / maxFactorMultiplierCube);

		sqrt = new double[kMax + 20];
		sqrtInv = new double[kMax + 20];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.pow(i, .5);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the square root and the inverse
			sqrtInv[i] = 1.0 / sqrtI;
		}
		System.out.println("sqrt tables for max trial factor " + maxFactorMultiplier + " built bytes used : " + sqrt.length);
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (n > 1l << 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		// with this implementation the lehman part is not slower then the trial division
		// we only apply the multiplier if we want to cut down the numbers to be tested
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.cbrt(n));
		//		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
		smallFactoriser.setMaxFactor(maxTrialFactor);
		// factor out all small factors if
		if (doTrialFirst) {
			final long nAfterTrial = smallFactoriser.findFactors(n, primeFactors);
			if (primeFactors == null && nAfterTrial != n)
				return nAfterTrial;
			n = nAfterTrial;
		}

		// if number is already factorized return immediately without calling lehman
		if (n == 1)
			return n;

		// re-adjust the maximal factor we have to search for. If small factors were found by trial division,
		// which is quite often the case for arbitrary numbers, this cuts down the runtime dramatically.
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.cbrt(n));
		// effectively kMax is reduced by maxFactorMultiplier^2 ->
		final int kMax = (int) (Math.ceil(maxTrialFactor / maxFactorMultiplierCube));
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		// we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
		final long nPow2Third = ((long) maxTrialFactor) * ((long) maxTrialFactor);
		final double nPow1Sixth = (nPow2Third / 4) / sqrtN;
		// surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
		// for k=2 we know that x has to be even and results in a factor more often
		final double sqrt4N = 2 * sqrtN;
		double xRange = Double.MAX_VALUE;
		for (int k = 1; k <= kMax ; k++) {
			final double sqrt4kn = sqrt[k] * sqrt4N;
			// adding a small constant to avoid rounding issues and rounding up is ~ 5% slower then
			// using the downcast and adding a constant close to 1. Took the constant from the yafu code
			//			final long xBegin = (long) Math.ceil(sqrt4kn);
			final long xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
			// use only multiplications instead of division here
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify.
			long xStep;
			long xFirst;
			final long kn = k * n;
			if (k % 2 == 0)
			{
				xStep = 2;
				xFirst = xBegin | 1l;
			}
			else
			{
				// this extra case gives ~ 10 %
				if (kn % 4 == 3) {
					xStep = 8;
					xFirst = (xBegin + ((7 - kn - xBegin) & 7));
				} else
				{
					xStep = 4;
					xFirst = xBegin + ((k + nMod4 - xBegin) & 3);
				}
			}
			xRange = nPow1Sixth * sqrtInv[k];
			// if we do not cast to long we have to compare long with double ->  50 % slower
			final long xEnd =  (long) (sqrt4kn + xRange);
			// a maxFactorMultiplier > 2 ensures that xEnd >= xBegin + 2
			for (long x = xFirst; x <= xEnd; x += xStep)
			{
				final long xSquare = x * x;
				final long right = xSquare - kn * 4;
				// instead of taking the square root (which is a very expensive operation)
				// and squaring it, we do some mod arguments to filter out non squares
				if (isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					if (y * y == right) {
						final long factor = PrimeMath.gcd(n, x - y);
						if (factor != 1) {
							if (findsPrimesOnly())
							{
								primeFactors.add(factor);
								if (n != factor)
									primeFactors.add(n / factor);
								return 1;
							}
							else
								// if we have not done trial division first, the factor might be composite -> we just return the factor
								return factor;
						}
					}
				}
			}
		}
		// if we have not found a factor we still have to do the trial division phase
		if (!doTrialFirst)
			n = smallFactoriser.findFactors(n, primeFactors);

		return n;
	}

	protected static boolean isSquare(final long number) {
		if (isSquareMod_1024[(int) (number & 1023)]) {
			final long y = (long) Math.sqrt(number);
			return y * y == number;
		}
		return false;
	}

	@Override
	public boolean findsPrimesOnly(){
		return doTrialFirst;
	}

	@Override
	public boolean returnsCompositeOnly(){
		return false;
	}

	@Override
	public String toString() {
		return "LehmanFactorFinderRange{" +
				"maxFactorMultiplier=" + maxFactorMultiplier +
				'}';
	}
}
