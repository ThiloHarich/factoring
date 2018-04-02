package factoring.fermat.lehman;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

import java.math.BigInteger;
import java.util.Collection;

/**
 * This is a version of the lehman factorizationByFactors, which is a variant of the fermat
 * factorizationByFactors.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster then the java version of the yafu lehman factorizationByFactors.
 * By storing the square roots of the multiplier k the range can be done faster.
 * It also uses a version of trial division, where the multiple inverse of the primes are stored.
 * So instead of a division a multiplication is needed to find out if a number is dividable
 * by a prime.
 * In the lehman algorithm in 15 out of 16 cases (for k) only one value of x has to be considered, but
 * the calculation of the lower and upper rage has to be done all the time.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 *
 * Like in the YAFU implementation we get no speed when using smooth multipliers first (like lehman has suggested it).
 * The Hart variant always just one x per multiplier k, this eliminates the determination of the
 * upper bound, but using it gives no extra speed.
 *
 * Open questions, possible improvements :
 * - can we get rid of storing the square roots? how can we calculate them efficiently?
 * - When we know sqrt(x) we can conclude sqrt(i*x) but seems to give no speed
 * - In most of the cases (more then 93%) for k only one or none the value x^2 -n has to be calculated.
 * Can we find them without trying if the value sqrt(k*n) fits in the range?
 *
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanBigFact extends FindPrimeFact {

	static double ONE_THIRD = 1.0/3;

//	float balanceTrial = 1.0f;
//	float balanceTrialCube;

	// This is a constant that is below 1 for rounding up double values to long
	protected static final double ROUND_UP_DOUBLE = 0.9999999665;
	// TODO see if it works for values above the precision/scale of the algorithm as well
	// If not check BigInteger
	// to initialize the square roots we have to determine the maximal k

	double [] sqrt;
	double [] sqrtInv;

	final TrialInvFact smallFactoriser;
	int maxTrialFactor;

	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 48 bits are not allowed since the precision of long does not allow more here
	 * @param bits
	 */
	public LehmanBigFact(int bits, float balanceTrial) {
//		maxTrialFactor = (int) Math.ceil(balanceTrial * Math.pow(1L << bits, ONE_THIRD));
		maxTrialFactor = (int) Math.ceil(Math.pow(1L << bits, ONE_THIRD));
		// using the trial division algorithm more doe not help
//		this.balanceTrial = balanceTrial;
//		balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
//        smallFactoriser = new TrialPrimesDynamicFact(maxTrialFactor);
        smallFactoriser = new TrialInvFact(maxTrialFactor);

		initSquares();
	}

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
//		float balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
//		int kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
		int kMax = (int) (Math.ceil(maxTrialFactor));

		sqrt = new double[kMax + 1];
		sqrtInv = new double[kMax + 1];
		for(int i = 1; i< sqrt.length; i++)
		{
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the  root and the inverse
			sqrtInv[i] = 1.0d / sqrtI;
		}
		System.out.printf(" sqrt table[0..%d] built: ", sqrt.length);
		for(int i=0; i<Math.min(5, maxTrialFactor); i++){
			System.out.printf("%f ", sqrt[i]);
		}
		System.out.printf(" ... %f %f\n", sqrt[sqrt.length-2],sqrt[sqrt.length-1]);
	}

	@Override
	public long findPrimeFactors(long nOrig, Collection<Long> primeFactors) {
		// with this implementation the lehman part is not slower then the trial division
		// we do not have to use a multiplier for the maximal factor were we apply the
		// trial division phase
//		double maxTrialFactor =  Math.ceil(balanceTrial * Math.pow(nOrig, ONE_THIRD));
		double maxTrialFactor =  Math.ceil(Math.pow(nOrig, ONE_THIRD));
		smallFactoriser.setMaxFactor((int) maxTrialFactor);
		// factor out all small factors
		final long n = smallFactoriser.findPrimeFactors(nOrig, primeFactors);

		if (n<maxTrialFactor)
			return n;

		if (PrimeMath.isSquare(n)){
			final long x = PrimeMath.sqrt(n);
			if (x*x == n) {
				primeFactors.add(x);
				return x;
			}
		}
		// re-adjust the maximal factor we have to search for. If factors were found, which is quite
		// often the case for arbitrary numbers, this cuts down the runtime dramatically.
//		maxTrialFactor =  balanceTrial * Math.pow(n, ONE_THIRD);
		maxTrialFactor =  Math.pow(n, ONE_THIRD);
		final int kMax = (int) (Math.ceil(maxTrialFactor));
//		int kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
		final int multiplier = 4;
		final BigInteger n4 = BigInteger.valueOf(n * multiplier);
		final int multiplierSqrt = 2;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		// TODO is division by 4 more efficient if we use int?
		final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);

		// surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
		// for k=2 we know that x has to be even
		for (int k = 1; k <= kMax; k++) {
			final double sqrt4kn = multiplierSqrt * sqrt[k] * sqrtN;
			// adding a small constant to avoid rounding issues and rounding up is much slower then
			// using the downcast and adding a constant close to 1. Took the constant from the yafu code
			int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
			// use only multiplications instead of division here
			// TODO use a float here
			final double xRange = nPow1Sixth * sqrtInv[k];
			// since it is much bigger as SQRT (Long.MaxValue) we have to take a long
			// for k > kMax / 16 xRange is less then 1 unfortunately i can not use this fact to
			// speed up the runtime. Since the range for x is very small applying mod arguments to the x values
			// makes not much sense.
			final long xEnd = (long) (sqrt4kn + xRange);
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify. But reducing the numbers by a factor
			// 2 or 4 the runtime is only reduced by  something like 10%. This might be due to the Hotspot
			// Behavior of java?

			//			final int xStep = 1;
			int xStep = 2;
			if (k % 2 == 0) {
				xBegin |= 1;
			}
			//			else{
			//				xBegin &= (Integer.MAX_VALUE -1);
			//			}
			else{
				//				if (k*16 < kMax) {
				// TODO use the k%8 cases as well
				xStep = 4;
				xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
				//				}
			}
			long right = 0;
			if (xBegin <= xEnd) {
				BigInteger xBig = BigInteger.valueOf(xBegin);
				final BigInteger x2 = xBig.multiply(xBig);
				BigInteger kn4 = BigInteger.valueOf(k).multiply(n4);
				right = x2.subtract(kn4).longValue();
			}
			for(long x = xBegin; x <= xEnd; x+= xStep) {
				// in java the trick to replace the multiplication with an addition does not help
				// TODO use a less restrictive is square check and apply the error shift
				// instead of taking the square root (which is a very expensive operation)
				// and squaring it, we do some mod arguments to filter out non squares
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor != 1) {
						primeFactors.add(factor);
						return n / factor;
						// we know that the remaining factor has to be a prime factor
						// but this gives no speedup for ramdom numbers
						//						if (n != factor)
						//							factors.add(n / factor);
						//						return 1;
					}
				}
				right += x*xStep + xStep * xStep;
			}
		}
		return n;
	}
}
