package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;
/**
 * This is a version of the lehman factorization, which is a variant of the fermat
 * factorization.
 * It is more twice as fats as the java version of the yafu lehman factorization.
 * In 15 out of 16 cases (for k) only one value of x has to be considered, but
 * the calculation of the rages has to be done all the time.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 * adressing a small array that fits in the level 1 cache
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanNoSqrtFact extends FindPrimeFact {

	static double ONE_THIRD = 1.0/3;
	// TODO see if it works for values above the precision/scale of the algorithm as well
	// If not check BigInteger
	// to initialize the square roots we have to determine the maximal k
	static int K_MAX = (int) (Math.ceil(Math.pow((4.0*Long.MAX_VALUE) /3.0, ONE_THIRD)));

	// if we know sqrt(k) sqrt(k+1) = sqrt(k * (k+1)/k) = sqrt(k) * sqrt((k+1)/k) = sqrt(k) * sqrt(1+1/k)
	// ~  sqrt(k) * (1+1/2k)= sqrt(k) + sqrt(k)^-1 / 2 so if we know sqrt(k) and sqrt(k)^-1 it is just an addition
	//  ~  sqrt(k) * (1+ 1/2k + 1/8 k-^2) =  sqrt(k) * (1+ 1/2k + (1/2 k)^2/2) -> we need a division and a multiplication
	// but better then calculating the sqrt
	// the error is < 1/8 * 1/k^2 < n^-2/3 , since k < n^1/3

	static float [] SQRT = new float[K_MAX + 1];
	static float [] SQRT_INV = new float[K_MAX + 1];
	// precalculate the square of all possible multipliers. This takes at most n^1/3
	static {
		for(int i = 1; i< SQRT.length; i++)
		{
			final float sqrtI = (float) Math.sqrt(i);
			SQRT[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the  root and the inverse
			SQRT_INV[i] = 1.0f / sqrtI;
		}
		System.out.printf(" sqrt table[0..%d] built: ", SQRT.length);
		for(int i=0; i<5; i++){ System.out.printf("%f ", SQRT[i]); }
		System.out.printf("%f ... %f %f\n", SQRT[5],SQRT[SQRT.length-2],SQRT[SQRT.length-1]);
	}

	@Override
	public long findPrimeFactors(long n, Collection<Long> factors) {
		final TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
		// with this implementation the lehman part is not slower then the trial division
		// we do not have to use a multiplier for the maximal factor were we apply the
		// trial division phase
		double maxTrialFactor =  Math.ceil(Math.pow(n, ONE_THIRD));
		smallFactoriser.setMaxFactor((int) maxTrialFactor);
		// factor out all small factors
		n = smallFactoriser.findPrimeFactors(n, factors);

		if (n<maxTrialFactor)
			return n;

		if (PrimeMath.isSquare(n)){
			final long x = PrimeMath.sqrt(n);
			if (x*x == n) {
				factors.add(x);
				return x;
			}
		}
		// re-adjust the maximal factor we have to search for. If factors were found, which is quite
		// often the case for arbitrary numbers, this cuts down the runtime dramatically.
		maxTrialFactor =  Math.pow(n, ONE_THIRD);
		final int kMax = (int) (Math.ceil(maxTrialFactor));
		final int multiplier = 4;
		final long n4 = n * multiplier;
		final int multiplierSqrt = 2;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		// TODO is division by 4 more efficient if we use int?
		final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);

		// surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
		// for k=2 we know that x has to be even
		for (int k = 1; k <= kMax; k++) {
			final double sqrt4kn = multiplierSqrt * SQRT[k] * sqrtN;
			// adding a small constant to avoid rounding issues and rounding up is much slower then
			// using the downcast and adding a constant close to 1. Took the constant from the yafu code
			int xBegin = (int) (sqrt4kn + 0.9999999665);
			// use only multiplications instead of division here
			// TODO use a float here
			final float xRange = nPow1Sixth * SQRT_INV[k];
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
			for(long x = xBegin; x <= xEnd; x+= xStep) {
				// in java the trick to replace the multiplication with an addition does not help
				final long x2 = x * x;
				final long right = x2 -  k * n4;
				// TODO use a less restrictive is square check and apply the error shift
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor != 1) {
						factors.add(factor);
						return n / factor;
						// we know that the remaining factor has to be a prime factor
						// but this gives no speedup for ramdom numbers
						//						if (n != factor)
						//							factors.add(n / factor);
						//						return 1;
					}
				}
			}
		}
		return n;
	}
}
