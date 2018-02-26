package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;
/**
 * This is a version of the lehman factorization, which is a variant of the fermat
 * factorization.
 * It is twice as fats as the java version of lehman factorization written by Warren D. Smith used in yafu.
 *
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we precalculate the square roots for the multipliers k and the inversion of it.
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 *
 * It first uses a factor of 1 and tries to find factors we have a range of n^1/6 here.
 * The height (kMax) of the calculation is much higher we can not expect much out of it
 */
public class LehmanNoSqrtFact2 extends FindPrimeFact {

	static double ONE_THIRD = 1.0/3;
	// to initialize the square roots we have to determine the maximal k
	static int K_MAX = (int) (Math.ceil(Math.pow(Long.MAX_VALUE, ONE_THIRD)));
	static float [] SQRT = new float[K_MAX + 1];
	static float [] SQRT_INV = new float[K_MAX + 1];
	// precalculate the square of all possible multipliers.
	static {
		for(int i = 1; i< SQRT.length; i++)
		{
			final float sqrtI = (float) Math.sqrt(i);
			SQRT[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the  root and the inverse
			SQRT_INV[i] = 1.0f / sqrtI;
		}
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
		// readjust the maximal factor we have to search for. If factors were found, which is quite
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
		// TODO use float avoid division
		final double nPow1Sixth = (nPow2Third / 4) / sqrtN;

		// first do the buttom level k=1 with a (8 times) bigger range TODO use a specific
		// algorithm here
		// then use the even k, which where every second number is a candidate
		// then use the odd k, which where every forth number is a candidate
		// This gives ~ 8% of speed
		int k = 1;
		double sqrt4kn = multiplierSqrt * SQRT[k] * sqrtN;
		// to avoid rounding issues we subtract a small number here.
		int xBegin = (int) (Math.ceil(sqrt4kn - 0.001));
		// use only multiplications instead of division here
		// TODO use a float here
		double xRange = nPow1Sixth * SQRT_INV[k];
		// take a bigger range here maximal 8*128 = 1024 we might use 9*8 = 72 or
		long xEnd = (long) (sqrt4kn + 4*xRange);
		// instead of using a step 1 here we use the mod argument from lehman
		// to reduce the possible numbers to verify. But reducing the numbers by a factor
		// 2 or 4 the runtime is only reduced by  something like 10%. This might be due to the Hotspot
		// Behavior of java?
		int xStep = 4;
		xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
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
		for (k = 2; k <= kMax; k += 2) {
			sqrt4kn = multiplierSqrt * SQRT[k] * sqrtN;
			// to avoid rounding issues we subtract a small number here.
			xBegin = (int) (Math.ceil(sqrt4kn - 0.001));
			// use only multiplications instead of division here
			// TODO use a float here
			xRange = nPow1Sixth * SQRT_INV[k];
			// since it is much bigger as SQRT (Long.MaxValue) we have to take a long
			// for k > kMax / 16 xRange is less then 1 unfortunately i can not use this fact to
			// speed up the runtime. Since the range for x is very small applying mod arguments to the x values
			// makes not much sense.
			xEnd = (long) (sqrt4kn + xRange);
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify. But reducing the numbers by a factor
			// 2 or 4 the runtime is only reduced by  something like 10%. This might be due to the Hotspot
			// Behavior of java?
			xStep = 2;
			xBegin |= 1;
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
		for (k = 3; k <= kMax; k += 2) {
			sqrt4kn = multiplierSqrt * SQRT[k] * sqrtN;
			// to avoid rounding issues we subtract a small number here.
			xBegin = (int) (Math.ceil(sqrt4kn - 0.001));
			// use only multiplications instead of division here
			// TODO use a float here
			xRange = nPow1Sixth * SQRT_INV[k];
			// since it is much bigger as SQRT (Long.MaxValue) we have to take a long
			// for k > kMax / 16 xRange is less then 1 unfortunately i can not use this fact to
			// speed up the runtime. Since the range for x is very small applying mod arguments to the x values
			// makes not much sense.
			xEnd = (long) (sqrt4kn + xRange);
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify. But reducing the numbers by a factor
			// 2 or 4 the runtime is only reduced by  something like 10%. This might be due to the Hotspot
			// Behavior of java?
			xStep = 4;
			xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
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
