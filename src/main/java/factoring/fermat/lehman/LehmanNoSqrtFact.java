package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;
/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanNoSqrtFact extends FindPrimeFact {

	static float BALANCE_TRIAL = 1.3f;
	static double ONE_THIRD = 1.0/3;
	final static int xBeginMod4 = Integer.MAX_VALUE - 3;

	static float balanceTrialCube = BALANCE_TRIAL * BALANCE_TRIAL;
	static int kMax = (int) (Math.ceil(Math.pow(Long.MAX_VALUE, ONE_THIRD)) / balanceTrialCube);

	static float [] sqrt = new float[kMax + 1];
	static float [] sqrtInv = new float[kMax + 1];
	static {
		for(int i=1; i<sqrt.length; i++)
		{
			final float sqrtI = (float) Math.sqrt(i);
			sqrt[i] = sqrtI;
			sqrtInv[i] = 1.0f / sqrtI;
		}
	}

	@Override
	public long findPrimeFactors(long n, Collection<Long> factors) {
		final TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
		double maxTrialFactor =  Math.ceil(BALANCE_TRIAL * Math.pow(n, ONE_THIRD));
		smallFactoriser.setMaxFactor((int) maxTrialFactor);
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
		maxTrialFactor =  BALANCE_TRIAL * Math.pow(n, ONE_THIRD);
		kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
		final int multiplier = 4;
		final long n4 = n * multiplier;
		final int multiplierSqrt = 2;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		// TODO is division by 4 more efficient if we use int?
		// TODO use float avoid division
		final double nPow1Sixth = (nPow2Third / 4) / sqrtN;

		for (int k = 1; k <= kMax; k++) {
			//            long kn4 = k * n4;
			final double sqrt4kn = multiplierSqrt * sqrt[k] * sqrtN;
			int xBegin = (int) (Math.ceil(sqrt4kn - 0.001));
			// use only multiplications instead of division here
			final double xRange = nPow1Sixth * sqrtInv[k];
			// since it is much bigger as sqrt (Long.MaxValue) we have to take a long
			//            double xEnd = sqrt4kn + xRange;
			final long xEnd = (long) (sqrt4kn + xRange);
			int xStep;
			if (k % 2 == 0) {
				xStep = 2;
				xBegin |= 1;
			}
			else{
				xStep = 4;
				xBegin &= xBeginMod4;
				//                xBegin -= xBegin % 4;
				xBegin |= (k + nMod4) % 4;
				//                if (xBegin < sqrtKn)
				//                    xBegin += 4;
			}
			for(long x = xBegin; x <= xEnd; x+= xStep) {
				final long x2 = x * x;
				final long right = x2 -  k * n4;
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor != 1) {
						factors.add(factor);
						return n / factor;
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
