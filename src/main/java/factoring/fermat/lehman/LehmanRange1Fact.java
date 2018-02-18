package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanRange1Fact extends FindPrimeFact {

	static final float FACTOR_TRIAL_DIVISION = 1.3f;
	static final double ONE_THIRD = 1.0/3;

	static float BALANCE_TRIAL_CUBE = FACTOR_TRIAL_DIVISION * FACTOR_TRIAL_DIVISION;
	static int K_MAX = (int) (Math.ceil(Math.pow(Long.MAX_VALUE, ONE_THIRD)) / BALANCE_TRIAL_CUBE);

	static float [] SQRT = new float[K_MAX + 1];
	static float [] SQRT_INV = new float[K_MAX + 1];
	static {
		for(int i=1; i<SQRT.length; i++)
		{
			final float sqrtI = (float) Math.sqrt(i);
			SQRT[i] = sqrtI;
			SQRT_INV[i] = 1.0f / sqrtI;
		}
	}
	static final int X_BEGIN_MOD_4 = Integer.MAX_VALUE - 3;

	// we store all the variables we need as fields to make the methods easy
	private Collection<Long> factors;
	private long n;
	int multiplier;
	private int multiplierSqrt;
	private long n4;
	private int nMod4;
	private double nPow1Sixth;
	private double sqrtN;
	int k = 1;
	int kMax = K_MAX;

	@Override
	public long findPrimeFactors(long nIn, Collection<Long> factorsIn) {
		factors = factorsIn;
		n = nIn;
		final TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
		double maxTrialFactor =  Math.ceil(FACTOR_TRIAL_DIVISION * Math.pow(n, ONE_THIRD));
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
		// readjust the maximal factor
		// TODO we can approximate the max factor here
		maxTrialFactor =  FACTOR_TRIAL_DIVISION * Math.pow(n, ONE_THIRD);
		kMax = (int) (Math.ceil(maxTrialFactor / BALANCE_TRIAL_CUBE));
		//		final int k4Range1 = (int) (FACTOR_TRIAL_DIVISION * FACTOR_TRIAL_DIVISION * FACTOR_TRIAL_DIVISION * kMax / 16);
		multiplier = 4;
		n4 = n * multiplier;
		multiplierSqrt = 2;
		sqrtN = Math.sqrt(n);
		nMod4 = (int) (n % 4);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		// TODO is division by 4 more efficient if we use int?
		// TODO use float avoid division
		nPow1Sixth = (nPow2Third / 4) / sqrtN;
		//		final int nMod3 = (int) (n % 3);

		//        x^2 - 4kn = y^2 mod 3
		//        x^2 - kn = y^2 mod 3
		//        x^2 = y^2 mod 3 , k*n == 0 -> k = 0 mod 3 since n != 0 mod 3 -> all solutions
		k = 1;
		long factor = getFactor();
		if (factor > 0)
			return factor;
		factor = getFactorOneX();
		if (factor > 0)
			return factor;
		//        long factor = getFactor(3, kMax, 3);
		//        if (factor > 0) return factor;
		//        x^2 - 1 = y^2 mod 3 , k*n == 1 -> x=1,2 mod 3 -> k = n^-1 mod 3 -> k = n mod 3
		//        mod 9 there are only 2 possible solutions as well
		//        factor = getFactor(nMod3, kMax, 3, 1, 3);
		//        if (factor > 0) return factor;
		//        factor = getFactor(nMod3, kMax, 3, 2, 3);
		//        if (factor > 0) return factor;
		//        x^2 - 2 = y^2 mod 3 , k*n == 2 -> x=0 mod 3   -> k = 2* n^-1 mod 3  -> k = 2n mod 3
		//        factor = getFactor(3 - nMod3, kMax, 3, 0, 3);
		//        if (factor > 0) return factor;

		return n;
	}

	public long getFactor() {
		double xRange;
		do {
			final long kn4 = k * n4;
			final double sqrtKn = SQRT[k] * sqrtN;
			final double sqrt4kn = multiplierSqrt * sqrtKn;
			int xBegin = (int) (Math.ceil(sqrt4kn - 0.001));
			// use only multiplications instead of division here
			xRange = nPow1Sixth * SQRT_INV[k];
			final long xEnd = (long) (sqrt4kn + xRange);
			int xStep;
			if (k % 2 == 0) {
				xStep = 2;
				xBegin |= 1;
			}
			else{
				xStep = 4;
				xBegin &= X_BEGIN_MOD_4;
				//                xBegin -= xBegin % 4;
				xBegin |= (k + nMod4) % 4;
				//				if (xBegin < sqrtKn)
				//					xBegin += 4;
			}
			for(long x = xBegin; x <= xEnd; x += xStep) {
				final long x2 = x * x;
				final long right = x2 - kn4;
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor > 1) {
						factors.add(factor);
						if (n != factor)
							factors.add(n / factor);
						return 1;
					}
				}
			}
			k++;
		}
		while(xRange > 1);
		return -1;
	}
	public long getFactorOneX() {
		for (; k <= kMax; k++) {
			final double sqrtKn = SQRT[k] * sqrtN;
			final double sqrt4kn = multiplierSqrt * sqrtKn;
			final long x = (long) (Math.ceil(sqrt4kn - 0.001));
			// check x mod 2 and 4
			//			final int kMod2 = k % 2;
			if ((x&1) != (k & 1))
				//				if ((x%2==0) ^ (kMod2 ==0))
			{
				//                if (kMod2==0 ||  (k + nMod4) % 4 == (x % 4) )
				{
					final long x2 = x * x;
					final long kn4 = k * n4;
					final long right = x2 - kn4;
					if (PrimeMath.isSquare(right)) {
						final long y = (long) Math.sqrt(right);
						final long factor = PrimeMath.gcd(n, x - y);
						if (factor != 1) {
							factors.add(factor);
							if (n != factor)
								factors.add(n / factor);
							return 1;
						}
					}
				}
			}
		}
		return -1;
	}

	//	public int nextX(int xStep, int nMod4, int k, int x) {
	//		if (k % 2 == 0) {
	//			if (x % 2 == 0)
	//				x += xStep;
	//		}
	//		else{
	//			// Instead of doing the while loop we may use the inverse
	//			final int knMod = (k + nMod4) % 4;
	//			while ((x % 4) != knMod)
	//				x += xStep;
	//		}
	//		return x;
	//	}
}
