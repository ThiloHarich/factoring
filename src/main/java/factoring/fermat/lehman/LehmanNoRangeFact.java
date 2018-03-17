package factoring.fermat.lehman;

import java.util.Collection;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

/**
 *
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanNoRangeFact extends FermatFact {

	double balanceTrial = 1.5;

	@Override
	public long findFactors(long n, Collection<Long> factors) {
		final TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
		int maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
		smallFactoriser.setMaxFactor(maxTrialFactor);
		n = smallFactoriser.findPrimeFactors(n, factors);

		if (n<maxTrialFactor)
			return n;

		maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
		if (PrimeMath.isSquare(n)){
			final long x = PrimeMath.sqrt(n);
			if (x*x == n) {
				factors.add(x);
				return x;
			}
		}
		final int multiplier = 4;
		final int multiplierSqrt = (int) Math.sqrt(multiplier);
		final double balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
		final int kMax = (int) (maxTrialFactor / balanceTrialCube);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		final double xRange = (nPow2Third / multiplierSqrt) / Math.sqrt(n*multiplier);
		final double xMax = Math.sqrt(n)+1 + xRange;
		final double rightMax = xMax * xMax - n;
		for (int k = 1; k <= kMax; k++) {
			final long k4n = k * n * multiplier;
			final double sqrtKN = Math.sqrt(k4n);
			long x = (long) (Math.ceil(sqrtKN));
			long right;
			do {
				final long x2 = x * x;
				right = x2 - k4n;
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor != 1) {
						factors.add(factor);
						return n / factor;
					}
				}
				x++;
			}while (right < rightMax);
		}

		return n;
	}
}
