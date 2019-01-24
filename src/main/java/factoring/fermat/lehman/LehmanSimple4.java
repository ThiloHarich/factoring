
package factoring.fermat.lehman;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.trial.TrialInvFact;

/**
 * This Algorithm is based on the lehman algorithm and the Hart variant of it.
 * It somehow tries to combine both algorithms.
 * First it tries to find solution with a k = 2*3*3*5*7*k' = 630 * k'.
 * This ensures that the created a = ceil (sqrt (k*n)) | 1 produces numbers a^2 - kn which
 * are squares mod 2,9,5 and 7. This is there is a high chance that this numbers are squares and the
 * algorithm finds a factor.
 * If this part of the algorithm does not find a factor we investigate in numbers for an odd k.
 * Here we choose k=1. in each step we look into numbers of both kinds.
 *
 * We can not ensure that this algorithm always finds a factor.
 */
public class LehmanSimple4 extends FactorAlgorithmBase {
	private static final Logger LOG = Logger.getLogger(Lehman_TillSimple.class);

	/** This is a constant that is below 1 for rounding up double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final Gcd63 gcdEngine = new Gcd63();

	static TrialInvFact smallFact= new TrialInvFact((int) (1L << (48/3)));;

	private double[] sqrt;

	public LehmanSimple4() {
		initSqrts();
	}

	private void initSqrts() {
		// precompute sqrts for all possible k. Requires ~ (tDivLimitMultiplier*2^15) entries.
		final int kMax = (int) (40*Math.cbrt(1L<<50) + 1);
		sqrt = new double[kMax + 1];
		for (int i = 1; i < sqrt.length; i++) {
			sqrt[i] = Math.sqrt(i);
		}
	}


	@Override
	public String getName() {
		return "Lehman_TDivLast()";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	public long findSingleFactor(long N) {
		final int kLimit = (int)  Math.ceil(Math.cbrt(N));
		final long fourN = N<<2;
		final long four3N =3*fourN;
		final double sqrt4N = Math.sqrt(fourN);
		final long MN = fourN * 2*3*3*5*7;// = 630
		final double sqrtMN = Math.sqrt(MN);

		long aForK1 = (long) (sqrt4N  + ROUND_UP_DOUBLE);
		long aForK3 = (long) (Math.sqrt(four3N)  + ROUND_UP_DOUBLE);
		long aStep;
		// n % 4 == 3
		if ((N & 3) == 3) {
			aStep = 8;
			aForK1 += ((7 - N - aForK1) & 7);
			aForK3 += ((1 + N - aForK1) & 3);
		} else
		{
			aStep = 4;
			aForK1 += ((1 + N - aForK1) & 3);
			aForK3 += ((7 - N - aForK1) & 7);
		}
		final long aStepForK3 = 8 - aStep;
		// we have to increase kLimit here, because we only take good but seldom k's
		for (int k=1; k < 40*kLimit; aForK1 += aStep, aForK3 += aStepForK3) {
			long gcd = 1;
			for (int i=0; i < 8; i++) {
				// investigate in k = 0 mod multiplier only like the Hart variant does it
				long  a, test, b;
				// this is always the same code a loop/method (in java) is slower then copy the code here
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE) | 1;
				test = a*a - k++ * MN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = gcdEngine.gcd(a+b, N);
				}
				//---- optinally remove this code an increase the limit of the i.loop by a factor of 4
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE) | 1;
				test = a*a - k++ * MN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = gcdEngine.gcd(a+b, N);
				}
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE) | 1;
				test = a*a - k++ * MN;
				// This sqrt is fast, while the sqrt(k) of smaller size - but double precision? -  above is slow ???????
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = gcdEngine.gcd(a+b, N);
				}
				a = (long) (sqrtMN * sqrt[k] + ROUND_UP_DOUBLE) | 1;
				test = a*a - k++ * MN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					gcd = gcdEngine.gcd(a+b, N);
				}
				// --- end remove code
			}
			// Here k is 1 or 3 depending on N and we increase 'a' by aStep.
			// This phase should save us if the first phase does not find anything
			long test = aForK1*aForK1 - fourN;
			long b = (long) Math.sqrt(test);
			if (b*b == test) {
				gcd= gcdEngine.gcd(aForK1+b, N);
			}
			test = aForK3*aForK3 - four3N;
			b = (long) Math.sqrt(test);
			if (b*b == test) {
				gcd= gcdEngine.gcd(aForK1+b, N);
			}
			if (gcd > 1 && gcd < N)
				return gcd;
		}
		//
		// Check via trial division whether N has a nontrivial divisor d <= cbrt(N).
		return smallFact.findFactors(N, null);
	}
}
