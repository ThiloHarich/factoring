package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanNoSqrt4Fact extends FermatFact {


    static double [] sqrt = new double[65536];
    static {
        for(int i=1; i<sqrt.length; i++)
        {
            sqrt[i] = Math.sqrt(4*i);
        }
    }
    double balanceTrial = 1.5;
    double firstTrialPortion = .1;
    @Override
    public long findFactors(long n, Collection<Long> factors) {
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<maxTrialFactor)
            return n;

        // readjust the maximal factor
        maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        if (PrimeMath.isSquare(n)){
            long x = PrimeMath.sqrt(n);
            if (x*x == n) {
                factors.add(x);
                return x;
            }
        }
        int multiplier = 4;
        int multiplierSqrt = 2;
        double nPow2Third = maxTrialFactor * maxTrialFactor;
        double balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
        int kMax = (int) (maxTrialFactor / balanceTrialCube);
        double sqrtN = Math.sqrt(n);

        for (int k = 1; k <= kMax; k++) {
            long kn = k * n;
            long k4n = kn * multiplier;
            // here we do not need to calculate the sqrt for k anymore since we have
            // it saves ~ 15% speed, we might store it only for the lower numbers and
            double sqrtKN = sqrt[k] * sqrtN;
            long xBegin = (long) (Math.ceil(sqrtKN));
            double xRange = (nPow2Third / multiplierSqrt) / sqrtKN;
            double xEnd = sqrtKN + xRange;
            for(long x = xBegin; x <= xEnd; x++) {
                long x2 = x * x;
                long right = x2 - k4n;
                if (PrimeMath.isSquare(right)) {
                    long y = (long) Math.sqrt(right);
                    long factor = PrimeMath.gcd(n, x - y);
                    if (factor != 1) {
                        factors.add(factor);
                        return n / factor;
                    }
                }
            }
        }

        return n;
    }
}
