package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanOneSqrtFact extends FermatFact {

    double balanceTrial = 1.5;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<maxTrialFactor)
            return n;

        maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        if (PrimeMath.isSquare(n)){
            long x = PrimeMath.sqrt(n);
            if (x*x == n) {
                factors.add(x);
                return x;
            }
        }
        int multiplier = 4;
        int multiplierSqrt = (int) Math.sqrt(multiplier);
        double balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
        // effectively the outer loop (height) is cut down by a balanceTrial^2
        // and the inner loop (width) is multiplied by balanceTrial^2
        // by doing this, it is much more important to reduce the candidates of each level.
        int kMax = (int) (maxTrialFactor / balanceTrialCube);
        double nPow2Third = maxTrialFactor * maxTrialFactor;
//        int kMax = (int)Math.pow(n, 1.0/3);
        for (int k = 1; k <= kMax; k++) {
            long k4n = k * n * multiplier;
            // we move the multiplier out of the calculation of the square root.
            // so we can reuse the result for calculating the end of the loop
            // unfortunately this means we are restricted to squares as multipliers
            // trick found https://github.com/DarkenCode/yafu/blob/master/factor/LehmanClean.c
            double sqrtKN = Math.sqrt(k4n);
            long xBegin = (long) (Math.ceil(sqrtKN));
            // instead of directly calculating x range = 1/(4*(r+1))* sqrt (n/k) = n^1/6 / (4 * sqrt (k))
            // we use lowest factor = sqrt(n/r+1)
            // ->  x range = (lowest factor)^2 /  sqrt(8k*n)
            // = n/((r+1) * sqrt(4k*n)) = sqrt(n)/((r+1) * 4 * sqrt(k))
            // (4 * sqrt(k*n)), which has the same value so
            // only the sqrt for the begin of the loop is needed, and can be used for the end
            // of the loop as well. Surprisingly this gives a speedup of ~ 2.5
            // This means that calculating the square root of a small number k is much more
            // time consuming then the calculation the square root of the big numbers k*n do not get it.
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
