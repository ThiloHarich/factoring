package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * This implementation only calculates the expensive end of the for loop by the sqrt,
 * if the range is greater then the choosen block size.
 * After this boarder we will always consider block size values for the given level i.
 * This ensures that calculating the begin of the loop by sqrt(i*n) is not more time consuming
 * then calculating is x^2 - i*n is a square root.
 * Due to the fact that we hit all the lehman candidates we can be sure to stop at level which is
 * equal to the biggest number of the trial division phase.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanBlockFact extends FermatFact {

    int multiplier = 4;
    double blockSize = 1;
    int multiplierSqrt;

    public LehmanBlockFact(int multiplier, double blockSize) {
        this.multiplier = multiplier;
        this.blockSize = blockSize;
        multiplierSqrt = (int) Math.sqrt(multiplier);
        if (multiplierSqrt * multiplierSqrt != multiplier)
            throw new IllegalArgumentException("The multiplier has to be a square");
    }

    public long findFactors(long n, Collection<Long> factors) {
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int maxTrialFactor = (int) Math.ceil(Math.pow(n, 1.0/3));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<maxTrialFactor)
            return n;
        maxTrialFactor = (int) Math.ceil(Math.pow(n, 1.0/3));

        if (PrimeMath.isSquare(n)){
            long x = PrimeMath.sqrt(n);
            if (x*x == n) {
                factors.add(x);
                return x;
            }
        }
        double xRange = Double.MAX_VALUE;
        double nPow2Third = maxTrialFactor * maxTrialFactor;
        double nSqrt = Math.sqrt(n);
        for (int k = 1; k <= maxTrialFactor; k++) {
            long kn = k * n;
            long k4n = kn * multiplier;
            // we move the multiplier out of the calculation of the square root.
            // so we can reuse the result for calculating the end of the loop
            // unfortunately this means we are restricted to squares as multipliers
            // trick found https://github.com/DarkenCode/yafu/blob/master/factor/LehmanClean.c
            double kNSqrt = k * nSqrt;
            long xBegin = (long) (Math.ceil(multiplierSqrt * kNSqrt));
            // n^1/6 / (4 * sqrt (k)) = n^2/3 / (4 * sqrt(k*n))
            // only the sqrt for the begin of the loop is needed, and can be used for the end
            // of the loop as well. It looks a little bit like heron
            // here we still have the division which is an expensive operation. This is going to be avoided
            // in the other implementations.
            if (xRange > blockSize)
                xRange = (nPow2Third / 4) / kNSqrt;
            double xEnd = multiplierSqrt * kNSqrt + xRange;
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
