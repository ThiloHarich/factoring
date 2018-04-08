package factoring.fermat.lehman.playground;

import factoring.FactorizationOfLongs;
import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.playgound.TrialWithPrimesFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanSquaresFact  implements FactorizationOfLongs {

    double balanceTrial = 1.9;
    int multiplier = 4;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        n = smallFactoriser.findFactors(n, factors);

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
        int multiplierSqrt = (int) Math.sqrt(multiplier);
        double balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
        // effectively the outer loop (height) is cut down by a balanceTrial^2
        // and the inner loop (width) is multiplied by balanceTrial^2
        // by doing this, it is much more important to reduce the candidates of each level.
        int kMax = (int) (maxTrialFactor / balanceTrialCube);
        double nPow2Third = maxTrialFactor * maxTrialFactor;
//        int K_MAX = (int)Math.pow(n, 1.0/3);
        int step;
        int nMod4 = (int) (n % 4);
        final double v = nPow2Third / multiplierSqrt;
        long factor = getFactor(1, kMax, 1, n, factors, nMod4, v);
        if (factor > 0) return factor;
//        factor = getFactor(1, K_MAX, 3, n, factors, nMod4, v);
//        if (factor > 0) return factor;
//        factor = getFactor(2, K_MAX, 3, n, factors, nMod4, v);
//        if (factor > 0) return factor;
        return n;
    }

//    1, 4, 9, 16, 25, ...
//    2, 8,18, 32, 50
//    3,
    public long getFactor(int kBegin, int kEnd, int kStep, long n, Collection<Long> factors, int nMod4, double v) {
//        int xStep;
        for (int k = kBegin; k <= kEnd; k += kStep) {
            long k4n = k * n * multiplier;
            double sqrtKN = Math.sqrt(k4n);
            double xRange = v / sqrtKN;

            for (int i = 1, k2 = k; k2 <= kEnd; k2 = k*i*i, i++) {
                sqrtKN *= i;
                xRange /= i;
                // if k is a square we can reuse the above steps
                long xBegin = (long) (Math.ceil(sqrtKN));
                int xEnd = (int) (sqrtKN + xRange);
                if (k2 % 2 == 0) {
//                xStep = 2;
                    xBegin |= 1;
                } else {
                    // TODO define the case
//                if (k%4 == 1)
//                {
//                    xStep = 8;
//                    xBegin -= xBegin % 4;
//                    xBegin += ((k+3)%8 + nMod4) % 4;
//                }
//                else
                    {
//                    xStep = 4;
                        xBegin -= xBegin % 4;
                        xBegin += (k2 + nMod4) % 4;
                        if (xBegin < sqrtKN)
                            xBegin += 4;
                    }
                }
                long x = xBegin;
                while (x <= xEnd) {
                    long x2 = x * x;
                    long right = x2 - k4n;
                    if (PrimeMath.isSquare(right)) {
                        long y = (long) Math.sqrt(right);
                        long factor = PrimeMath.gcd(n, x - y);
                        if (factor != 1) {
                            factors.add(factor);
                            if (n / factor != 1)
                                factors.add(n / factor);
                            return 1;
                        }
                    }
                    if (k2 % 2 == 0) {
                        x += 2;
                    } else {
                        x += 4;
                    }
                }
            }
        }
        return -1;
    }
}
