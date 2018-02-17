package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanFact extends FermatFact {

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        //
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int limit = (int) Math.pow(n, 1.0/3);
        smallFactoriser.setMaxFactor(limit);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<= limit) {
            n = smallFactoriser.findPrimeFactors(n, factors);
        }
        else
        {
            if (PrimeMath.isSquare(n)){
                long x = PrimeMath.sqrt(n);
                if (x*x == n) {
                    factors.add(x);
                    return x;
                }
            }
            long multiplier = 4;
            double nPow6Div4 = Math.pow(n, 1.0/6) / 4;
            for (long k = 1; k <= limit; k++) {
                long k4n = multiplier * k * n;
                double sqrtKN = Math.sqrt(k4n);
                long xBegin = (long) Math.ceil(sqrtKN);
                // TODO only calculate the range when it might changes
                double xRange = nPow6Div4 / Math.sqrt(k);
                for(long x = xBegin; x <= sqrtKN + xRange; x++) {
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
        }
        return n;
    }
}
