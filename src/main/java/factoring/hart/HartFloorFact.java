package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.variant.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class HartFloorFact extends FermatFact {

    TrialFactMod smallFactoriser = new TrialFactMod();
    int factor = 1;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        int limit =  (int) Math.ceil(Math.pow(n, .33));
        smallFactoriser.setLimit(limit);
        n = smallFactoriser.findFactors(n, factors);

        if (n<= limit) {
            n = smallFactoriser.findFactors(n, factors);
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
            int range = 8;
            int multiplier = 24;
            for (long i = 1; i <= 3 *limit/range; i++) {
                long in = i * multiplier * n;
                long sqrtIN = PrimeMath.sqrt(in);
                long xBegin = sqrtIN;
                for(long x = xBegin; x< xBegin + range; x++) {
                    long x2 = x * x;
                    long right = x2 - in;
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
