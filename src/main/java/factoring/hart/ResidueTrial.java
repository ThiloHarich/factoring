package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class ResidueTrial extends FermatFact {

    TrialFactMod smallFactoriser = new TrialFactMod();
    int factor = 1;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        int limit =  (int) Math.ceil(Math.pow(n, .33));
        smallFactoriser.setLimit(limit);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<= limit) {
            n = smallFactoriser.findPrimeFactors(n, factors);
        }
        else
        {

            for (long i = 1; i <= limit*limit; i++) {
                long in = i * n;
                long x = PrimeMath.sqrt(in) + 1;
                long x2 = x*x;
//                if (xArray < in)
//                    xArray = (xArray+1)*(xArray+1);
                long right = x2 - in;
                if (PrimeMath.isSquare(right)) {
                    long y = (long) Math.sqrt(right);
                    long factor = PrimeMath.gcd(n, x-y);
                    factors.add(factor);
                    n = n / (x-y);
                    return n;
                }
            }
        }
        return n;
    }
}
