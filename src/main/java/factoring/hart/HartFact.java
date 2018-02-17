package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class HartFact extends FermatFact {

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
            if (PrimeMath.isSquare(n)){
                long x = PrimeMath.sqrt(n);
                if (x*x == n) {
                    factors.add(x);
                    return x;
                }
            }
            for (long i = 1; i <= 4*limit; i++) {
                long in = i * n;
                // TODO
                long x = PrimeMath.sqrt(in);
                long x2 = x*x;
                long right = x2 - in;
                if (right < 0)
                    System.out.println();
                if (PrimeMath.isSquare(right)) {
                    long y = (long) Math.sqrt(right);
                    long factor = PrimeMath.gcd(n, x-y);
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
