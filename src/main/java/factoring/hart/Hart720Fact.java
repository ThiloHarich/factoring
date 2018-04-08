package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.variant.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class Hart720Fact extends FermatFact {

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
            int multiplier = 4*4 * 3*3 * 5*5; // = 4^2 * 3^2 * 5
//            multiplier = 1;
            for (long i = 1; i <= (3*limit)/2; i++) {
                long in = i * multiplier * n;
                if (in < 0)
                    System.out.println();
                long sqrtIN = PrimeMath.sqrt(in);
                for(long x = sqrtIN; x<= sqrtIN+1; x++) {
                    // TODO analyze sqrt vs. division. We might use sqrt(i-i)
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
