package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.variant.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class Hart24Fact extends FermatFact {

    TrialFactMod smallFactoriser = new TrialFactMod();

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
            int multiplier = 1;
            int mod = 1;
//            int nMod = mod(-n, mod);
//            int multiplied = mod(nMod * n, mod);
//            if (multiplied == mod - 1) {
//                multiplier = nMod;
//            }
            for (long i = 1; i <= 8*limit; i++) {
                long in = i * n;
                // TODO
                long x = PrimeMath.sqrt(in);
                if (x * x == n) {
                    factors.add(x);
                    return x;
                }
                x++;
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
                if (i < limit/4)
                {
                    in = i * n;
                    long sqrtIN = PrimeMath.sqrt(in);
                    long xBegin = ((sqrtIN + 2) /4) * 4;
                    for (x = xBegin; x <= xBegin+4; x += 1) {
                        // TODO analyze sqrt vs. division. We might use sqrt(i-i)
                        if (x * x == n) {
                            factors.add(x);
                            return x;
                        }
                        x2 = x * x;
                        right = x2 - in;
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
        }
        return n;
    }
}
