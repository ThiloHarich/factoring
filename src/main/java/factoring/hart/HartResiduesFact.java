package factoring.hart;

import factoring.fermat.FermatFact;
import factoring.fermat.residue.FermatResiduesArray;
import factoring.math.PrimeMath;
import factoring.trial.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class HartResiduesFact extends FermatFact {

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
            FermatResiduesArray resFact = new FermatResiduesArray(8);
            for (long i = 1; i <= 6 *limit/resFact.xRangeM1; i++) // xRange = 96
            {
                long in = i * n;
                long factor = resFact.findFactors(in, factors, n);
                    if (factor != in)
                        return factor;
                }
            }
        return n;
    }
}
