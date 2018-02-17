package factoring.fermat.residue;

import factoring.fermat.FermatFact;
import factoring.trial.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class TrialResiduesFact extends FermatFact {

    TrialFactMod smallFactoriser = new TrialFactMod();
    int factor = 1;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        // TODO estimate the speedup should be something like
        // exp ( ln(n')/ ln(ln(n')))
        int estSpeedup = 10;
        int limit =  (int) Math.ceil(Math.pow(n, .5))/estSpeedup;
        limit = Math.max(19, limit);
        smallFactoriser.setLimit(limit);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<= limit) {
            n = smallFactoriser.findPrimeFactors(n, factors);
        }
        else
        {
//            if (n < 1600)
            FermatResiduesRec fermatRes = new FermatResiduesRec(limit);
//            long factor = fermatRes.findPrimeFactors(n, factors);
//            FermatFact fermatRes = new FermatFact(limit);
            return fermatRes.findFactors(n, factors, n);
        }
        return n;
    }
}
