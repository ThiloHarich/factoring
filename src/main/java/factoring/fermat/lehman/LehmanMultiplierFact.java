package factoring.fermat.lehman;

import factoring.Factorizer;
import factoring.FindPrimeFact;
import factoring.fermat.FermatFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 11.01.2018.
 */
public class LehmanMultiplierFact extends FermatFact {
    @Override
    public long findFactors(long n, Collection<Long> factors) {
        FermatFact fact24 = new LehmanBlockFact(4,1);
        FindPrimeFact fact1 = new LehmanNoSqrtFact();
        long factor = fact24.findFactors(n, factors);
        if (factor != n)
            return factor;
        return fact1.findPrimeFactors(n, factors);
    }
}
