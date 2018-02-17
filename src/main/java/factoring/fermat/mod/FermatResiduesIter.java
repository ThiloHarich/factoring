package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.FermatResidueIter;

import java.util.Collection;

/**
 * It factorizes the number n, by
 *
 * For x^2 + 1 = y^2 mod 24 there is only two solution for x (0, 12). 24 = 2^3 * 3
 * We need a n' = n * m = -1 mod 24. m = -n mod 24 does this since:
 * 5*-5   =  -25 = -1 mod 24
 * 7*-7   =  -49 = -1 mod 24
 * 11*-11 = -121 = -1 mod 24
 *
 * Created by Thilo Harich on 29.11.2017.
 */
@Deprecated
public class FermatResiduesIter extends FermatFact {

    public FermatResiduesIter() {
    }

    public long findFactors(long n, Collection<Long> factors) {
        if (n < 200)
            return super.findFactors(n, factors);

//        FermatResiduesRec residues = FermatResiduesRec.create(n);
        FermatResidueIter residues = FermatResidueIter.create(n);

        return residues.findFactors(n, factors);
    }
}
