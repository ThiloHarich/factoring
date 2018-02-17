package factoring.fermat.residue;

import factoring.fermat.FermatFact;
import factoring.math.FermatResidueMergeByInversion;
import factoring.math.FermatResiduesRecursive;

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
public class FermatResiduesSieve extends FermatFact {

    public FermatResiduesSieve() {
    }

    public long findFactors(long n, Collection<Long> factors) {
        for (int prime : FermatResiduesRecursive.residuesArray) {
            if (n % prime == 0) {
                factors.add((long) prime);
                return n / prime;
            }
        }
        if(n < 200)
            return super.findFactors(n, factors);

        FermatResiduesRecursive residues = FermatResidueMergeByInversion.create(n);

        for (int i=0; residues.xArray[i] >= 0; i++) {
            long factor = residues.nextFermatResidue.findFactors(n, factors, residues.xArray[i], residues.mod);
            if (factor != n)
                return factor;
        }
        return n;
    }
}
