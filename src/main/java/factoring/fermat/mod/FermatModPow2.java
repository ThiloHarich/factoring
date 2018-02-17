package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.*;

import java.util.Collection;

/**
 * Since we look for solutions of the form:
 * xArray^2 - n = y^2
 * we Define the set squaresMask(m) of all pairs (l,r) where r = l^2 mod m
 * squaresMask(m) - n is defined as (l,r) where r = l^2 -n
 * We need disjunction for the second element of
 * squaresMask(m) - n and squaresMask(m). Candidates for xArray
 * are the first elements of the pairs above
 * so we look for
 *
 * (xArray+i*m)^2 = xArray^2 + 2*i*m +1 = xArray^2 + 1 mod 2*i*m
 * Created by Thilo Harich on 05.08.2017.
 */
public class FermatModPow2 extends FermatFact {

    public long findFactors(long n, Collection<Long> factors) {
        if(n < 200)
            return super.findFactors(n, factors);

        long xBound = 2 + n / 6;
//        int mod = QuadraticDiophantineModBig.getMod(xBound);
        int mod = 9*5*16;

        // TODO
        // TODO store the different squares for the residueClasses
//        13 * 17^2 * 101 * 173
//        int mod = Math.max(4, 1 << PrimeMath.log2(xBound/64));
//        int mod = 128;

//        QuadraticDiophantineModBit squares = new QuadraticDiophantineModBit(mod);
        QuadraticDiophantineModBit squares = QuadraticDiophantineModBig.create(mod);
        int[] xCandidatesMod = squares.xArray((int) -n);
        long sqrtN = (long) Math.ceil(Math.sqrt(n));
        long sqrtNMod = (sqrtN / mod) * mod;
        int i = 0;
        while (xCandidatesMod[i] >= 0) {
            long x = sqrtNMod + xCandidatesMod[i];
            while (x <= xBound) {
                long right = x * x - n;
                if (PrimeMath.isSquare(right)) {
                    long y = PrimeMath.sqrt(right);
                    factors.add((x + y));
                    return x - y;
                }
                x+= mod;
            }
            i++;
        }
        return n;
    }
}
