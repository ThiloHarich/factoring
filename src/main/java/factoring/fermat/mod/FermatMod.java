package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.*;

import java.util.Collection;

/**
 * Since we look for solutions of the form:
 * xArray^2 - n = y^2
 * We look at solutions modulo some mod.
 * These solutions (for xArray) were calculated in a class {@link QuadraticDiophantineModBit}.
 * Then we only iterate over these solutions xArray.
 * Additionally we try to avoid the calculation for the sqrt to find out if xArray^2 - n is a square y^2.
 * This is done in the method {@link PrimeMath#isSquare(long)}.
 * Created by Thilo Harich on 05.08.2017.
 */
public class FermatMod extends FermatFact {

    // TODO mod should be calculated depending on n
    // we might add factor after factor
//    int mod = 9*8*7*5;
//    int mod = 3*4*5;
    byte mod = 27;
QuadraticDiophantineModBit squares = new QuadraticDiophantineModBit(mod);
//    QuadraticDiophantineModBit squares = new QuadraticDiophantineModBit(3,4,5);
//    QuadraticDiophantineModBit squares = new QuadraticDiophantineModBit(mod);


    @Override
    public long findFactors(long n, Collection<Long> factors) {
        if(n < 200)
            return super.findFactors(n, factors);

        long xBound = 2 + n / 6;
        byte nMod = (byte) PrimeMath.mod(-n, mod);
        int[] xCandidatesMod = squares.xArray(nMod);
        int mod = squares.getMod();
        long sqrtN = (long) Math.ceil(Math.sqrt(n));
        long sqrtNMod = (sqrtN / mod) * mod;
        int i = 0;
        while (xCandidatesMod[i] >= 0) {
            long x = sqrtNMod + xCandidatesMod[i];
            while (x <= xBound) {
                long right = x * x - n;
                if (SquaresMod.isSquare(right)) {
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
