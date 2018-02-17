package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.math.QuadraticDiophantineModBit;
import factoring.math.SquaresModArray;

import java.util.Collection;

/**
 * Since we look for solutions of the form:
 * xArray^2 - n = y^2
 * We look at solutions modulo some mod.
 * These solutions (for xArray) were calculated in a class {@link QuadraticDiophantineModBit}.
 * Then we only iterate over these solutions xArray.
 * The modulus will be extended during the sieving process.
 * Additionally we try to avoid the calculation for the sqrt to find out if xArray^2 - n is a square y^2.
 * This is done in the method {@link PrimeMath#isSquare(long)}.
 * Created by Thilo Harich on 05.08.2017.
 */
public class FermatModDynamic extends FermatFact {


    // we might drop some values and replace them with an exponentiation
//    private int[] modulus = {9,8};
    private int[] modulus = {3,4,5};
//    private int[] modulus = {9,8,7,25,11,13,17,19};

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        int modIndex = 0;
        int modCurrent = modulus[modIndex];
        int modProd = modCurrent;

        // first sieve with a step of 1, until we reach a starting point for sieving with the first mod
        long sqrtN = (long) Math.ceil(Math.sqrt(n));
        long sieveStart = ((sqrtN-1) / modCurrent) * modCurrent + modCurrent;
        long xBound = Math.min(2 + n / 6, sieveStart);
        for (long x = sqrtN; x <= xBound; x+=1) {
            long right = x * x - n;
            if (PrimeMath.isSquare(right)) {
                long y = PrimeMath.sqrt(right);
                factors.add((x + y));
                return n / (x + y);
            }
        }
        // for the first mod we calculate the sieve locations
        // 43 mod 60, xArray=301183, srtNMod = 8160
        int nMod = SquaresModArray.mod(-n, modCurrent);
//        int modIndex = 0;
        QuadraticDiophantineModBit squaresMod = new QuadraticDiophantineModBit(modulus[modIndex++]);
        int[] xCandidatesMod = squaresMod.xArray(nMod);
        long x = 0;
        // we sieve with the old sieve as long as we have just reached the new mod
        xBound = 2 + n / 6;
        while (x <= xBound) {
            int i = 0;
            int modNext = modIndex >= modulus.length ? Integer.MAX_VALUE : modulus[modIndex];
            while (xCandidatesMod[i] >= 0) {
                x = sieveStart + xCandidatesMod[i];
                int round = 0;
                // for the actual mod we sieve for rounds times, and then increase the modolus
                while (x <= xBound && round <= modNext) {
                    long right = x * x - n;
                    if (PrimeMath.isSquare(right)) {
                        long y = PrimeMath.sqrt(right);
                        factors.add((x + y));
                        return n / (x + y);
                    }
                    // we might integrate the merge step here, since we basically do the same
                    x += modCurrent;
                    round++;
                }
                i++;
            }
            sieveStart += modCurrent*modNext;
            // if we have a new mods left apply it
            if (modNext < Integer.MAX_VALUE) {
//                if (modIndex<modulus.length) {
                QuadraticDiophantineModBit squaresMods2 = new QuadraticDiophantineModBit(modNext);
                nMod = PrimeMath.mod(-n, modulus[modIndex]);
                squaresMods2.xMask(nMod);
                int[] xMerged = squaresMods2.merge(xCandidatesMod, squaresMod.xLength(), modCurrent);
                squaresMod = squaresMods2;
                xCandidatesMod = xMerged;
                modProd *= modNext;
                modIndex++;
            }
        }
        return n;
    }
}
