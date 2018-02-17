package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.math.SquaresMod;

import java.util.Collection;

import static factoring.math.PrimeMath.mod;

/**
 * Created by Thilo Harich on 21.12.2017.
 *
 * Falls n < 24^2 kein multiplier
 * sonst
 * n = sqrt(n) * sqrt(n) * m   ;   m < sqrt(n)
 * f1 = (sqrt(n) + sqrt(n) * m) / 2
 * f1 = sqrt(n) * (m+1) /2
 *
 *  n = m * s  ,
 * n = s' * m   , s' >= s > m
 * f1 = (n/m + m)/2
 *
 * n = s' * s   , s' >= s
 * f1 = (n/s + s)/2
 *
 * n = s * s * s * m
 * f1 = (n/(s*m) +  s*m) / 2 =
 *
 *
 *
 */
public class Fermat24 extends FermatFact {

    // all factors below 24 to be sure the multiplier is the lowest factor
//    int [] smallFactors = {3,5,7,11,13,17,19,23};
    int [] smallFactors = {3};

    @Override
    public long findFactors(long n, Collection<Long> factors) {

        for (int factor : smallFactors)
        if (n%factor == 0)
        {
            if (n/factor > 1)
                factors.add(n/factor);
            return factor;
        }
        // now we can assume that we can always find an inverse of n mod 24
        int multiplier = 1;
        int step = 1;
        // Ensure that the lower factor is at least smallestFactor * multiplier ->
        if (n >= 24)
        {
            int nMod = mod(-n, 24);
            int multiplied = mod(nMod * n, 24);
            if (multiplied == 24 - 1 && nMod < 12) {
                multiplier = nMod;
                step = 12;
            }
        }
        long nMultiplied = multiplier*n;
        long factorLow = 5;
        long factorHigh = nMultiplied/factorLow;
        long xBegin = PrimeMath.sqrt(nMultiplied);
        xBegin = (xBegin / 24) * 24;
        // worst case
        long xEnd = (factorHigh + factorLow) / 2;
        for (long x = xBegin; x <= xEnd; x+=step) {
            long right = x * x - nMultiplied;
            if (SquaresMod.isSquare(right)) {
                long y = PrimeMath.sqrt(right);
                factorLow = x-y;
                // handle the multiplier
                if (factorLow % multiplier == 0 && factorLow > multiplier)
                    factorLow = factorLow/multiplier;
                if (n%factorLow == 0) {
                    long factor1 = n / factorLow;
                    if (factor1 != 1)
                        factors.add(factor1);
                    return factorLow;
                }
            }
        }
        return n;
    }
}