package factoring.fermat.error;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;

import java.util.Collection;

/**
 * Created by Thilo Harich on 22.02.2018.
 */
public class ErrorShiftFact extends FermatFact{
    static boolean [] squaresMod = null;
    static int mod = 256;

    public static void initializeSquares() {
        if (squaresMod == null) {
            squaresMod = new boolean[mod];
            for (int i = 0; i < squaresMod.length; i++) {
                final int square = (i * i) % mod;
                squaresMod[square] = true;
            }
        }
    }

    public long findFactors(long n, Collection<Long> factors) {
        // This part is needed for the mod variants. It also improves performance
        for (int factor : smallFactors)
            if (n%factor == 0)
            {
                if (n/factor > 1)
                    factors.add(n/factor);
                return factor;
            }

        long sqrtN = (long) Math.ceil(Math.sqrt(n));
        long xEnd = (n / minFactor + minFactor) / 2;
        for (long x = sqrtN; x <= xEnd; x++) {
            long right = x*x - n;
            if (isProbableSquare(right)) {
                long y = PrimeMath.sqrt(right);
                if (y*y == right) {
                    long factorHigh = x + y;
                    factors.add(factorHigh);
                    return x - y;
                }
                else
                {
                    x = x + y;
                    right = x*x - n;
                    if (isProbableSquare(right)) {
                        y = PrimeMath.sqrt(right);
                        if (y * y == right) {
                            long factorHigh = x + y;
                            factors.add(factorHigh);
                            return x - y;
                        }
                    }
                }
            }
        }
        return n;
    }
    public static boolean isProbableSquare(long n) {
        initializeSquares();
        final int nMod = (int) (n & (mod-1));
        if (!squaresMod[nMod])
            return false;
        return true;
    }
}
