package factoring.math;

/**
 * Created by Thilo Harich on 23.11.2017.
 */
public class SquaresMod {

    private static final long MATH_SQRT_WORKS = 10354000L;
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

    /**
     * This method check if the value n is a square y^2.
     * All squares modulu mod2Pow are precomputed and stored.
     * Thus it can be checked easily if this number might be a prime
     * modulu mod2Pow, by looking up the square. Since mod2Pow is a power of 2
     * the modulo (%) operation is fast.
     *
     * @param n
     * @return
     */
    public static boolean isSquare(long n)
    {
        if (n < 0)
            return false;
        if (!isProbableSquare(n))
            return false;
        long sqrt = sqrt(n);
        return sqrt*sqrt == n;
    }

    public static boolean isProbableSquare(long n) {
        initializeSquares();
        final int nMod = (int) (n & (mod-1));
        if (!squaresMod[nMod])
            return false;
        return true;
    }

    public static long sqrt(long n) {
        long sqrt = (long)Math.sqrt(n);
        if (n < MATH_SQRT_WORKS)
        {
            return sqrt;
        }
        long sqrt2=sqrt, tmp;
        do{
            sqrt = sqrt2;

            final long nDivSqrt = (long) ((0.0 +n)/sqrt);
            sqrt2 = (sqrt + nDivSqrt)/2;
        }
        while (sqrt2 - sqrt > 0);
        return sqrt;
    }


}
