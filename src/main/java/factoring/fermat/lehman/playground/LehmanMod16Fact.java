package factoring.fermat.lehman.playground;

import java.util.Collection;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.variant.TrialInvFact;

/**
 * This is a version of the lehman factorizationByFactors, which is a variant of the fermat
 * factorizationByFactors.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster then the java version of the yafu lehman factorizationByFactors.
 * By storing the square roots of the multiplier k the range can be done faster.
 * It also uses a version of trial division, where the multiple inverse of the primes are stored.
 * So instead of a division a multiplication is needed to find out if a number is dividable
 * by a prime.
 * In the lehman algorithm in most of the cases i.e. k > n^1/3 / 16 the upper bound for x is less the
 * The lower bound plus 1. In this case at most one value of x has to be considered, but
 * the calculation of the lower and upper rage has to be done all the time.
 * Here we ignore the fact that we always increase x by 2 or 4.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 * <p>
 * Like in the YAFU implementation we get no speed when using smooth multipliers (for k) first.
 * This is surprising since most of the implementations use small multipliers, since they should
 * increase the chance that a created number is a square.
 * <p>
 * The Hart variant always just one x per multiplier k, this eliminates the determination of the
 * upper bound, but using it gives no extra speed. Again why?
 * <p>
 * We need 2^42/3 = 2^14 = 16364 locations to store the squares and the inversive. This does not fit in the L1 cache,
 * but it is accessed in a sequential way -> still efficient
 * <p>
 * Open questions, possible improvements :
 * - can we get rid of storing the square roots? how can we calculate them efficiently?
 * - In most of the cases (more then 93%) for k only one or none the value x^2 -n needs to be calculated.
 * <p>
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanMod16Fact implements FactorizationOfLongs {

    static double ONE_THIRD = 1.0 / 3;

    // This is a constant that is below 1 for rounding up double values to long
    protected static final double ROUND_UP_DOUBLE = 0.9999999665;

    double[] sqrt;
    float[] sqrtInv;
    byte[][] nextX;
    // a fast way to do the trial division phase
    final TrialInvFact smallFactoriser;
    int maxTrialFactor;

    /**
     * create a instance that can handle factorizations of numbers up to 2^bits
     * This uses memory 2^(bits/3)
     * more then 42 bits are not allowed since the precision of long does not allow more here
     *
     * @param bits
     */
    public LehmanMod16Fact(int bits) {
        // in contrast to the YAFU impl it makes no sense to do more trial divisions
        // and do not waste time in the slow lehman phase. This might be due to
        // the more performant impl of the lehman phase
        maxTrialFactor = (int) Math.ceil(Math.pow(1L << bits, ONE_THIRD));
        smallFactoriser = new TrialInvFact(maxTrialFactor);
        initSquares();
    }

    protected void initSquares() {
        // precalculate the square of all possible multipliers. This takes at most n^1/3
        final int kMax = maxTrialFactor;

        sqrt = new double[kMax + 10];
        sqrtInv = new float[kMax + 10];
        nextX = new byte[16][16];
        for (int i = 1; i < sqrt.length; i++) {
            final double sqrtI = Math.sqrt(i);
            sqrt[i] = sqrtI;
            // Since a multiplication is
            // faster then the division we also calculate the square root and the inverse
            sqrtInv[i] = (float) (1.0 / sqrtI);
        }
        System.out.printf(" sqrt table[0..%d] built: ", sqrt.length);
        for (int i = 0; i < Math.min(5, maxTrialFactor); i++) {
            System.out.printf("%f ", sqrt[i]);
        }
        System.out.printf(" ... %f %f\n", sqrt[sqrt.length - 2], sqrt[sqrt.length - 1]);

        final int mod = 16;
        final boolean[] squares = new boolean[16];
        for (int i = 0; i < mod; i++)
            squares[(i * i) % 16] = true;

        for (int kn = 0; kn < mod; kn++) {
            int lastX = 16;
            for (int x = mod - 1; x >= 0; x--) {
                final int right = PrimeMath.mod(x * x - kn, 16);
                nextX[kn][x] = (byte) (lastX - x - 1);
                if (squares[right]) {
                    lastX = x;
                }
            }
            nextX[kn][mod - 1] = (byte) (lastX);
        }
    }

    @Override
    public long findFactors(long nOrig, Collection<Long> primeFactors) {
        // with this implementation the lehman part is not slower then the trial division
        // we do not have to use a multiplier for the maximal factor were we apply the
        // trial division phase
        maxTrialFactor = (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        // factor out all small factors
        final long n = smallFactoriser.findFactors(nOrig, primeFactors);

        if (n < maxTrialFactor)
            return n;

        if (PrimeMath.isSquare(n)) {
            final long x = PrimeMath.sqrt(n);
            if (x * x == n) {
                primeFactors.add(x);
                return x;
            }
        }
        // re-adjust the maximal factor we have to search for. If factors were found, which is quite
        // often the case for arbitrary numbers, this cuts down the runtime dramatically.
        // TODO maybe we can do even better here?
        maxTrialFactor = (int) Math.pow(n, ONE_THIRD);
        final int kMax = maxTrialFactor;
        final long n4 = 4 * n;
        final double sqrtN = Math.sqrt(n);
        final int nMod16 = (int) (n % 16);
        // we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
        final int nPow2Third = maxTrialFactor * maxTrialFactor;
        final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
        // surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
        // for k=2 we know that x has to be even and results in a factor more often
        final double sqrt4N = 2 * sqrtN;
        for (int k = 1; k <= kMax; k++) {
            final double sqrt4kn = sqrt[k] * sqrt4N;
            // adding a small constant to avoid rounding issues and rounding up is much slower then
            // using the downcast and adding a constant close to 1. Took the constant from the yafu code
            int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
            // use only multiplications instead of division here
            final float xRange = nPow1Sixth * sqrtInv[k];
            final int xEnd = (int) (sqrt4kn + xRange);
            // instead of using a step 1 here we use the mod argument from lehman
            // to reduce the possible numbers to verify.
            // TODO do a mod argument here.

            final int knMod16 = k * nMod16 % 16;
            int xMod16 = xBegin % 16;
            xBegin += nextX[knMod16][xMod16];

            long x = xBegin;
            while (x <= xEnd) {
                // here we start using long values
                // in java the trick to replace the multiplication with an addition does not help
                final long x2 = x * x;
                final long right = x2 - k * n4;
                // instead of taking the square root (which is a very expensive operation)
                // and squaring it, we do some mod arguments to filter out non squares
                if (PrimeMath.isSquare(right)) {
                    final long y = (long) Math.sqrt(right);
                    final long factor = PrimeMath.gcd(n, x - y);
                    if (factor != 1) {
                        primeFactors.add(factor);
                        return n / factor;
                        // we know that the remaining factor has to be a prime factor
                        // but this gives no speedup for ramdom numbers
                        //						if (n != factor)
                        //							factors.add(n / factor);
                        //						return 1;
                    }
                }
                xMod16 = (int) (x + 1 % 16);
                x += nextX[knMod16][xMod16];
            }
        }
        return n;
    }
}
