package factoring.fermat.lehman.playground;

import factoring.FactorizationOfLongs;
import factoring.fermat.lehman.LehmanYafuFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

import java.util.Collection;

/**
 * Here we first look at k's with many factor 3 and 5
 * Here is a order we use Mod 15:
 * 15,  9,  3, 12, 6,
 * 10,  4, 13, 7,  1,
 *  5, 14,  8, 2, 11,
 *
 *
 *
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanFactorFinderMod60 implements FactorizationOfLongs {

    static double ONE_THIRD = 1.0 / 3;
    // to be fast to decide if a number is a square we consider the
    // number modulo some values. First filer/modulus is a power of 2
    // since we can easily calculate the remainder by using bit and operation
    // here also low mod like 64 work fine as well
    static int mod_1024 = 1024;
    // using an additional second filter gives around 10% of speedup
    // Since we do not want to exceed the running time of  then n^1/3 ~ 20000 we
    // restrict the mod to the product of small primes below these limit
    static int mod_9_5_7_11 = 3 * 3 * 5 * 7 * 11; // 3465
    // we might extend this limit if we know we have to factorize more number we might extend
    //	static int mod_9_5_7_11 = 3*3 * 5 * 7 * 11 * 13; // 45045

    static int [] invMod9 = new int[9];


    private static boolean[] isSquareMod_1024;
    // for the other mod arguments we have to do a real expensive "%" operation
    // so we might have time to do a bit more expensive but space saving representation
    // of the squares
    private static boolean[] isSquareMod_9_5_7_11;

    float maxTrialFactorMultiplier = 1.0f;
    float maxFactorMultiplierCube;

    // This is a constant that is below 1 for rounding up double values to long
    protected static final double ROUND_UP_DOUBLE = 0.9999999665;

    double[] sqrt;
    float[] sqrtInv;
    // a fast way to do the trial division phase
    final TrialInvFact smallFactoriser;
    int maxTrialFactor;

    static {
        for (int i=1; i<9; i++){
            invMod9[i] = (int) PrimeMath.invert(i, 9);
        }
        isSquareMod_1024 = new boolean[mod_1024];
        //		isSquareMod9_5_7_11 = new BitSet(mod9_5_7_11);
        isSquareMod_9_5_7_11 = new boolean[mod_9_5_7_11];
        for (int i = 0; i < mod_9_5_7_11; i++) {
            if (i < mod_1024)
                //				isSquareMod_1024[(i*i) % mod_1024] = true;
                isSquareMod_1024[(i * i) & 1023] = true;
            //			isSquareMod_9_5_7_11[(i*i) % mod_9_5_7_11] = true;
            isSquareMod_9_5_7_11[(i * i) % 3456] = true;
        }
        System.out.println("sqrt mod table built                       bytes used : " + (isSquareMod_1024.length + isSquareMod_9_5_7_11.length));
    }


    /**
     * create a instance that can handle factorizations of numbers up to 2^bits
     * This uses memory 2^(bits/3)
     * more then 42 bits are not allowed since the precision of long does not allow more here
     *
     * @param bits                  the maximum number of bits the algorithm is going to work.
     * @param maxTrialFactorMultiplier Defines the size of the factors the algorithm is first looking for.
     *                              The algorithm is working best for number where the second highest factor is
     *                              maxTrialFactorMultiplier * n^1/3. If maxTrialFactorMultiplier is 1 or lower we first do trial division
     *                              for numbers below n^1/3 then do the lehman factorizationByFactors above n^1/3.
     *                              For maxTrialFactorMultiplier >= 1 we first inspect numbers above maxTrialFactorMultiplier * n^1/3 with the lehman
     *                              factorization. After completing this phase we do
     *                              trial division for numbers below maxTrialFactorMultiplier * n^1/3<br>
     *                              Use<br>
     *                              maxTrialFactorMultiplier <= 1<br>
     *                              - if you do not know anything about the numbers to factorize<br>
     *                              - if you know the second highest factor can not exceed n^maxTrialFactorMultiplier < n^1/3<br>
     *                              here all found factors are prime factors and will be stored in the primeFactors Collection.
     *                              maxTrialFactorMultiplier = 3.2 (this is the optimal value) <br>
     *                              - if you know for most of the numbers the maximal factors will exceed 3.2*n^1/3<br>
     *                              In the last case {@link #findFactors(long, Collection)} might return a composite number
     */
    public LehmanFactorFinderMod60(int bits, float maxTrialFactorMultiplier) {
        if (bits > 41)
            throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
        this.maxTrialFactorMultiplier = maxTrialFactorMultiplier < 1 ? 1 : maxTrialFactorMultiplier;
        maxTrialFactor = (int) Math.ceil(this.maxTrialFactorMultiplier * Math.pow(1L << bits, ONE_THIRD));
        maxFactorMultiplierCube = this.maxTrialFactorMultiplier * this.maxTrialFactorMultiplier * this.maxTrialFactorMultiplier;
        smallFactoriser = new TrialInvFact(maxTrialFactor);
        initSquares();
    }

    protected void initSquares() {
        // precalculate the square of all possible multipliers. This takes at most n^1/3
        final int kMax = (int) (maxTrialFactor / maxFactorMultiplierCube);

        sqrt = new double[kMax + 10];
        sqrtInv = new float[kMax + 10];
        for (int i = 1; i < sqrt.length; i++) {
            final double sqrtI = Math.sqrt(i);
            sqrt[i] = sqrtI;
            // Since a multiplication is
            // faster then the division we also calculate the square root and the inverse
            sqrtInv[i] = (float) (1.0 / sqrtI);
        }
        System.out.println("sqrt tables for max trial factor " + maxTrialFactorMultiplier + " built bytes used : " + sqrt.length);
    }

    @Override
    public long findFactors(long n, Collection<Long> primeFactors) {
        if (n > 1l << 41)
            throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
        // with this implementation the lehman part is not slower then the trial division
        // we only apply the multiplier if we want to cut down the numbers to be tested
        maxTrialFactor = (int) Math.ceil(maxTrialFactorMultiplier * Math.pow(n, ONE_THIRD));
        //		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
        if (maxTrialFactorMultiplier <= 1) {
            smallFactoriser.setMaxFactor(maxTrialFactor);
        }
        else {
            // factor out the 3, needed for the mod 3 tricks
            smallFactoriser.setMaxFactor(3);
        }
        // factor out all small factors if
        long nAfterTrial = smallFactoriser.findFactors(n, primeFactors);
        if (primeFactors == null && nAfterTrial != n)
            return nAfterTrial;
        n = nAfterTrial;

        // if number is already factorized return immediately without calling lehman
        if (n == 1)
            return n;

        // re-adjust the maximal factor we have to search for. If small factors were found by trial division,
        // which is quite often the case for arbitrary numbers, this cuts down the runtime dramatically.
        // TODO maybe we can do even better here?
        maxTrialFactor = (int) Math.ceil(maxTrialFactorMultiplier * Math.pow(n, ONE_THIRD));
        // effectively kMax is reduced by maxFactorMultiplier^2 ->
        final int kMax = (int) (Math.ceil(maxTrialFactor / maxFactorMultiplierCube));
        final long n4 = 4 * n;
        final double sqrtN = Math.sqrt(n);
        final int nMod4 = (int) (n % 4);
        // we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
        final long nPow2Third = ((long) maxTrialFactor) * ((long) maxTrialFactor);
        final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
        final double sqrt4N = 2 * sqrtN;
        int nMod3 =  (int)  (n % 3);
        // we want kn = 4kn = 1 mod 3, 2*2 = 1, 1*1 = 1
        int [] ksMod3 = {3,3-nMod3,nMod3};

        for (int i=15; i > 0; i -= 5) {
            int kBegin = i;
            for (int j = 0; j < 5; j++, kBegin -= 6)
            {
                kBegin = kBegin<0 ? kBegin + 15 : kBegin;
//                int kBegin = ksMod3[i];
                int xStepMultiplier = 1;
//                int xStepMultiplier = i == 1 ? 3 : 1;
                for (int k = kBegin; k <= kMax; k += 15) {
                    final double sqrt4kn = sqrt[k] * sqrt4N;
                    // adding a small constant to avoid rounding issues and rounding up is much slower then
                    // using the downcast and adding a constant close to 1. Took the constant from the yafu code
                    int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
                    // use only multiplications instead of division here
                    final float xRange = nPow1Sixth * sqrtInv[k];
                    final int xEnd = (int) (sqrt4kn + xRange);
                    // instead of using a step 1 here we use the mod argument from lehman
                    // to reduce the possible numbers to verify.
                    int xStep = 2;
                    if (k % 2 == 0) {
                        xBegin |= 1;
                    } else {
                        if (k * n % 4 == 3) {
                            xStep = 8;
                            xBegin = xBegin + PrimeMath.mod(7 - k * n - xBegin, 8);
                        } else {
                            xStep = 4;
                            xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
                        }
                    }
                    // for kn = 1 mod 3 we know x = 0 mod 3 -> try to adjust xBegin
                    while (xStepMultiplier == 3 && xBegin <= xEnd && xBegin % 3 != 0) {
                        xBegin += xStep;
                    }
                    // since the upper limit is the lower limit plus a number lower then 1
                    // in most of the cases checking the upper limit is very helpful
                    for (long x = xBegin; x <= xEnd; x += xStep * xStepMultiplier) {
                        // here we start using long values
                        // in java the trick to replace the multiplication with an addition does not help
                        final long x2 = x * x;
                        final long right = x2 - k * n4;
                        // instead of taking the square root (which is a very expensive operation)
                        // and squaring it, we do some mod arguments to filter out non squares
                        //				if (PrimeMath.isSquare(right)) {
                        if (isProbableSquare(right)) {
                            final long y = (long) Math.sqrt(right);
                            if (y * y == right) {
                                final long factor = PrimeMath.gcd(n, x - y);
                                if (factor != 1) {
                                    if (primeFactors == null || maxTrialFactorMultiplier > 1)
                                        return factor;
                                        // when maxFactorMultiplier == 1 we have done the trial division first -> the factor
                                        // has to be a prime factor, n/factor is of size < n^2/3 and can not be a composite factor
                                    else {
                                        primeFactors.add(factor);
                                        if (n != factor)
                                            primeFactors.add(n / factor);
                                        return 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // if we have not found a factor we still have to do the trial division phase
        if (maxTrialFactorMultiplier > 1)
            smallFactoriser.setMaxFactor(maxTrialFactor);
            n = smallFactoriser.findFactors(n, primeFactors);

        return n;
    }

    protected static boolean isProbableSquare(final long number) {
        if (isSquareMod_1024[(int) (number & 1023)]) {
            // using mod9_5_7_11 instead of hard 3456 coded number causes double run time
            // the % is expensive, but saves ~ 10% of the time, since the sqrt takes even more time
            // another mod filter gives not gain, in the YAFU impl it is used
            if (isSquareMod_9_5_7_11[(int) (number % 3456)]) {

                final long y = (long) Math.sqrt(number);
                return y * y == number;
            }
        }
        return false;
    }
    @Override
    public boolean returnsOnlyPrimeFactors(){
        return maxTrialFactorMultiplier <= 1;
    }

    @Override
    public String toString() {
        return "LehmanFactorFinder{" +
                "maxFactorMultiplier=" + maxTrialFactorMultiplier +
                '}';
    }
}
