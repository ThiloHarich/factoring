package factoring.fermat.lehman;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

import java.util.Collection;

/**
 * /**
 * Mod 15:
 * 15, 10, 5               -5 
 * 14, 9, 4  -> 9,4,14     -5
 * 13, 8, 3  -> 3,13,8     -5
 * 12, 7, 2  -> 12,7,2     -5
 * 11, 6, 1  -> 6,1,11     -5
 * for (int i=15; i> 0; i-=5)
 * 	for (int j=0; j<6; j++,i-=6)
 * 		for(l=0; l<3;l++,l i-=5)
 * 		if (i<0)
 * 			i+=15;
 * 
 * 15(5*3) , 9 (3*3), 3, 12, 6  -6 mod
 * 10, 5,                       -5
 * 14, 8, 2, 11                 -6 mod
 * 4, 13, 7, 1                  -6 mod
 * 
 * 15(5*3) , 9 (3*3), 3, 12, 6  -6 mod
 * 10, 4, 13, 7, 1,             -6 mod
 * 5, 14,  8, 2, 11,            -6 mod
 * 
 * for (int i=15; i> 0; i-=5)
 * 	for (int j=0; j<6; j++,i-=6)
 * 		if (i<0)
 * 			i+=15;

 * 
 * rest
 * In this algorithm we use an additional multiplier of 3 for k.
 * This means we first consider k' = 3*k and then the other values.
 * For k = 1 mod 3 we know that x must be dividable by 3.
 * 0^2 = 0 mod 3
 * 1^2 = 1 mod 3
 * 2^2 = 1 mod 3
 * x^2 - 4kn = y^2 mod 3
 * x^2 = y^2 + kn mod 3
 * assume n != 0 mod 3 by checking if 3 divides n before applying the lehman algorithm.
 * kn = 0 mod 3 -> every x is possible -> k = 0 mod 3
 * kn = 1 mod 3 -> y^2 = 0, x^2 = 1 -> x = 1,2 mod 3 and k = n Mod 3
 * kn = 2 mod 3 -> y^2 = 1, x^2 = 0 -> x = 0 mod 3 and k = -n Mod 3
 * If we analyze the equation mod 9 we see that for kn = 1 we have 3 (out of 9) solutions for x,
 * but for kn = 2 we only have 2 solutions.
 * So we first analyze kn = 0 then kn = 2 and then kn = 1
 * in the kn = 2 case we can multiply the step for x by 3
 *
 *
 * When k*n % 4 == 3 we can increase x by 8. This comes out if you investigate
 * x^2 - 4kn = y^2 mod 32.
 * This is a version of the lehman factorizationByFactors, which is a variant of the fermat
 * factorizationByFactors.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster then the java version of the yafu lehman factorizationByFactors.
 * If can be used for up to 41 bit. If is faster then other algorithms like SQUAFU for up to around 30 bits
 * for hard numbers with two or three (big) factors or random numbers. For finding one factor and for
 * finding all factors.
 *
 * The main performance features are:
 * By storing the square roots of the multiplier k the range for x (of x^2 - 4kn = y^2) can be calculated faster.
 * It also uses a version of trial division, where the multiple inverse of the primes are stored.
 * So instead of a division a multiplication is needed to find out if a number is dividable
 * by a prime.
 * Beside the arguments modulus 2 and 4 for calculating the step for x it investigates the
 * fermat equation for modulus 3 and 8. This can be seen as first considering k with a multiplier three.
 *
 * In the lehman algorithm in most of the cases i.e. k > n^1/3 / 16 the upper bound for x is less the
 * The lower bound plus 1. In this case at most one value of x has to be considered, but
 * the calculation of the lower and upper rage has to be done all the time.
 * Here we ignore the fact that we always increase x by 2 or 4.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 * <p>
 * The Hart variant uses just one x per multiplier k, this eliminates the determination of the
 * upper bound, but using it gives no extra speed. Why?
 * <p>
 * We need 2^42/3 = 2^14 = 16364 locations to store the squares and the inversive. This does not fit in the L1 cache,
 * but it is accessed in a sequential way -> still efficient
 * <p>
 * Open questions, possible improvements :
 * - can we get rid of storing the square roots? how can we calculate them efficiently?
 * - In most of the cases (more then 93%) for k only one or none the value x^2 -n needs to be calculated.
 * - can we extend the mod 3 observations to mod 9 or 15? Since the the fermat formula uses k*n it seems to be
 * too complicated to use it.
 * <p>
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanFactorFinderStep8Mod3 implements FactorizationOfLongs {

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

    float maxFactorMultiplier = 1.0f;
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
     * @param maxFactorMultiplierIn Defines the size of the factors the algorithm is first looking for.
     *                              The algorithm is working best for number where the lowest factor is
     *                              maxFactorMultiplier * n^1/3 or higher. If maxFactorMultiplier is 1 or lower we first do trial division
     *                              for numbers below n^1/3 then do the lehman factorizationByFactors above n^1/3. In this case the running time is
     *                              bounded by max(maximal factor, c*n^1/3 / log(n)).
     *                              For maxFactorMultiplier >= 1 we first inspect numbers above maxFactorMultiplier * n^1/3 with the lehman
     *                              factorizationByFactors. After completing this phase we do
     *                              trial division for numbers below maxFactorMultiplier * n^1/3<br>
     *                              Use<br>
     *                              maxFactorMultiplier <= 1<br>
     *                              - if you do not know anything about the numbers to factorize<br>
     *                              - if you know the maximal factors can not exceed n^maxFactorMultiplier < n^1/3<br>
     *                              here all found factors are prime factors and will be stored in the primeFactors Collection.
     *                              maxFactorMultiplier = 3 (this is the optimal value) <br>
     *                              - if you know for most of the numbers the maximal factors will exceed 3*n^1/3<br>
     *                              In the last case {@link #findFactors(long, Collection)} might return a composite number
     */
    public LehmanFactorFinderStep8Mod3(int bits, float maxFactorMultiplierIn) {
        if (bits > 41)
            throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
        maxFactorMultiplier = maxFactorMultiplierIn < 1 ? 1 : maxFactorMultiplierIn;
        maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(1L << bits, ONE_THIRD));
        maxFactorMultiplierCube = maxFactorMultiplier * maxFactorMultiplier * maxFactorMultiplier;
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
        System.out.println("sqrt tables for max trial factor " + maxFactorMultiplier + " built bytes used : " + sqrt.length);
    }

    @Override
    public long findFactors(long n, Collection<Long> primeFactors) {
        if (n > 1l << 41)
            throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
        // with this implementation the lehman part is not slower then the trial division
        // we only apply the multiplier if we want to cut down the numbers to be tested
        maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
        //		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
        if (maxFactorMultiplier <= 1) {
            smallFactoriser.setMaxFactor(maxTrialFactor);
        }
        else {
            // fcator out the 3
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
        maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
        // effectively kMax is reduced by maxFactorMultiplier^2 ->
        final int kMax = (int) (Math.ceil(maxTrialFactor / maxFactorMultiplierCube));
        final long n4 = 4 * n;
        byte n4Mod9 = (byte) (n4 % 9);
        byte nMod9 = (byte) (n % 9);
        final double sqrtN = Math.sqrt(n);
        final int nMod4 = (int) (n % 4);
        // we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
        final long nPow2Third = ((long) maxTrialFactor) * ((long) maxTrialFactor);
        final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
        // surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
        // for k=2 we know that x has to be even and results in a factor more often
        final double sqrt4N = 2 * sqrtN;
//        int kMod9 = 1;
        byte kn4Mod9 = 0;
        int nMod3 =  (int)  (n % 3);
        // we want kn = 4kn = 1 mod 3, 2*2 = 1, 1*1 = 1
//        int inv5 = invMod9[5];
//        int inv8 = invMod9[8];
//        int nInvMod9 = invMod9[nMod9];
//        invMod9[0] = 9;
//        int nInvI = nInvMod9;
//        for (int i=1; i<9; i++)
//        {
//            nInvI += nInvMod9;
//            nInvI = nInvI > 9 ? nInvI - 9 : nInvI;
//        }
//        int [] ksMod9 = {9,(6*nInvMod9)%9,(3*nInvMod9)%9,(1*nInvMod9)%9,(4*nInvMod9)%9,(7*nInvMod9)%9,(2*nInvMod9)%9,(5*nInvMod9)%9,(8*nInvMod9)%9};
//        int [] ksMod9 = {9,6,3,8,7,5,4,2,1};
        // first use k = 0 mod 3 (each x possible), then k=1 mod 3 ( 2 of 3 x possible), k=2 mod 3 -> x=0 mod 3
        int [] ksMod3 = {3,3-nMod3,nMod3};
//        for (int j=1; j <= 2; j++)
            for (int i=0; i <3; i++)
            {
                int kBegin = ksMod3[i];
//                int kBegin = ksMod9[i];
                // for kn = 1 mod 3 x = 0 mod 3 -> we can multiply the x step by 3
                int xStepMultiplier = i == 1 ? 3 : 1;
//        for (int k = 1; k <= kMax; k++) {
            float xRange = nPow1Sixth * sqrtInv[kBegin];
            for (int k = kBegin; k <= kMax; k+=3) {
                final double sqrt4kn = sqrt[k] * sqrt4N;
                // adding a small constant to avoid rounding issues and rounding up is much slower then
                // using the downcast and adding a constant close to 1. Took the constant from the yafu code
                int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
                // use only multiplications instead of division here
                xRange = xRange <= 3 ? 3 : nPow1Sixth * sqrtInv[k];
                final int xEnd = (int) (sqrt4kn + xRange);
                // instead of using a step 1 here we use the mod argument from lehman
                // to reduce the possible numbers to verify.
                int xStep = 2;
                if (k % 2 == 0 ) {
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
                while (xStepMultiplier == 3 && xBegin % 3 != 0)
                    xBegin += xStep;
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
                                if (primeFactors == null || maxFactorMultiplier > 1)
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
        // if we have not found a factor we still have to do the trial division phase
        if (maxFactorMultiplier > 1)
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
        return maxFactorMultiplier <= 1;
    }

    @Override
    public String toString() {
        return "LehmanFactorFinder{" +
                "maxFactorMultiplier=" + maxFactorMultiplier +
                '}';
    }
}
