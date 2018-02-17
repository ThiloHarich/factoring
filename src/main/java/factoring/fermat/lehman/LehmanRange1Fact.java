package factoring.fermat.lehman;

import factoring.FindPrimeFact;
import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanRange1Fact extends FindPrimeFact {

    final static float balanceTrial = 1.4f;
    final static double ONE_THIRD = 1.0/3;

    static float balanceTrialCube = balanceTrial * balanceTrial;
    static int kMax = (int) (Math.ceil(Math.pow(Long.MAX_VALUE, ONE_THIRD)) / balanceTrialCube);

    static float [] sqrt = new float[kMax + 1];
    static float [] sqrtInv = new float[kMax + 1];
    static {
        for(int i=1; i<sqrt.length; i++)
        {
            float sqrtI = (float) Math.sqrt(i);
            sqrt[i] = sqrtI;
            sqrtInv[i] = 1.0f / sqrtI;
        }
    }
    final static int xBeginMod4 = Integer.MAX_VALUE - 3;
    
    // we store all the variables we need as fields to make the methods easy
    private Collection<Long> factors;
    private long n;
    int multiplier;
    private int multiplierSqrt;
    private long n4;
    private int nMod4;
    private double nPow1Sixth;
    private double sqrtN;
    int k = 1;

    @Override
    public long findPrimeFactors(long nIn, Collection<Long> factorsIn) {
        this.factors = factorsIn;
        this.n = nIn;
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        double maxTrialFactor =  Math.ceil(balanceTrial * Math.pow(n, ONE_THIRD));
        smallFactoriser.setMaxFactor((int) maxTrialFactor);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<maxTrialFactor)
            return n;

        if (PrimeMath.isSquare(n)){
            long x = PrimeMath.sqrt(n);
            if (x*x == n) {
                factors.add(x);
                return x;
            }
        }
        // readjust the maximal factor
        // TODO we can approximate the max factor here
        maxTrialFactor =  balanceTrial * Math.pow(n, ONE_THIRD);
        kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
        int k4Range1 = (int) (balanceTrial * balanceTrial * balanceTrial * kMax / 16);
        multiplier = 4;
        n4 = n * multiplier;
        multiplierSqrt = 2;
        sqrtN = Math.sqrt(n);
        nMod4 = (int) (n % 4);
        double nPow2Third = maxTrialFactor * maxTrialFactor;
        // TODO is division by 4 more efficient if we use int?
        // TODO use float avoid division
        nPow1Sixth = (nPow2Third / 4) / sqrtN;
        int nMod3 = (int) (n % 3);

//        x^2 - 4kn = y^2 mod 3
//        x^2 - kn = y^2 mod 3
//        x^2 = y^2 mod 3 , k*n == 0 -> k = 0 mod 3 since n != 0 mod 3 -> all solutions
        long factor = getFactor();
        if (factor > 0) return factor;
        factor = getFactorOneX();
        if (factor > 0) return factor;
//        long factor = getFactor(3, kMax, 3);
//        if (factor > 0) return factor;
//        x^2 - 1 = y^2 mod 3 , k*n == 1 -> x=1,2 mod 3 -> k = n^-1 mod 3 -> k = n mod 3
//        mod 9 there are only 2 possible solutions as well
//        factor = getFactor(nMod3, kMax, 3, 1, 3);
//        if (factor > 0) return factor;
//        factor = getFactor(nMod3, kMax, 3, 2, 3);
//        if (factor > 0) return factor;
//        x^2 - 2 = y^2 mod 3 , k*n == 2 -> x=0 mod 3   -> k = 2* n^-1 mod 3  -> k = 2n mod 3
//        factor = getFactor(3 - nMod3, kMax, 3, 0, 3);
//        if (factor > 0) return factor;

        return n;
    }

    public long getFactor() {
        double xRange;
        do {
            long kn4 = k * n4;
            double sqrtKn = sqrt[k] * sqrtN;
            double sqrt4kn = multiplierSqrt * sqrtKn;
            int x = (int) (Math.ceil(sqrt4kn - 0.001));
            // use only multiplications instead of division here
            xRange = nPow1Sixth * sqrtInv[k];
            long xEnd = (long) (sqrt4kn + xRange);
            // TODO use
            x = nextX(1, nMod4, k, x);
            while( x <= xEnd) {
                long x2 = x * x;
                long right = x2 - kn4;
                if (PrimeMath.isSquare(right)) {
                    long y = (long) Math.sqrt(right);
                    long factor = PrimeMath.gcd(n, x - y);
                    if (factor != 1) {
                        factors.add(factor);
                        if (n / factor != 1)
                            factors.add(n / factor);
                        return 1;
                    }
                }
                if (k % 2 == 0) {
                    x += 2;
                }
                else{
                    x += 4;
                }
            }
            k++;
        }
        while(xRange > 1);
        return -1;
    }
    public long getFactorOneX() {
        for (; k <= kMax; k++) {
            double sqrtKn = sqrt[k] * sqrtN;
            double sqrt4kn = multiplierSqrt * sqrtKn;
            int x = (int) (Math.ceil(sqrt4kn - 0.001));
            // check x mod 2 and 4
            int kMod2 = k % 2;
            if ((x%2==0) ^ (kMod2 ==0)) {
//                if (kMod2==0 ||  (k + nMod4) % 4 == (x % 4) )
                {
                    long x2 = x * x;
                    long kn4 = k * n4;
                    long right = x2 - kn4;
                    if (PrimeMath.isSquare(right)) {
                        long y = (long) Math.sqrt(right);
                        long factor = PrimeMath.gcd(n, x - y);
                        if (factor != 1) {
                            factors.add(factor);
                            if (n / factor != 1)
                                factors.add(n / factor);
                            return 1;
                        }
                    }
                }
            }
        }
        return -1;
    }

    public int nextX(int xStep, int nMod4, int k, int x) {
        if (k % 2 == 0) {
            if (x % 2 == 0)
                x += xStep;
        }
        else{
            // Instead of doing the while loop we may use the inverse
            int knMod = (k + nMod4) % 4;
            while ((x % 4) != knMod)
                x += xStep;
        }
        return x;
    }
}
