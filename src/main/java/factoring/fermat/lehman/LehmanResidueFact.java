package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.fermat.residue.FermatResiduesArray;
import factoring.math.PrimeMath;
import factoring.trial.TrialWithPrimesFact;

import java.util.Collection;

/**
 * This is a variant of the Lehman factorizationByFactors were we only take into account solutions for x modulo
 * certain powers of 2. The original lehman algorithm does this only for modulo 2 and 4. This
 * variant can be applied also to higher mods. The possible solutions for such mods are stored in an
 * array. Lehman has solved the equations for 2 and 4 explicitly.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanResidueFact extends FermatFact {

    double balanceTrial = 1;
    int mod =9;
    FermatResiduesArray residues = new FermatResiduesArray(mod);


    @Override
    public long findFactors(long n, Collection<Long> factors) {
        TrialWithPrimesFact smallFactoriser = new TrialWithPrimesFact();
        int maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        smallFactoriser.setMaxFactor(maxTrialFactor);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<maxTrialFactor)
            return n;

        maxTrialFactor =  (int) Math.ceil(balanceTrial * Math.pow(n, 1.0/3));
        if (PrimeMath.isSquare(n)){
            long x = PrimeMath.sqrt(n);
            if (x*x == n) {
                factors.add(x);
                return x;
            }
        }
        int multiplier = 4;
        int multiplierSqrt = (int) Math.sqrt(multiplier);
        double balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
        // effectively the outer loop (height) is cut down by a balanceTrial^2
        // and the inner loop (width) is multiplied by balanceTrial^2
        // by doing this, it is much more important to reduce the candidates of each level.
        int kMax = (int) (maxTrialFactor / balanceTrialCube);

        double nPow2Third = maxTrialFactor * maxTrialFactor;
        for (int k = 1; k <= kMax; k++) {
            long kn = k * n;
            long k4n = kn * multiplier;
            residues.init(k4n, new int[]{mod});
            int k4nMod = (int) (k4n % mod);
            // we move the multiplier out of the calculation of the square root.
            // so we can reuse the result for calculating the end of the loop
            // unfortunately this means we are restricted to squares as multipliers
            // trick found https://github.com/DarkenCode/yafu/blob/master/factor/LehmanClean.c
            double sqrtKN = Math.sqrt(k4n);
            long xBegin = (long) (Math.ceil(sqrtKN));
            // if the xSol are all solutions for a mod x in ascending order,7
            // find the first index which has a value within the range we want to sieve. This is x >= xBegin % mod
            int xBeginIndex = 0;
            int xBeginMod = (int) (xBegin % mod);
            for (; residues.x4Mod[k4nMod][xBeginIndex] >= 0 && residues.x4Mod[k4nMod][xBeginIndex] < xBeginMod; xBeginIndex++) ;

            double xRange = (nPow2Third / multiplierSqrt) / sqrtKN;
            double xEnd = sqrtKN + xRange;
//            long xBeginFloor = xBegin & (Integer.MAX_VALUE - mod + 1);
            long xBeginFloor = xBegin - xBeginMod;
            for (int xIndex = xBeginIndex; residues.x4Mod[k4nMod][xIndex] >= 0; xIndex++) {
                long x = xBeginFloor + residues.x4Mod[k4nMod][xIndex];
                while (x <= xEnd) {
                    long x2 = x * x;
                    long right = x2 - k4n;
                    if (PrimeMath.isSquare(right)) {
                        long y = (long) Math.sqrt(right);
                        long factor = PrimeMath.gcd(n, x - y);
                        if (factor != 1) {
                            factors.add(factor);
                            return n / factor;
                        }
                    }
                    x += mod;
                }
            }
            // now we do the same loop as above but ensure that the values are above xBegin
            xBeginFloor += mod;
            for (int xIndex = 0; xIndex < xBeginIndex; xIndex++) {
                long x = xBeginFloor + residues.x4Mod[k4nMod][xIndex];
                while (x <= xEnd) {
                    long x2 = x * x;
                    long right = x2 - k4n;
                    if (PrimeMath.isSquare(right)) {
                        long y = (long) Math.sqrt(right);
                        long factor = PrimeMath.gcd(n, x - y);
                        if (factor != 1) {
                            factors.add(factor);
                            return n / factor;
                        }
                    }
                    x += mod;
                }
            }
        }
        return n;
    }
}
