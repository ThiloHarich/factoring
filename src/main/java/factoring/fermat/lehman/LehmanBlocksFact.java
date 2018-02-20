package factoring.fermat.lehman;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialFactMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanBlocksFact extends FermatFact {

    TrialFactMod smallFactoriser = new TrialFactMod();
    int factor = 1;

    @Override
    public long findFactors(long n, Collection<Long> factors) {
        int limit =  (int) Math.ceil(Math.pow(n, 1.0/3));
        smallFactoriser.setLimit(limit);
        n = smallFactoriser.findPrimeFactors(n, factors);

        if (n<= limit) {
            n = smallFactoriser.findPrimeFactors(n, factors);
        }
        else
        {
            if (PrimeMath.isSquare(n)){
                long x = PrimeMath.sqrt(n);
                if (x*x == n) {
                    factors.add(x);
                    return x;
                }
            }
            long multiplier = 4;
            double n6 = Math.pow(n, 1.0/6) / 4;
            int i =1;
            for (int block = 1; i <= limit; block++) {
                // we calculate the range of the loop not per level
                // since we have to calculate the SQRT this will take 3 times more time than calculating
                // checking the number. We calculate the end of the loop per block.
                // a block is defined such that we have nearly the same number of operations, which is
                // n^1/6/4 < 16 for number n < 2^48. For the lehman algorithm this is the maximal number
                // for n since 2^48 * 2^(48/3) = 2^64.
                long xRange = (long) (n6 / block) + 1;
                for(i=block*block; i<(block+1)*(block+1); i++) {
                    long in = multiplier * i * n;
                    long xBegin = PrimeMath.sqrt(in);
                    for (long x = xBegin; x <= xBegin + xRange; x++) {
                        long x2 = x * x;
                        long right = x2 - in;
                        if (PrimeMath.isSquare(right)) {
                            long y = (long) Math.sqrt(right);
                            long factor = PrimeMath.gcd(n, x - y);
                            if (factor != 1) {
                                factors.add(factor);
                                return n / factor;
                            }
                        }
                    }
                }
            }
        }
        return n;
    }
}
