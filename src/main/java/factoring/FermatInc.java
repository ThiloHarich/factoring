package factoring;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class FermatInc implements Factorizer {

    public Collection<Integer> storeFactors(long n, Collection<Integer> factors) {
        // first make the number odd, this is required by the fermat factorizations
        while ((n & 1) == 0)
        {
            factors.add(2);
            n = n >> 1;
        }
        List<Long> candidates = new ArrayList<>();
        candidates.add(n);
        do {
            Long f = candidates.remove(candidates.size()-1);
            long factor = findFactor(f);
            if (factor != -1) {
                addProbableFactor(factors, candidates, factor);
                addProbableFactor(factors, candidates, (f / factor));
            }
            else
            {
                // the last factor must be a prime
//                if(n > 1 && !BigInteger.valueOf(n).isProbablePrime(10))
//                    System.err.println("no factor found for " + n);
                factors.add((int)n);
                n=1;
            }
        }
        while (candidates.size() > 0);
        return factors;
    }

    public void addProbableFactor(Collection<Integer> factors, List cadidates, long factor) {
        if (factor <= 1)
            return;
        if (BigInteger.valueOf(factor).isProbablePrime(10))
        {
            factors.add((int)factor);
        }
        else {
            cadidates.add(factor);
        }
    }

    @Override
    public int findFactor(long n) {
        long sqrtN = (long) Math.ceil(Math.sqrt(n));
//        (3+x/3)/2 = 3/2 + n/6
        for (long x = sqrtN; x <= 2 + n/6; x++) {
            long right = x*x - n;
//            long sqrtR = MyMath.sqrt(right);
//            long error = n + sqrtR*sqrtR;
//            right = n + (x-error)*(x-error);

            if(MyMath.isSquare(right))
            {
                long y = (long) Math.sqrt(right);
                return (int) (x+y);
            }

        }
        return -1;
    }
}
