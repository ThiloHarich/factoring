package factoring;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFact implements Factorizer {

//    public Collection<Integer> storeFactors(long n, Collection<Integer> factors) {
//        // first make the number odd, this is required by the fermat factorizations
//        while ((n & 1) == 0)
//        {
//            factors.add(2);
//            n = n >> 1;
//        }
//        int factor = 2;
//        do {
//            factor = findFactor(n, factor);
//            if (factor != -1) {
//                factors.add(factor);
//                n = n / factor;
//            }
//            else
//            {
//                // the last factor must be a prime
////                if(n > 1 && !BigInteger.valueOf(n).isProbablePrime(10))
////                    System.err.println("no factor found for " + n);
//                factors.add((int)n);
//                n=1;
//            }
//        }
//        while (n > 1);
//        return factors;
//    }

    @Override
    public int findFactor(long n) {
//        for (int i = 2; i <= Math.sqrt(n); i++) {
        for (int i = 3; i*i <= n; i+=2) {
//            for (int i = 3; i*i <= n; i+=2) {
            if (n%i == 0)
                return i;
        }
        return -1;
    }

    public int findFactor(long n, int tested) {
        int lower = tested == 2 ? 3 : tested + 2;
        for (int i = lower; i*i <= n; i+=2) {
            if (n%i == 0)
                return i;
        }
        return -1;
    }
}
