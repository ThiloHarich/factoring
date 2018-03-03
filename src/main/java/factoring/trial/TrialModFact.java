package factoring.trial;

import com.google.common.primitives.Bytes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;


/**
 * This implementation uses a small set of primesInv, calculates all
 * possible primesInv in the range of the product of all the primesInv by sieving.
 * Then this list of possible primesInv (which also contains non primesInv) is used
 * to see if the number n is dividable by.
 * TODO do a proper Wheel factorization implementation.
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialModFact extends TrialFact {

    static int[] primes;
    static int[] sievePrimes = {2,3,5,7};
    static int range = 1;

    static {
        // set up the
        Arrays.stream(sievePrimes).forEach(p -> range *= p);
        boolean[] nonPrimes=  new boolean [range];


        nonPrimes[0] = true;
        for (int i=0; i < sievePrimes.length; i++) {
            if (!nonPrimes[sievePrimes[i]]) {
                for (int j = sievePrimes[i]; j < range; j += sievePrimes[i]) {
                    nonPrimes[j] = true;
                }
            }
        }
        List<Integer> primesList = new ArrayList();
        for (int i = 0; i < range; i++) {
            if (!nonPrimes[i]) {
                primesList.add(i);
            }
        }
        // convert to int
        primes = primesList.stream().mapToInt(i->i).toArray();
    }
    int factor;

    public void setLimit(int limit) {
        this.limit = limit;
    }

    int limit = Integer.MAX_VALUE;

    @Override
    public long findPrimeFactors(long n, Collection<Long> factors) {
        n = findSmallFactors(n, factors);

        // adjust limit, if not set
        int sqrtN = (int) Math.sqrt(n) + 1;
        int maxFactor = Math.min(sqrtN, limit);

        for (int i = 0; i < primes.length; i++) {
            factor = primes[i];
            while (factor <= maxFactor) {
                if (n % factor == 0 && factor != 1) {
                    factors.add((long) factor);
//                    n = n/factor;
                }
                factor += range;
            }
        }
        // TODO find if the bigger factors can be split up.
        return n;
    }

    public long findSmallFactors(long n, Collection<Long> factors) {
        // first check if the sieve primesInv divide n
        for (int i = 1; i < sievePrimes.length; i++) {
            while (n % sievePrimes[i] == 0) {
                factors.add((long) sievePrimes[i]);
                n = n / sievePrimes[i];
            }
        }
//        for (int i = 0; factor * factor < n && i < primesInv.length; i++) {
//            factor = primesInv[i];
//            while (n % factor == 0) {
//                factors.add((long) factor);
//                n = n / factor;
//            }
//        }
        return n;
    }
}
