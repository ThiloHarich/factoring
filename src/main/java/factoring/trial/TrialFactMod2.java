package factoring.trial;

import com.google.common.primitives.Bytes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;


/**
 * Is a copy of Mod fermat factorizer.
 * But uses the possible mods and always add the mod to one possible mod.
 * Since this step is constant it might be iterated faster, but is not.
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFactMod2 extends TrialFact {

    static byte[] primeDist;
    static int[] sievePrimes = {2,3,5};
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
        List<Byte> primesDist = new ArrayList();
        int lastPrime = -1;
        for (int i = 0; i < range; i++) {
            if (!nonPrimes[i]) {
                primesDist.add((byte)(i - lastPrime));
                lastPrime = i;
            }
        }
        // convert to int
//        primeDist = primesDist.stream().mapToInt(i->i).toArray();
        primeDist = Bytes.toArray(primesDist);
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

            for (int i = 0; i < primeDist.length; i++) {
                factor += primeDist[i];
                int factor2 = factor;
                while (factor2 <= maxFactor) {
                    while (n % factor2 == 0) {
                        factors.add((long) factor2);
                        n = n / factor2;
                    }
                    factor2 += range;
                }
            }
        return n;
    }

    public long findSmallFactors(long n, Collection<Long> factors) {
        // first check if the sieve primes divide n
        for (int i = 1; i < sievePrimes.length; i++) {
            while (n % sievePrimes[i] == 0) {
                factors.add((long) sievePrimes[i]);
                n = n / sievePrimes[i];
            }
        }
        factor = 1;
//        int factor = maxFactorTested;
        // find the right place to to start with the sieve
        // 1 + primeDist[1] = 1 + 6 = 7 = next factor > sievePrimes
        // 1 + primeDist[1] = 1 + 10 = 11 = next factor > sievePrimes
        for (int i = 1; factor * factor < n && i < primeDist.length; i++) {
            factor += primeDist[i];
            while (n % factor == 0) {
                factors.add((long) factor);
                n = n / factor;
            }
        }
        return n;
    }
}
