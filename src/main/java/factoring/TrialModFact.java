package factoring;

import com.google.common.primitives.Bytes;
import org.apache.commons.lang3.ArrayUtils;

import java.lang.reflect.Array;
import java.util.ArrayList;

import java.util.Arrays;
import java.util.List;


/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialModFact implements Factorizer {

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

    @Override
    public int findFactor(long n) {
        for (int i=1; i < sievePrimes.length; i++) {
            if (n % sievePrimes[i] == 0)
                return sievePrimes[i];
        }
        int prime = 1;
        // find the right place to to start
        // 1 + primeDist[1] = 1 + 6 = 7 = next prime > sievePrimes
        // 1 + primeDist[1] = 1 + 10 = 11 = next prime > sievePrimes
        for (int i = 1; prime * prime < n && i < primeDist.length; i++) {
            prime += primeDist[i];
            if (n % prime == 0)
                return prime;
        }

        while (prime * prime <= n) {
            for (int i = 0; i < primeDist.length; i++) {
                prime += primeDist[i];
                if (prime * prime > n)
                    break;
                if (n % prime == 0)
                    return prime;
            }
        }
        return -1;
    }
}
