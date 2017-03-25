package factoring;

import com.google.common.primitives.Bytes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Thilo Harich on 06.03.2017.
 */
public class TrialUnroll implements Factorizer{

    /**
     * define the sequence of primes here
     */
    static int[] sievePrimes = {2,3,5};

    static int range = 1;
    static byte[] primeDist;
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

    public static void main (String[] args)
    {
        System.out.println("while (prime * prime < n) {");
        for (byte step : primeDist) {
            System.out.println("    prime += " + step + ";");
            System.out.println("    if (prime * prime > n)");
            System.out.println("        break;");
            System.out.println("    if (n % prime == 0)");
            System.out.println("        return prime;");
        }
        System.out.println("}");
    }
    @Override
    public int findFactor(long n) {
        for (int i=1; i < sievePrimes.length; i++) {
            if (n % sievePrimes[i] == 0)
                return sievePrimes[i];
        }
        int prime = 1;
        // find the right place to to start
        // 1 + primeDist[1] = 1 + 4 = 5 = next prime > sievePrimes 2,3
        // 1 + primeDist[1] = 1 + 6 = 7 = next prime > sievePrimes 2,3,5
        // 1 + primeDist[1] = 1 + 10 = 11 = next prime > sievePrimes 2,3,5,7
        for (int i = 1; prime * prime < n && i < primeDist.length; i++) {
            prime += primeDist[i];
            if (n % prime == 0)
                return prime;
        }

        /**
         * run main method and put in generated code here
         */
        while (prime * prime < n) {
            prime += 2;
            if (n % prime == 0)
                return prime;
            prime += 6;
            if (prime * prime > n)
                break;
            if (n % prime == 0)
                return prime;
            prime += 4;
            if (n % prime == 0)
                return prime;
            prime += 2;
            if (prime * prime > n)
                break;
            if (n % prime == 0)
                return prime;
            prime += 4;
            if (n % prime == 0)
                return prime;
            prime += 2;
            if (prime * prime > n)
                break;
            if (n % prime == 0)
                return prime;
            prime += 4;
            if (n % prime == 0)
                return prime;
            prime += 6;
            if (prime * prime > n)
                break;
            if (n % prime == 0)
                return prime;
        }
        return -1;
    }

}
