package factoring.trial;

import com.google.common.primitives.Bytes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Created by Thilo Harich on 06.03.2017.
 */
public class TrialUnroll extends TrialFactMod {

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
        System.out.println("while (factor <= maxFactor) {");
        for (byte step : primeDist) {
            System.out.println("    factor += " + step + ";");
            System.out.println("    if (factor > maxFactor)");
            System.out.println("        break;");
            System.out.println("    if (n % factor == 0)");
            System.out.println("    {");
            System.out.println("        factors.add((long)factor);");
            System.out.println("        n = n/factor;");
            System.out.println("    }");
        }
        System.out.println("}");
    }
    @Override
    public long findPrimeFactors(long n, Collection<Long> factors) {
        n = findSmallFactors(n, factors);

        // adjust limit, if not set
        int sqrtN = (int) Math.sqrt(n) + 1;
        int maxFactor = Math.min(sqrtN, limit);

        /**
         * run main method and put in generated code here
         */
        while (factor <= maxFactor) {
            factor += 2;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 6;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 4;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 2;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 4;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 2;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 4;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
            factor += 6;
            if (factor > maxFactor)
                break;
            if (n % factor == 0)
            {
                factors.add((long)factor);
                n = n/factor;
            }
        }
        return n;
    }

}
