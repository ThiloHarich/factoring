package factoring;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public interface Factorizer {

    default String printFactors(long n) {
        TreeMultiset<Integer> factors = findAllFactors(n);

        List<String> s = new ArrayList<String>();
        for (Multiset.Entry<Integer> entry : factors.entrySet()) {
            int exponent = entry.getCount();
            String part = "" + entry.getElement();
            part += exponent == 1 ? "" : "^" + exponent;
            s.add(part);
        }
        String join = StringUtils.join(s, " * ");
        return join;
    }

    default TreeMultiset<Integer> findAllFactors(long n) {
        TreeMultiset<Integer> factors = TreeMultiset.create();
        return (TreeMultiset<Integer>) storeFactors(n, factors);
    }

    default List<Integer> findAllFactorsList(long n) {
        List<Integer> factors = new ArrayList<>();
        return (List<Integer>) storeFactors(n, factors);
    }

    /**
     * This should not be part of the interface
     * @param n
     * @param factors
     * @return
     */
    default Collection<Integer> storeFactors(long n, Collection<Integer> factors) {
        // first make the number odd, this is required by the fermat factorizations
        while ((n & 1) == 0)
        {
            factors.add(2);
            n = n >> 1;
        }
        do {
            if (n == 1)
                return factors;
//            int factor = findFactor(n);
            int factor = findFactor(n);
            if (factor != -1) {
                factors.add(factor);
                n = n / factor;
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
        while (n > 1);
        return factors;
    }


    /**
     * just for speed testing
     * @param n
     * @return
     */
    default int countAllFactors(long n) {
        int count = 0;
        // first make the number odd, this is required by the fermat factorizations
        while ((n & 1) == 0)
        {
            count++;
            n = n >> 1;
        }
        do {
            long factor = findFactor(n);
            if (factor != -1) {
                count++;
                n = n / factor;
            }
            else
            {
                // the last factor must be a prime
//                if(n > 1 && !BigInteger.valueOf(n).isProbablePrime(10))
//                    System.err.println("no factor found for " + n);
                count++;
                n=1;
            }
        }
        while (n > 1);
        return count;
    }


    int findFactor(long n);
}
