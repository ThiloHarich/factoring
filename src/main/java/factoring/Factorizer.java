package factoring;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public interface Factorizer {

    default String printFactors(long n) {
        TreeMultiset<Long> factors = findAllPrimeFactors(n);

        List<String> s = new ArrayList<String>();
        for (Multiset.Entry<Long> entry : factors.entrySet()) {
            int exponent = entry.getCount();
            String part = "" + entry.getElement();
            part += exponent == 1 ? "" : "^" + exponent;
            s.add(part);
        }
        String join = StringUtils.join(s, " * ");
        return join;
    }

    default TreeMultiset<Long> findAllFactors(long n) {
        TreeMultiset<Long> factors = TreeMultiset.create();
        return (TreeMultiset<Long>) storeFactors(n, factors);
    }

    // TODO is this needed?
    default List<Long> findAllFactorsList(long n) {
        List<Long> factors = new ArrayList();
        return (List<Long>) storeFactors(n, factors);
    }
    default List<Long> findAllPrimeFactorsList(long n) {
        return  findAllFactorsList(n);
    }

    /**
     *
     * @param n
     * @return
     */
//    default List<Integer> findAllPrimeFactorsList(long n) {
//        List<Integer> primeFactors = new ArrayList<>();
//        List<Integer> factors = findAllFactorsList(n);
//        while (factors.size() > 0) {
//            int factor = factors.remove(factors.size() - 1);
//            List<Integer> newFactors = new ArrayList<>();
//            storeFactors(factor, newFactors);
//            if (newFactors.get(newFactors.size()-1) == factor)
//                primeFactors.add(factor);
//            else
//                factors.addAll(newFactors);
//        }
//        return primeFactors;
//    }

    default TreeMultiset<Long> findAllPrimeFactors(long n) {
        return findAllFactors(n);
    }
//    default TreeMultiset<Integer> findAllPrimeFactors(long n) {
//        TreeMultiset<Integer> primeFactors = TreeMultiset.create();
//        TreeMultiset<Integer> factors = findAllFactors(n);
//        while (factors.size() > 0) {
//            Multiset.Entry<Integer> entry = factors.pollLastEntry();
//            int factor = entry.getElement();
//            if (factor*factor > n || factor < 30)
//                primeFactors.add(factor, entry.getCount());
//            else {
//                TreeMultiset<Integer> newFactors = TreeMultiset.create();
//                storeFactors(factor, newFactors);
//                if (newFactors.pollLastEntry().getElement() == factor)
//                    primeFactors.add(factor, entry.getCount());
//                else {
//                    for (int newFactor : newFactors) {
//                        primeFactors.add(newFactor, entry.getCount());
//                    }
//                }
//            }
//        }
//        return primeFactors;
//    }


    /**
     * Returns a Collection of prime factors, starting with an empty set of factors.
     * @param n
     * @param factors
     * @return
     */
    Collection<Long> storeFactors(long n, Collection<Long> factors);


//    /**
//     * just for speed testing
//     * @param n
//     * @return
//     */
//    default int countAllFactors(long n) {
//        int count = 0;
//        // first make the number odd, this is required by the fermat factorizations
//        while ((n & 1) == 0)
//        {
//            count++;
//            n = n >> 1;
//        }
//        do {
//            findPrimeFactors(n, null);
//                // the last factor must be a prime
////                if(n > 1 && !BigInteger.valueOf(n).isProbablePrime(10))
////                    System.err.println("no factor found for " + n);
//                count++;
//                n=1;
//        }
//        while (n > 1);
//        return count;
//    }


}
