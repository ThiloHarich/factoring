package factoring.squfof;

import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import de.tilman_neumann.math.base.bigint.primes.probable.BPSWTest;
import de.tilman_neumann.math.factor.squfof.SquFoF31;
import de.tilman_neumann.math.factor.squfof.SquFoF63;
import factoring.Factorizer;
import factoring.math.PrimeMath;

import java.util.Collection;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class Squfof31 implements Factorizer {

    SquFoF31 squfo = new SquFoF31();
    public Squfof31(){
    }
    /**
     * In contrast to the trial division we can not be sure that a found factor is a prime
     * so we have to factorize it again
     * @param n
     * @return
     */
    public TreeMultiset<Long> findAllPrimeFactors(long n) {
        TreeMultiset<Long> primeFactors = TreeMultiset.create();
        TreeMultiset<Long> factors = findAllFactors(n);
        while (factors.size() > 0) {
            Multiset.Entry<Long> entry = factors.pollLastEntry();
            long factor = entry.getElement();
//            if (factor*factor > n || factor < 30)
//                primeFactors.add(factor, entry.getCount());
//            else
                {
                TreeMultiset<Long> newFactors = TreeMultiset.create();
                storeFactors(factor, newFactors);
                // if no factor found it is a prime factor
                if (newFactors.lastEntry().getElement() == factor) {
                    primeFactors.add(factor, entry.getCount());
                    newFactors.pollLastEntry();
                }
                else {
                    for (long newFactor : newFactors) {
                        factors.add(newFactor, entry.getCount());
                    }
                }
            }
        }
        return primeFactors;
    }

    /**
     *
     * @param n
     * @param factors
     * @return
     */
    public Collection<Long> storeFactors(long n, Collection<Long> factors) {
        // first make the number odd, this is required by the fermat factorizations
        while ((n & 1) == 0)
        {
            factors.add(2l);
            n = n >> 1;
        }
        if (n>1) {
            long remainder = findFactors(n, factors);
            if (remainder != 1)
                factors.add(remainder);
        }
        return factors;
    }

    public long findFactors(long n, Collection<Long> factors) {
        BPSWTest primeTest = new BPSWTest();
        if (primeTest.isProbablePrime(n))
            return n;
        else{
            long factor = squfo.findSingleFactor(n);
            factors.add (factor);
            return n/factor;
        }
    }
}
