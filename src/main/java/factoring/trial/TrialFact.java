package factoring.trial;

import factoring.FindPrimeFact;

import java.util.Collection;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFact extends FindPrimeFact {


    @Override
    public long findPrimeFactors(long n, Collection<Long> factors) {
//        for (int i = 2; i <= Math.sqrt(n); i++) {
        long initialN = n;
        for (int factor = 3; factor*factor <= initialN; factor = factor + 2) {
            while (n%factor == 0) {
                factors.add((long)factor);
                n /= factor;
            }
        }
        return n;
    }

}
