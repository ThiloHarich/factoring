package factoring.trial.playgound;

import java.util.Collection;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFactCond {


    public long findPrimeFactors(long n, Collection<Long> primeFactors) {
//        for (int i = 2; i <= Math.sqrt(n); i++) {
        long initialN = n;
        long nDivI = n/3;
        for (int i = 3; i <= nDivI; i+=2) {
//            for (int i = 3; i*i <= n; i+=2) {
            nDivI = initialN/i;
//            if (n%i == 0) {
            if (nDivI * i == n) {
                primeFactors.add((long)i);
                n = n/i;
            }
        }
        return n;
    }

}
