package factoring;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFactCond implements Factorizer {


    @Override
    public int findFactor(long n) {
//        for (int i = 2; i <= Math.sqrt(n); i++) {
        long nDivI = n/3;
        for (int i = 3; i <= nDivI; i+=2) {
//            for (int i = 3; i*i <= n; i+=2) {
            nDivI = n/i;
            if (n%i == 0)
                return i;
        }
        return -1;
    }

}
