package factoring;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialNewton implements Factorizer {


    @Override
    public int findFactor(long n) {
        double inversIApprox = 1./2;
        long nDivI = n/3 + 1;
        //since we have n/i we can compare with i; we do not have to call i*i <= n or i <= sqrt(n), which is more expensive
        for (int i = 3; i < nDivI; i+= 2) {
            // try to approximate 1/i from 1/(i-1) with the newton method to get around the expensive n/i calculation
            inversIApprox = inversIApprox * (2 - i * inversIApprox);
            // This provides better precision
            double error = 1 - i * inversIApprox;
            //					inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
            nDivI = Math.round(n*inversIApprox);
            double nDivIMuliplyI = nDivI*i;

            // if the error is to high we can fix it by doing an other iteration or may just in-/decrease nDivIMuliplyIAppox by i
            while (Math.abs(nDivIMuliplyI - n) >= i)
            {
                // do an other iteration
                inversIApprox = inversIApprox * (2 - i * inversIApprox);
                //						inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
                error = 1 - i * inversIApprox;
                nDivI = Math.round(n*inversIApprox);
                nDivIMuliplyI = nDivI*i;
            }
            if (nDivIMuliplyI == n)
                return i;
        }
        return -1;
    }

}
