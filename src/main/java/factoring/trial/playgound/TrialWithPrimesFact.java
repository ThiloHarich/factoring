package factoring.trial.playgound;

import factoring.FactorizationOfLongs;

import java.util.Collection;

/**
 * This implementation is generating a list of all primesInv up to a limit and will then check if the
 * number is dividable. Here we us a limit of 65536=2^16.
 * We can only factorize numbers up to 2^32.
 * When calling it with bigger numbers only prime factors below
 * 2^16 were added to the factors. {@link #findFactors(long, Collection)} then might return a
 * composite number. This strictly means it can not be applied for factorizationByFactors big numbers,
 * but it can be applied for finding the prime factors below 2^16 and applying other algorithms
 * to factor the returned remainder. This is exactly what the lehman and hart algorithms need.
 * We have choosen 2^16 because when factorizationByFactors long numbers by the lehman method, they have to be
 * lower n =  2^48 = (2^16)^3, which means we can find all primesInv below n^1/3 with this algorithm.
 *
 * Create a new instance, unless you want to have control over the search interval by using {@link #setMaxFactor(int)}.
 *
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialWithPrimesFact implements FactorizationOfLongs {

    static int[] prime = new int [6543]; //the 6542 primesInv up to 65536=2^16, then sentinel 65535 at end

    static
    {
        // TODO we do not need 2
        int i,j,k;
        prime[0]=2;
        prime[1]=3;
        prime[2]=5;
        k=3;
        for(i=7; i<65536; i+=2){
            boolean isPime = true;
            for(j=0; prime[j]* prime[j] <= i && isPime; j++){
                if(i%prime[j]==0)
                    isPime = false;
            }
            if (isPime) {
                prime[k] = i;
                k++;
            }
        }
        assert(k==6542);
        prime[k] = 65535; //sentinel
        System.out.printf("Prime table[0..%d] built: ", k);
        for(i=0; i<20; i++){ System.out.printf("%d,", prime[i]); }
        System.out.printf("%d,...,%d,(%d)\n", prime[20],prime[6541],prime[6542]);
    }

    private int maxFactor = 65535;
    // we keep the index so that we can continue with the factorizationByFactors after applying
    // other algorithms like lehman, not so nice
    int primeIndex = 0;

    public void setMaxFactor (int maxFactor) {
        if (maxFactor > 65535)
            throw new IllegalArgumentException("the maximal factor has to be lower then 65536");
        this.maxFactor = maxFactor;
    }

    @Override
    public long findFactors(long n, Collection<Long> primeFactors) {
        for (; prime[primeIndex] < maxFactor; primeIndex++) {
//            if (maxFactorIndex == 6542)
//                System.out.println();
            while (n%prime[primeIndex] == 0) {
                primeFactors.add((long)prime[primeIndex]);
                n /= prime[primeIndex];
            }
        }
        return n;
    }
}
