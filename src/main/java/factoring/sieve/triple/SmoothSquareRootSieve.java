package factoring.sieve.triple;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv31Barrett;
import de.tilman_neumann.util.SortedMultiset;

import java.math.BigInteger;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertTrue;

/**
 * Starting with the square t^2 , t = ceil(sqrt(n/s)) we can build (t-c)*(t+c) = t^2 - c^2.
 * With c^2 we can reduce the size of s*(t-c)*(t+c) - n to O(n^1/4). But t-c and t+c are no square and do not have to be smooth.
 * We can have a lookup table with parameters t and c (and the size of the factor base) giving back smooth numbers in O(log (n)).
 * If the size of s is within a small? range, the size of t is basically bound to n and such to a fixed factor base.
 * Since the prob for one or two values being smooth is sub polynomial , the size of (t-c)*(t+c) - n is O(n^1/4 + e)
 * with a (asymptotically) small e.
 *
 * If t itself has a smooth part c we know that t-c and t+c also have the factor c.
 * But if we know that c is of size O(n^e) and c divides t, we know that t + ic has the factor c.
 * We only need to factor (t+i*c)/c = O(n^(1/2-e))
 * r = s*(t-c)*(t+c) - n
 * r = s*(t^2 - c^2) - n
 * 0 = s*t^2 - s*c^2 - n mod p
 * s*c^2 = s*t^2 - n  mod p
 * c^2 = t^2 - n*s^-1 mod p
 * 1) we look for a smooth part's u (of size sqrt(t)) for t = ceil (sqrt( n/s)) u ~ n^e , e.g. n^1/8
 * 2) we calculate c = sqrt((s * t * t - n)/s) , i' = c/u
 * 3) find define a function  f_l(i) = t - u*(i' - i) and f_h = t + u*(i' + i)    where u*i ~ c , f_h ~ n^(1/2 -e)
 * 4) build r(i) = s*(f_l(i))*(f_h(i)) - n  r ~ n^(1/4 + e)*i
 * 5) we sieve over f_l, f_h and r  , f_l * f_h * r = n^(1 - 2e + 1/4 + e) = n^(5/4 -e)*i
 * for e= 1/8 f_l,f_h have the same size n^3/8, r has size n^3/8 *i
 * 6) factorize the product   (t'-d)*(t'+d)*r over the factor base
 * 6) if it is smooth factorize t-c', t+c', r
 * 7) combine the relations to a solution x^2 - n = y^2 mod n
 *
 * we do 3 by sieving with the primes of the smooth part of t
 *
 * 1) since t does not has to be completely smooth like in the QS this should be fast enough
 * 2) we have to find two numbers of size sqrt(t) < n^1/4
 * 3) r is of size O(n^1/4)
 * This is instead of looking for one number of size n^1/2 we have to find 3 numbers of size n^1/4.
 *
 */
public class SmoothSquareRootSieve extends FactorAlgorithm {

    	static int [] primes = {2,     3,     5,     7,    11,    13,    17,    19,    23,    29
    	,    31,    37,    41,    43,
	   47,    53,    59,    61,    67
    	,    71,    73,    79,    83,    89,    97,   101,   103,   107,
	  109,   113,   127,   131,   137,   139,   149,   151,   157,   163,   167,   173,   179,   181};
//	  191,   193,   197,   199,   211,   223,   227,   229,   233,   239,   241,   251,   257,   263,
//	  269,   271,   277,   281,   283,   293,   307,   311,   313,   317,   331,   337,   347,   349,
//	  353,   359,   367,   373,   379,   383,   389,   397,   401,   409,   419,   421,   431,   433,
//	  439,   443,   449,   457,   461,   463,   467,   479,   487,   491,   499,   503,   509,   521};

    @Override
    public String getName() {
        return "SmoothRootSieve";
    }

    @Override
    public BigInteger findSingleFactor(BigInteger N) {
        // we only need to factorize over the factor base.
        FactorAlgorithm smallFactoriser = new TDiv31Barrett();
        long n = N.longValue();

        int maxPrime = primes[primes.length -1];
        int maxPrimeLength = length(maxPrime);
        SortedMultiset<BigInteger> factors = null;
        int sqrtNDivS = -1;
        int s = 0;
        int relations = 0;
        while (relations < primes.length) {
            boolean smoothRoot = false;
            while (! smoothRoot)
                {
                s++;
                sqrtNDivS = (int) Math.ceil(Math.sqrt((double)(n/s)));
                // TODO we only need the factors up to a bound. Which?
                factors = smallFactoriser.factor(BigInteger.valueOf(sqrtNDivS));
                int smallFactors = 1;
                for (Map.Entry<BigInteger, Integer> factorEntry: factors.entrySet()) {
                    if (factorEntry.getKey().intValue() <= maxPrime)
                    for (int i = 0; i < factorEntry.getValue(); i++) {
                        smallFactors *= factorEntry.getKey().intValue();
                    }
                }
//                if (smallFactors * smallFactors >= sqrtNDivS)
                if (factors.totalCount() >= 3)
                {
                    // the bigger the size of factors the better they are.
                    // we might calculate sum ( 1/p ), o in factors
                    smoothRoot = true;
                }
            }
            int c = (int) Math.round(Math.sqrt((s * sqrtNDivS * sqrtNDivS - n)/s));
            int sieveIntervalHalf = (int) (3* Math.sqrt(c));
//            int sieveIntervalHalf = (int) (Math.sqrt(c) /2);
            System.out.println("New s : " + s + " = " + factors);
            System.out.println("best dist : " + c);
            System.out.println("sieve interval : " + 2*sieveIntervalHalf );
            //  |                 |  ---            |
            //  0        sieveIntervalHalf  2*sieveIntervalHalf
            // -sieveIntervalHalf 0          sieveIntervalHalf
            int begin = c - sieveIntervalHalf;
            int[] smoothLength = new int[2*sieveIntervalHalf+1];
//        while(search)
            {
                // TODO only half of the primes can hit sqrtNDivS - c and sqrtNDivS + c.
                // but we might need to calculate a square root mod the prime for this
                for (Map.Entry<BigInteger, Integer> factorEntry: factors.entrySet()) {
                    if (factorEntry.getKey().intValue() <= maxPrime) {
                        int prime = 1;
                        for (int i = 0; i < factorEntry.getValue(); i++) {
                            prime *= factorEntry.getKey().intValue();
                            int offsetLower = ((int) Math.ceil((begin + 0.0) / prime)) * prime;
                            int primeLength = length(prime);
                            for (int hit = offsetLower; hit < c + sieveIntervalHalf - 1; hit += prime) {
                                smoothLength[hit - begin] += primeLength;
                                assertTrue(hit % prime == 0);
                            }
                        }
                    }
                }
                for (int i = 0; i < 2 * sieveIntervalHalf - 1; i++) {
                    // we allow 4*maxPrime per factor
                    final double length = length((int) Math.sqrt(sqrtNDivS))- 3;
                    if (smoothLength[i] >= length) {
                        final int cSmooth = begin + i;
                        // TODO we have to build lower * higher * right and only factor this
                        int lower = sqrtNDivS - cSmooth;
                        SortedMultiset<BigInteger> factorsSmall = smallFactoriser.factor(BigInteger.valueOf(lower));
                        final int higher = sqrtNDivS + cSmooth;
                        SortedMultiset<BigInteger> factorsHigh = smallFactoriser.factor(BigInteger.valueOf(higher));
                        System.out.println("at dist " + (higher-sqrtNDivS));
                        System.out.println("found smooth lower " + lower + " = " + factorsSmall);
                        System.out.println(" and  smooth higher" + higher + " = " + factorsHigh);
                        int right = (int) (s * lower * higher - n);
                        SortedMultiset<BigInteger> factorsRight = smallFactoriser.factor(BigInteger.valueOf(Math.abs(right)));
                        System.out.println(" and  smooth right " + right + " = " + factorsRight);
                        System.out.println();
                        if (factorsSmall.getBiggestElement().intValue() <= maxPrime &&
                        factorsHigh.getBiggestElement().intValue() <= maxPrime &&
                        factorsRight.getBiggestElement().intValue() <= maxPrime){
                            relations++;
                            System.out.println("------ new relation ------------");
                        }
                    }
                }
            }

        }
        return null;
    }

    public static int length (int x){
        return 31 - Integer.numberOfLeadingZeros(x);
    }
}
