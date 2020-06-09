package factoring.sieve.triple;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.util.SortedMultiset;
import factoring.math.Column;
import factoring.math.SquareFinder;

import java.math.BigInteger;
import java.util.*;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * This is a variant of the quadratic sieve.
 * It tries to reduce the size of the number(s) to be sieved on.
 * Where the quadratic sieve uses squares x^2 on the left side, subtracts n and then
 * sieve on the right side over y = x^2 - n. The resulting numbers are bigger then n^1/2
 * (by choosing x ~ sqrt(n)).
 * Here we have to sieve over 3 numbers each of size n^(3/8 + epsilon) < n^.376
 * This variant gives up on using squares on the left side. It uses a square s^2 and two
 * smooth parts p_1 ans p_2. Such the resulting number on the right side is s^2 * p-1 * p_2 - n.
 * It splits the numbers x > sqrt(n) from the quadratic sieve in a small 
 * (and such usually smooth) number s and a bigger number p. x = s*p.
 * The aim is to produce a smooth number near n.
 * Like in the quadratic sieve x^2 - n will be of size n^1/2.
 * we use a correction term to reduce the size of the numbers after
 * subtracting n.
 * Instead of x^2 = s^2 * p^2 we use 
 * x(s,i) = s^2*(p-i)*(p+i) = s^2 (p^2 - i^2) = s^2*p^2 - s^2*i^2
 * = x^2 - s^2 * i^2
 * x(s,i) - n = (x^2 - n) - s^2 * i^2
 * we consider x = ceil(sqrt(k*n))
 * We can find different s in the order of n^1/8 by sieving with primes smaller then n^1/8
 * - or by a pre calculated lookup table of primes below n^1/2 -
 * and multiplying the primes together.
 * -> p = x / s
 * for s and p we calculate p-i, p+i, x(s,i) - k*n
 * By choosing s ~ n^1/8 the 3 numbers above have size ~ n^3/8
 * With a pre calculated Lookup table of size a little bit bigger then n^3/8 we can easily iterate over
 * the smooth numbers p-i and p+i (in common if we want) and lookup if
 * x(s,i) - k*n is prime. We choose the lookup table a little bigger then n^3/8 such that
 * there is an i such that a smooth p-i and p+i exist.
 * 
 * We have to ensure to generate only different x(s,i).
 * Different s can lead to the same value x(s,i).
 *
 *  y(s,i,k) = x(s,i) - kn
 *  we choose i_opt such that x(s,i) - kn is minimal
 *  x(s,i) - kn = 0
 *  x^2 - s^2 * i^2 - kn = 0
 *  x^2 - kn  - s^2 * i^2 = 0
 *  i^2  = (x^2 - kn)/s^2
 *  i_opt  = sqrt(x^2 - kn)/s = sqrt(y)/s
 *
 *  = x^2 - n - s^2*(sqrt(x^2 - n)/s +/- t)^2 =
 *  * x^2 - n - s^2*((sqrt(x^2 - n)/s)^2 +/- 2* sqrt(x^2 - n)/s * t +  t^2) =
 *  * x^2 - n -  s^2((x^2-n)/s^2) +/- s^2*(2 * sqrt(x^2 - n)/s * t +  t^2) =
 *  * x^2 - n -  (x^2-n) +/- (2 *s*t* sqrt(x^2 - n)  +  s^2 * t^2 =
 *  * +/- (2+e) *s*t* sqrt(x^2 - n) , when s*t < e*sqrt(x^2 - n)
 * for x=ceil(sqrt(kn)) = sqrt(kn) + e , e < 1
 * sqrt(x^2 - kn)  =
 * sqrt(sqrt(n) + l^2) <=
 *  * sqrt((2+ e')*l)*n^1/4 , when l < e'*sqrt(n)
 *  * ->
 *  * y(s,t,l) < (2+e)*s*t * sqrt((2+ e')*l)*n^1/4  , for s*t*l < e'' * n
 *  * y(s,t,l) < (sqrt(8)+e'')*s*t *l^1/2 *n^1/4
 *  * y(s,t,l) < 3*s*t *l^1/2 *n^1/4
 *  * we want p  = x(s,t,l)
 *  * <-> n^1/2 / s = 3*s*t *l^1/2 *n^1/4
 *  * s^2 = n^1/4 / (3*t * l^1/2)
 *  * s ~ 0.6 * n^1/8 /(t^1/4 + l^1/4)
 * since x(s) - n = x^2 - s^2*i^2 - n = x^2 - n - s^2*i^2
 * for i = sqrt(x^2 - n)/s +/- t
 *
 * y(s,t) < (2+e)*s*t * sqrt(x^2 - n)  , for small t
 * We want y(s,t) < sieveBound
 * (2+e)*s*t * sqrt(x^2 - n)  < sieveBound
 * s > sieveBound / (2t * sqrt(x^2 - n))

 * -> p, y(s,t,l) ~ 1.7 * (t^1/4 + l^1/4) * n^3/8
 * we can choose s smaller as long as all values p + i = n^1/2 / s + sqrt(x^2 - n)/s are precalculated.
 * p + i < 1/s *(n^1/2 + sqrt((2+ e')*l)*n^1/4) <  1.8 * (t^1/4 + l^1/4) * n^3/8
 * y(s,t,l) < smoothBound
 * -> 3*s*t *l^1/2 *n^1/4 < smoothBound
 * s < smoothBound / (3*t *l^1/2 *n^1/4)
 * we also want p < smoothBound
 * Since p = n^1/2 / s < smoothBound
 * s > n^1/2 / smoothBound
 * Since the searchIntervals l and t are subpolynomial all the numbers
 * p-i, p+i and x(i) - n are of size n^(3/8 + epsilon), were epsilon can be as small as we want for some
 * n. Since all 3 numbers are a polynomial over i, we can sieve over all of them and pick i where
 * the sum of the length of the sieved factors for all of the values are below the threshold, and
 * determine the factors.
 * This version of the algorithm is based on long values (i.e. n < 2^63).
 * We will precalculate all smooth values below n^(3/8 + epsilon) in a BitSet.
 * So we can easily determine if the 3 numbers are smooth. Additionally we can easily get the next
 * sooth value p-i' from p-i. t will be small to keep the resulting number x(i) - n small.
 * By choosing s, p > n^3/8 , p < n^(3/8 + e)
 * for x = sqrt(n) + l , were l is a small search interval < n^e' for any e'
 * i = sqrt(x^2 - n)/s +/- t =
 * sqrt((sqrt(n) + l)^2 - n)/s +/- t =
 * sqrt(2*sqrt(n)*l + l^2)/s +/- t =
 * 2*sqrt(l) * n^1/4/s +/- t  for some n
 * -> i < 3*sqrt(l) * n^1/8
 * for p +/- i there are O(n^(1/8+e)) choices. They will be precalculated.
 */
public class TripleLookupAllSieve extends FactorAlgorithm {

    static int [] primes = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
            103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
            199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
            313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
            433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
            563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
            673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
            811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
            941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,
            1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,
            1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,
            1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,
            1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,
            1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,
            1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,
            1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
            1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,
            1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063/*,
			2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,
			2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,
			2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,
			2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,
			2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,
			2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,
			2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,
			2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,
			3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,
			3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,
			3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,
			3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,
			3527,3529,3533,3539,3541,3547,3557,3559,3571,3581,3583,3593,3607,3613,3617,
			3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733,
			3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,
			3877,3881,3889,3907,3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,
			4003,4007,4013,4019,4021,4027,4049,4051,4057,4073,4079,4091,4093,4099,4111,
			4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,4241,
			4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,
			4373,4391,4397,4409,4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,4507,
			4513,4517,4519,4523,4547,4549,4561,4567,4583,4591,4597,4603,4621,4637,4639,
			4643,4649,4651,4657,4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,4759,
			4783,4787,4789,4793,4799,4801,4813,4817,4831,4861,4871,4877,4889,4903,4909,
			4919,4931,4933,4937,4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,5009,
			5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101,5107,5113,5119,5147,
			5153,5167,5171,5179,5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,5281,
			5297,5303,5309,5323,5333,5347,5351,5381,5387,5393,5399,5407,5413,5417,5419,
			5431,5437,5441,5443,5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,5527,
			5531,5557,5563,5569,5573,5581,5591,5623,5639,5641,5647,5651,5653,5657,5659,
			5669,5683,5689,5693,5701,5711,5717,5737,5741,5743,5749,5779,5783,5791,5801,
			5807,5813,5821,5827,5839,5843,5849,5851,5857,5861,5867,5869,5879,5881,5897,
			5903,5923,5927,5939,5953,5981,5987,6007,6011,6029,6037,6043,6047,6053,6067,
			6073,6079,6089,6091,6101,6113,6121,6131,6133,6143,6151,6163,6173,6197,6199,
			6203,6211,6217,6221,6229,6247,6257,6263,6269,6271,6277,6287,6299,6301,6311,
			6317,6323,6329,6337,6343,6353,6359,6361,6367,6373,6379,6389,6397,6421,6427,
			6449,6451,6469,6473,6481,6491,6521,6529,6547,6551,6553,6563,6569,6571,6577,
			6581,6599,6607,6619,6637,6653,6659,6661,6673,6679,6689,6691,6701,6703,6709,
			6719,6733,6737,6761,6763,6779,6781,6791,6793,6803,6823,6827,6829,6833,6841,
			6857,6863,6869,6871,6883,6899,6907,6911,6917,6947,6949,6959,6961,6967,6971,
			6977,6983,6991,6997,7001,7013,7019,7027,7039,7043,7057,7069,7079,7103,7109,
			7121,7127,7129,7151,7159,7177,7187,7193,7207,7211,7213,7219,7229,7237,7243,
			7247,7253,7283,7297,7307,7309,7321,7331,7333,7349,7351,7369,7393,7411,7417,
			7433,7451,7457,7459,7477,7481,7487,7489,7499,7507,7517,7523,7529,7537,7541,
			7547,7549,7559,7561,7573,7577,7583,7589,7591,7603,7607,7621,7639,7643,7649,
			7669,7673,7681,7687,7691,7699,7703,7717,7723,7727,7741,7753,7757,7759,7789,
			7793,7817,7823,7829,7841,7853,7867,7873,7877,7879,7883,7901,7907,7919,7927,
			7933,7937,7949,7951,7963,7993,8009,8011,8017,8039,8053,8059,8069,8081,8087,
			8089,8093,8101,8111,8117,8123,8147,8161,8167,8171,8179,8191,8209,8219,8221,
			8231,8233,8237,8243,8263,8269,8273,8287,8291,8293,8297,8311,8317,8329,8353,
			8363,8369,8377,8387,8389,8419,8423,8429,8431,8443,8447,8461,8467,8501,8513,
			8521,8527,8537,8539,8543,8563,8573,8581,8597,8599,8609,8623,8627,8629,8641,
			8647,8663,8669,8677,8681,8689,8693,8699,8707,8713,8719,8731,8737,8741,8747,
			8753,8761,8779,8783,8803,8807,8819,8821,8831,8837,8839,8849,8861,8863,8867,
			8887,8893,8923,8929,8933,8941,8951,8963,8969,8971,8999,9001,9007,9011,9013,
			9029,9041,9043,9049,9059,9067,9091,9103,9109,9127,9133,9137,9151,9157,9161,
			9173,9181,9187,9199,9203,9209,9221,9227,9239,9241,9257,9277,9281,9283,9293,
			9311,9319,9323,9337,9341,9343,9349,9371,9377,9391,9397,9403,9413,9419,9421,
			9431,9433,9437,9439,9461,9463,9467,9473,9479,9491,9497,9511,9521,9533,9539,
			9547,9551,9587,9601,9613,9619,9623,9629,9631,9643,9649,9661,9677,9679,9689,
			9697,9719,9721,9733,9739,9743,9749,9767,9769,9781,9787,9791,9803,9811,9817,
			9829,9833,9839,9851,9857,9859,9871,9883,9887,9901,9907,9923,9929,9931,9941,
			9949,9967,9973,10007,10009,10037,10039,10061,10067,10069,10079,10091,10093*/};
    int smoothBound;
    BitSet smoothNumber;
    short[] maxFactor;
    double smoothExpLower = .1;
    double smoothExpHigher = .15;
    long maxN = 1l << 50;
    short[][] smoothPI;
    private int smoothBits;
    int bitsN = Integer.MAX_VALUE;

    public TripleLookupAllSieve(int bitsN) {
        this.bitsN = bitsN;
    }


    public long findSingleFactor(long n) {
        initSmoothNumbers();
        int factorBaseSize = (int) factorBaseSize(n);
        int factorBase []= Arrays.copyOf(primes, factorBaseSize);
        int factorBaseMax = factorBase [factorBaseSize-1];
//        final TDiv63Inverse tdiv = new TDiv63Inverse(factorBase[factorBase.length -1]);
        int sieveInterval = factorBase[factorBase.length -1]*4;
//        int smoothBound = (int) Math.pow(n, 0.5);
//        int[] numberLengthSqrtN = new int [sieveInterval];
        int relations = 0;
        long xi = 0;
        Integer[] factorBaseArr = IntStream.of(factorBase).boxed().toArray(Integer[]::new);
        SquareFinder finder = new SquareFinder(Arrays.asList(factorBaseArr), n);
        for (int k=1; relations < factorBaseSize; k++) {
            System.out.println("Next k : " + k);
            double sqrtN = sqrt(k*n);
            int x = (int) Math.ceil(sqrtN);
            int nPowOneDivEight = (int) Math.pow(n, .125);
            System.out.println("n^1/8 : " + nPowOneDivEight);
            Set<Integer> smooth = new HashSet<>();
            int searchIntervalP = 10;
//             this ensures that y is lower then smoothBond in the hole searchIntervalP
            double sBest = smoothBound / (2 * searchIntervalP * sqrt(x * x - n));
            // sLower s -> higher p - less smooth values, bigger range for i
            int sHigher = Math.max((int) (sBest * 2), 2);
            int[] factorList = getPrimeFactors(x, sHigher);
//                    int[] factorArray = IntStream.of(factorListRaw).filter(i -> i >sLower).toArray(new int []);
            System.out.println(" s : " + factorList + " y : " + (x*x - n));
//                    System.out.println(" s : " + factorArray + " y : " + (x*x - n));
            int operations;

            for (int l=0; factorList[l] > 1; l++){
                System.out.print(factorList[l] + ",");
//                    for (int s : smooth) {
                List smoothY = new ArrayList();
                operations = 0;
                int relationsPerS = 0;
                int s = factorList[l];
                long p = (int) (x / s);
                int y = (int) (s * s * p * p - n);
                long yDivSSquare = y / (s * s);
                if (y == 0){
                    relations++;
                    relationsPerS++;
                    List<Integer> pFactors = getFactorList((int) p, sHigher);
                    Multiset<Integer> factorsLeft = HashMultiset.create();
                    factorsLeft.add(s);
                    factorsLeft.addAll(pFactors);
                    assertEquals(s*s*p*p, multiply(factorsLeft));
                    finder.addFactors(factorsLeft, HashMultiset.create());

                } else {
                    long iOpt = (long) sqrt(yDivSSquare);
                    int yMax = 2 * (int) Math.pow(n, .375);
                    long iRange = iOpt == 0 ? 0 : (yMax - yDivSSquare) / (s * s * 2 * iOpt);
//                        iRange /= 4;
                    iRange = iRange < 0 ? 1 : iRange;
                    long i = iOpt - iRange;
                    i = i < 0 ? 0 : i;
                    i = nextSmoothI(factorBaseMax, p, i);
                    while (i < iOpt + iRange) {
                        int pMinI = (int) (p - i);
                        if (smoothNumber.get(pMinI) && maxFactor[pMinI] < factorBaseMax) {
                            xi = s * s * (p - i) * (p + i);
                            y = Math.abs((int) (xi - n));
                            if (smoothNumber.get(y) && maxFactor[y] < factorBaseMax) {
                                smoothY.add(y);
                                relations++;
                                relationsPerS++;
                                List<Integer> pMinIFactors = getFactorList(pMinI, sHigher);
                                List<Integer> pPlusIFactors = getFactorList((int) (p + i), sHigher);
                                List<Integer> yFactors = getFactorList(y, sHigher);
                                Multiset<Integer> factorsLeft = HashMultiset.create();
                                factorsLeft.add(s);
                                factorsLeft.addAll(pMinIFactors);
                                factorsLeft.addAll(pPlusIFactors);
                                Multiset<Integer> factorsRight = HashMultiset.create(yFactors);
                                assertEquals(s*s*(p+i)*(p-i), multiply(factorsLeft));
                                assertEquals(y, multiply(factorsRight));
                                finder.addFactors(factorsLeft, factorsRight);

                            }
                        }
                        operations++;
//                            iIndex++;
//                            i = smoothI[iIndex];
                        i = nextSmoothI(factorBaseMax, p, i + 1);
                    }
                }
                System.out.println(" s : " + s + " relations : " + relationsPerS + " operations : " + operations + " rate : " + (0.0 + relationsPerS)/operations + " y's : " + smoothY);
            }
        }
        List<Column> smoothMatrix = finder.initMatrix();
        List<Column> reducedMatrix = finder.reduceMatrix(smoothMatrix);
        do {
            smoothMatrix = reducedMatrix;
            reducedMatrix = finder.reduceMatrix(smoothMatrix);
        }
        while (reducedMatrix.size() < smoothMatrix.size());
        long factor = finder.doGaussElimination(reducedMatrix);

        return factor;
    }

    private long multiply(Multiset<Integer> factorsLeft) {
        return factorsLeft.elementSet().stream().reduce(1, (a,b)-> a*b);
    }

    private double factorBaseSize(long n) {
        return .5 * sqrt(exp(sqrt(log(n) * log(log(n)))));
    }

    private long nextSmoothI(int factorBaseMax, long p, long i) {
        int pPlusI = (int) (p + i);
        // the smooth numbers are defined over a bigger factor base
        pPlusI = smoothNumber.nextSetBit(pPlusI);
        while (maxFactor[pPlusI] > factorBaseMax){
            pPlusI = smoothNumber.nextSetBit(pPlusI + 1);
        }
        i = pPlusI - p;
        return i;
    }

    private int[] getPrimeFactors(int xDiv, int upperBound) {
        int[] factors = new int[smoothBits];
        int i = 0;
        int factor = Integer.MAX_VALUE;
        while (maxFactor[xDiv] > 1){
            factor = maxFactor[xDiv];
            xDiv = xDiv / factor;
            // add only different factors.
            if ((i == 0 || factors[i-1] != factor) && factor <= upperBound )
                factors[i++] = (short) factor;
        }
        if (factors[0] == 0) factors [0] = factor;
        return factors;
    }

    private List<Integer> getPrimeFactorList(int xDiv, int upperBound) {
        List<Integer> factors = new ArrayList<>();
        while (maxFactor[xDiv] > 1){
            int factor = maxFactor[xDiv];
            xDiv = xDiv / factor;
            if (factor <= upperBound)
                factors.add(factor);
        }
        return factors;
    }

    private List<Integer> getFactorList(int xDiv, int upperBound) {
        List<Integer> primeFactors = getPrimeFactorList(xDiv, upperBound);
        
        powerSet(primeFactors, 0, new ArrayList<>(), upperBound, 1);      
    }
    public static void powerSet(List<Integer> numbers, int index, ArrayList<Integer> currentFactors, ArrayList<ArrayList<Integer>> result, int upperBound, int prod) {
		if (index == numbers.size() && prod < upperBound)
			result.add(currentFactors);
		else {
			powerSet(numbers, index + 1, currentFactors, result, upperBound, prod); //do not take the current number
			ArrayList<Integer> withNewFactor = new ArrayList<>(currentFactors);
			withNewFactor.add(numbers.get(index));
			powerSet(numbers, index + 1, withNewFactor, result, upperBound, prod * numbers.get(index)); //take it
		}
	}

    private void findFactors(double lower, double higher, List<Integer> factorX, int prod, Set<Integer> shifts, int index) {
        if (lower < prod && prod < higher) {
            shifts.add(prod);
        }
        if (index == factorX.size())
            return;
        findFactors(lower, higher, factorX, prod, shifts, index + 1);
        findFactors(lower, higher, factorX, prod * factorX.get(index), shifts, index + 1);
    }

    private void initSmoothNumbers() {
        long n = 1l << 42;
        // 23 bits is ~ 2 MB -> in best case usable for 61 Bit numbers
        // 34 bits is ~ 4 GB -> in best case for 90 Bit numbers
//        final int smoothBits = 28;
        smoothBits = 24;
        if (bitsN/2 + 4 < smoothBits)
            smoothBits = bitsN/2 + 2;
//        smoothBits = 21;
        smoothBound = 1 << smoothBits;
//        smoothBound = Integer.MAX_VALUE - 1;
        smoothNumber = new BitSet(smoothBound);
        // this is only temporary we might do it in different runs.
//        byte[] numberLength = new byte [smoothBound];
        byte[] numberLength = new byte [smoothBound];
        int factorBaseSize = (int) (factorBaseSize(n));
        int factorBase []= Arrays.copyOf(primes, factorBaseSize);
        // a power of 2, a little bit bigger then the smoothBound

        for (int i = 0; i < factorBaseSize; i++) {
            int primeLength = 32 - Integer.numberOfLeadingZeros(factorBase[i] - 1);
            for (int pos = factorBase[i] ; pos < numberLength.length; pos += factorBase[i]) {
                numberLength[pos] -= primeLength;
                assertTrue((pos) % factorBase[i] == 0);
            }
        }
        int lengthThreshold;
        SortedMultiset<BigInteger> minFactorXLower;
        for (int j=smoothBound-1; j > 0 ; j--){
            lengthThreshold = 32 - Integer.numberOfLeadingZeros(j-1);
            int soothLength = (int) (-lengthThreshold) + 8;
            if (numberLength[j] <= soothLength) {
                smoothNumber.set(j);
            }
        }
        System.out.println("smooth numbers : " + smoothNumber.cardinality());
        System.out.println("ratio : " + ((double)smoothBound)/ smoothNumber.cardinality());
        // around the same size as smooth numbers BitSet
        maxFactor = new short [(int) smoothBound];
        for (int i = 0; i < factorBaseSize; i++) {
            for (int pos = factorBase[i] ; pos < numberLength.length; pos += factorBase[i]) {
                // we only need the factors for numbers which have at least n^1/8
//                if (smoothNumber.get((int) pos))
                {
                    maxFactor[pos] = (short) factorBase[i];
                    assertTrue((pos) % factorBase[i] == 0);
                }
            }
        }
        System.out.println("initialization ready ");
    }

    private void initPMinPlusI(long n) {
        double nPow3Div8 = pow(n, .375);
        int nPow1Div8 = (int) pow(n, .125);
        int sqrtSearchIntervalX = 20;
        smoothPI = new short[smoothBound][];
        short[] smoothPIndex = new short[smoothBound];
//        int p=1; // we might start even later
        int pMinI= (int) smoothBound - 2*nPow1Div8*sqrtSearchIntervalX; // we might start even later
        pMinI = smoothNumber.previousSetBit(pMinI);
        while (pMinI > 0){
            int iHalf = 0;
            while (iHalf < nPow1Div8*sqrtSearchIntervalX){
                if (smoothNumber.get(pMinI + 2*iHalf)){
                    if (smoothPI[pMinI+iHalf] == null) {
                        smoothPI[pMinI + iHalf] = new short[nPow1Div8 * sqrtSearchIntervalX];
                        for (int i = 0; i< smoothPI[pMinI + iHalf].length; i++) {
                            smoothPI[pMinI + iHalf][i] = Short.MAX_VALUE;
                        }
                        smoothPIndex[pMinI + iHalf] = 0;
                    }
                    smoothPI[pMinI+iHalf][smoothPIndex[pMinI+iHalf]] = (short) iHalf;
                    smoothPIndex[pMinI+iHalf] = (short) (smoothPIndex[pMinI+iHalf]+1);
                }
                iHalf = smoothNumber.nextSetBit(pMinI + 2*iHalf +1) - (pMinI + iHalf);
            }
            pMinI = smoothNumber.previousSetBit(pMinI-1);
        }
    }

    @Override
    public String getName() {
        return "TripleLookupSieve";
    }

    @Override
    public BigInteger findSingleFactor(BigInteger N) {
        return BigInteger.valueOf(findSingleFactor(N.longValue()));
    }

}
