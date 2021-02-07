package factoring.sieve.triple;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.SortedMultiset_BottomUp;
import factoring.math.Column;
import factoring.math.SquareFinder;
import factoring.trial.TDiv31Barrett;

import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static java.lang.Math.log;
import static org.junit.Assert.assertTrue;

/**
 * Starting with the square t^2 , t = ceil(sqrt(n/s)) we can build (t-c)*(t+c) = t^2 - c^2.
 * With c^2 we can reduce the size of s*(t-c)*(t+c) - n to O(n^1/4). But t-c and t+c are no square and do not have to be smooth.
 * We will have a bit representation of smooth numbers x and xMax - x up to xMax = sqrt(n), we can calculate smooth
 * number t-c and t+c by a small copy and an and of the smooth numbers.
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

public class SmoothNumbers extends FactorAlgorithm {

    private static final String SMOOTH_FILE_NAME = "smoothNumbers";
    static int [] primes = {-1, 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
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

    // for a range we hash the smooth numbers.
    // TODO only use the hash to reduce memory.
    // instead of 2* n^1/2 bits we only need n^1/2 / sqrt(exp(sqrt(log(n) * log(log(n)))))
    private BitSet smoothNegativeHash;
    private BitSet smoothNumber;
    int[] maxSmoothFactorIndex;

    Multimap<Integer,Integer> numbersForDist = HashMultimap.create();

    private int smoothBits;
    private int bitsN;
    private int smoothBound;

    private short[] maxFactor;

    TDiv31Barrett smallFactoriser = new TDiv31Barrett();

    private static int steps = 0;

    @Override
    public String getName() {
        return "Smooth numbers ";
    }

    @Override
    public BigInteger findSingleFactor(BigInteger N) {
        // TODO we need a static initializer here.
        initSmoothNumbers();
        // we only need to factorize over the factor base.
        long n = N.longValue();

//        int maxPrimeLength = length(maxPrime);
        long sqrtNDivS = -1;
        long s = 0;
        Set<Integer> relations = new HashSet<>();
        final int factorBaseSize = (int) factorBaseSize((long) n);
        System.out.println("factorbase size : " + factorBaseSize);
        int maxPrime = primes[factorBaseSize -1];
        System.out.println("max factor : " + maxPrime);
        int lastRelationSize = 0;
        int k = 0;
        Integer[] factorBaseArr = IntStream.of(primes).limit(factorBaseSize).boxed().toArray(Integer[]::new);
        SquareFinder finder = new SquareFinder(Arrays.asList(factorBaseArr), n);
        List<Double> fractions = new ArrayList<>();
        int allTries = 0;
        while (relations.size() < 2* factorBaseSize) {
            k++;
            s = k+1;
            while (relations.size() < 2*factorBaseSize && s > 1) {
                s--;
                // we consider s * (t^2 - c^2) - kn = r for many k,s to keep r small.
                Double fraction = (double)k / s;
                int index = Collections.binarySearch(fractions, fraction);
                index = index < 0 ? -index -1 : index;
                if (index >= fractions.size() || Math.abs(fractions.get(index) - fraction) > 0.001) {
                    fractions.add(fractions.size() == 0 ? 0 : index, fraction);
                    allTries += findRelations(k, n, s, relations, factorBaseSize, maxPrime, lastRelationSize, finder);
                    lastRelationSize = relations.size();
                    if (k != s) {
                        allTries += findRelations((int) s, n, k, relations, factorBaseSize, maxPrime, lastRelationSize, finder);
                        lastRelationSize = relations.size();
                    }
                }
            }
        }
        System.out.println("relations considered : " + allTries);
        List<Column> smoothMatrix = finder.initMatrix();
        List<Column> reducedMatrix = finder.reduceMatrix(smoothMatrix);
        do {
            smoothMatrix = reducedMatrix;
            reducedMatrix = finder.reduceMatrix(smoothMatrix);
        }
        while (reducedMatrix.size() < smoothMatrix.size());
        long factor = finder.doGaussElimination(reducedMatrix);

        System.out.println("Number of relations considered overall : " + steps);
        return BigInteger.valueOf(factor);
    }

    private int findRelations(int k, long n, long s, Set<Integer> relations, int factorBaseSize, int maxPrime, int lastRelationSize, SquareFinder finder) {
        n = k*n;
        long sqrtNDivS = (long) Math.ceil(Math.sqrt((double) (n / s)));
        long tInterval = 100;
        int tries = 0;
        long t = sqrtNDivS;
        for (long tDiff = 0; tDiff < tInterval; tDiff++, t++) {
            // TODO we only need the factors up to a bound. Which?
            // c is the minimal c for |(t^2 - c^2) - n|
            double c = sqrt((s * t * t - n) / s);
            // (t^2 - (c+l)^2) - n > -n^1/2
            // (t^2 - (c^2+2cl+l^2)) - n > -n^1/2
            //  2cl+l^2 < n^1/2
            //  l < n^1/4 and
            //  l < n^1/2 / 2c
//                int sieveIntervalHalf = (int) (3 * Math.sqrt(c));
            // we need two smooth numbers over the factor base
//            double sieveIntervalHalf = c/ (tDiff+1);
            double sieveIntervalHalf1 = Math.pow(n, .25);
            double sieveIntervalHalf2 = Math.sqrt(n) / (2*c);
            double smoothFactor = 4;
            double sieveIntervalHalf = min(smoothFactor * min(sieveIntervalHalf1, sieveIntervalHalf2), c);
//            double sieveIntervalHalf = factorBaseSize * factorBaseSize / 2;
//            if (sieveIntervalHalf > c)
//                sieveIntervalHalf = c;
//            int sieveIntervalHalf = (int) (Math.sqrt(c) /2);
            //  |                 |  ---            |
            //  0        sieveIntervalHalf  2*sieveIntervalHalf
            // -sieveIntervalHalf 0          sieveIntervalHalf
//    int beginDist = c - sieveIntervalHalf - 1;
//            int endDist = (int) Math.ceil(sieveIntervalHalf) + 1;
//        int tryDist = beginDist;
            int beginPos = (int) (t + c - sieveIntervalHalf);
            // TODO use double c for more precision
            // TODO upper
            List<Integer> smoothPos = smoothPosPre((int) t, (int)Math.round(c), (int)sieveIntervalHalf);
            for (int tryDist : smoothPos) {
//            while ((tryDist = nextSmoothPos((int) sqrtNDivS, tryDist, endDist)) >= 0) {
                if (smoothPos.size() < tries + 1 || smoothPos.get(tries) != tryDist) {
                    System.out.println();
                }
                tries++;
                steps++;
                final long lower = t - tryDist;
                final long higher = t + tryDist;
                int right = (int) (s * lower * higher - n);
                boolean print = true;
                if (print) {
                    SortedMultiset<BigInteger> factorsRight = new SortedMultiset_BottomUp<>();
                    checkSmooth(factorBaseSize, right, factorsRight);
                    SortedMultiset<BigInteger> factorsSmall = new SortedMultiset_BottomUp<>();
                    checkSmooth(factorBaseSize, (int) lower, factorsSmall);
                    SortedMultiset<BigInteger> factorsHigh = new SortedMultiset_BottomUp<>();
                    checkSmooth(factorBaseSize, (int) higher, factorsHigh);
                    System.out.println("try at dist " + (higher - t));
                    System.out.println("tries " + tries);
                    System.out.println("found smooth lower " + lower + " = " + factorsSmall);
                    System.out.println(" and  smooth higher" + higher + " = " + factorsHigh);
                    //                        SortedMultiset<BigInteger> factorsRight = smallFactoriser.factor(BigInteger.valueOf(Math.abs(right)));
                    System.out.println(" and  smooth right " + right + " = " + factorsRight);
                    System.out.println();
                }
//            boolean isRightSmooth1 = smallFactoriser.factor(right, factorBaseSize, factorsRight);
                // here we need the full range of factors!?
                boolean isRightSmooth2 = maxSmoothFactorIndex[Math.abs(right)] < factorBaseSize;
                boolean isRightSmooth3 = smoothNumber.get(Math.abs(right));
//            if ( (isRightSmooth1 != isRightSmooth2) || (isRightSmooth2 != isRightSmooth3)) {
//                System.out.println();
//            }
//            Assert.assertTrue(smoothNumber.get(right) == isRightSmooth);
                if (isRightSmooth2) {
                    addRelation(relations, finder, s, (int) lower, (int) higher, right, factorBaseSize, n, k);
                }
            }
            final int hits = relations.size() - lastRelationSize;
            int hitRate = hits == 0 ? Integer.MAX_VALUE : (tries / hits);
            SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
            smallFactoriser.factor((int) t, factorBaseSize, factors);
            System.out.println("k " + k + " s : " + s + "  t diff : " + tDiff + " : " + factors);
            System.out.println("best dist : " + c);
            System.out.println("sieve interval : " + 2 * sieveIntervalHalf);
            System.out.println("rel  " + hits + " for s : " + s + " hit rate 1 : " + hitRate);
            System.out.println("work done : " + relations.size() / (2.0 * factorBaseSize) + "%");
            System.out.println();
        }
        return tries;
    }

    private void addRelation(Set<Integer> relations, SquareFinder finder, long s, int lower, int higher, int right, int factorBaseSize, long n, int k) {
        // ensure that this relation is not a duplicate, and all numbers fit in the matrix
        // right should not be dividable by k
        // TODO if we can divide lower or higher by k we can still use the relation
        // TODO check deeper if we no already have such a relation e.g. adding removing some small factors.
        if (relations.contains(right) || k > 1 && ((right / k) * k == right) )
            return;
        // can only happen if k is high -> the generated numbers might be bigger then sqrt(n)
        if (maxSmoothFactorIndex[Math.abs(lower)] >= factorBaseSize) {
            System.out.println("lower number 't-c'  has a factor " + primes[maxSmoothFactorIndex[Math.abs(lower)]] + " outside the factor base.");
            return;
        }
        if (maxSmoothFactorIndex[Math.abs(higher)] >= factorBaseSize) {
            System.out.println("higher number 't+c' has a factor " + primes[maxSmoothFactorIndex[Math.abs(higher)]] + " outside the factor base.");
            return;
        }

        relations.add(right);
        Multiset<Integer> leftFactors = HashMultiset.create();
        Multiset<Integer> rightFactors = HashMultiset.create();
        addFactors((int) s, leftFactors);
        addFactors(lower, leftFactors);
        addFactors(higher, leftFactors);
        addFactors(right, rightFactors);
        int prodLeft = leftFactors.stream().reduce(1, (a,b) -> a*b);
        int prodRight = rightFactors.stream().reduce(1, (a,b) -> a*b);
        if (prodLeft % n != (prodRight +n) %n)
            System.out.println(prodLeft + " % " + n + " = " + prodRight + " but should be " +  prodLeft % n);
        finder.addFactors(leftFactors, rightFactors);
        System.out.println("------ new relation ------------");

    }

    private void checkSmooth(int factorBaseSize, int number, SortedMultiset<BigInteger> factors) {
        boolean isRightSmooth1 = smallFactoriser.factor(number, factorBaseSize, factors);
        boolean isRightSmooth2 = maxSmoothFactorIndex[Math.abs(number)] < factorBaseSize;
        boolean isRightSmooth3 = smoothNumber.get(Math.abs(number));
        if ( (isRightSmooth1 != isRightSmooth2) || (isRightSmooth2 != isRightSmooth3)) {
            System.out.println(" wrong factors : " + factors + " for number "+ number + " factor outside : " + primes[maxSmoothFactorIndex[Math.abs(number)]]);
            return;
        }
    }

    private void addFactors(int number, Multiset<Integer> factors) {
        if (number < 0)
            factors.add(-1);
        int x = Math.abs(number);
        while (x != 1){
            int factor = primes[maxSmoothFactorIndex[x]];
            factors.add(factor);
            // TODO replace division
            x = x / factor;
        }
    }

    /**
     * Gives back the minimal position c > dist where
     * pos + c and pos - c are smooth over a factor base of size
     * sqrt(exp(sqrt(2*log(pos) * log(2*log(pos))))) / 2.
     * and dist < = maxDist.
     * This is the most expensive step in this algorithm.
     * If we store for each pos the smooth pos-2c with 2c < sqrt(pos) in order
     * The we start at smooth pos+c looking for numbers c
     * If we take a hash h (pos) of the smooth numbers Mod m =2^k for the numbers pos with a
     * range sqrt(pos) then (pos+dist) = (pos-dist)
     *
     * 15 = 10 + 5 -> smooth
     *  5 = 10 - 5 -> smooth
     *
     * @param pos
     * @param dist
     * @return
     */
    public int nextSmoothPos(int pos, int dist, int maxDist){
        // we need to store this
        int tryDist = dist;
        while ((tryDist = nextSmoothPosDist(pos, tryDist)) <= maxDist){
            int pMinI = (int) (pos - tryDist);
            if (smoothNumber.get(pMinI)) {
                return tryDist;
            }
        }
        return -1;
    }
    // TODO delete
    public List<Integer> smoothPos(int pos, int dist){
//        final double smoothDist = factorBaseSize(pos) / 8;
        final double avgSmoothDist = .5;
        int mod = 1 << (31 - Integer.numberOfLeadingZeros((int) (dist / avgSmoothDist)));
        BitSet smoothLowerDistHash = new BitSet(mod);
        BitSet smoothUpperDistHash = new BitSet(mod);
        Multimap<Integer,Integer> numbersLowerDist = HashMultimap.create();

        int smoothPos = pos - 1;
        while ((smoothPos = smoothNumber.nextSetBit((int) smoothPos + 1)) <= pos + dist){
            int hashUpperDist = (smoothPos-pos) % mod;
            smoothUpperDistHash.set(hashUpperDist);
            numbersLowerDist.put(hashUpperDist, smoothPos);
        }
        smoothPos = pos -dist -1;
        while ((smoothPos = smoothNumber.nextSetBit((int) smoothPos + 1)) <= pos){
            int hashLowerDist = (pos - smoothPos) % mod;
            smoothLowerDistHash.set(hashLowerDist);
        }
        smoothLowerDistHash.and(smoothUpperDistHash);
        List<Integer> bothSmooth = new ArrayList<>();
//        while ((tryPos = smoothLowerHash.nextSetBit((int) pos)) <= dist) {
        int misses = 0;
        for (int smoothDist = smoothLowerDistHash.nextSetBit(0); smoothDist >= 0; smoothDist = smoothLowerDistHash.nextSetBit(smoothDist+1)) {
            for (int lowerNumber : numbersLowerDist.get(smoothDist)) {
                if (smoothNumber.get(lowerNumber) && smoothNumber.get(2 * pos - lowerNumber)) {
                    bothSmooth.add(lowerNumber - pos);
                }
                else{
                    misses++;
                }
                if (smoothDist == Integer.MAX_VALUE) {
                    break; // or (i+1) would overflow
                }
            }
        }
        return bothSmooth;
    }

    /**
     * Gives back the c+s  with t + c + s is smooth and t - c - s, -dist < s < dist is smooth.
     * We take the Binary representation of smooth numbers pos + i and
     * the Binary representation of smooth numbers pos + i.
     * Do an and, and iterate over the results. Since getting the next set bit is fast
     * when the dist is less then 64 (long bit length).
     * @param t
     * @param c
     * @param dist
     * @return
     */
    public List<Integer> smoothPosPre(int t, int c, int dist){
//        final double smoothDist = factorBaseSize(pos) / 8;
        // TODO dist should bew a multiple of 64 to improve performance of BitSet.get
        BitSet smoothPosDistHash = smoothNumber.get(t + c - dist, t + (c + dist));
        BitSet smoothNegDistHash = smoothNegativeHash.get(smoothBound-(t-c+dist), smoothBound - (t- c-dist));
        smoothPosDistHash.and(smoothNegDistHash);
        List<Integer> bothSmooth = new ArrayList<>();
//        while ((tryPos = smoothLowerHash.nextSetBit((int) pos)) <= dist) {
        int misses = 0;
        // TODO directly use the iterator; do not use the List
        for (int smoothDist = smoothPosDistHash.nextSetBit(0); smoothDist >= 0; smoothDist = smoothPosDistHash.nextSetBit(smoothDist+1)) {
            int higherNumber = t + c - dist + smoothDist;
                // since we do not use a hash at the moment this is always true.
                if (smoothNumber.get(higherNumber) && smoothNumber.get(2 * t - higherNumber)) {
                    bothSmooth.add(higherNumber-t);
                }
                else{
                    misses++;
                }
                if (smoothDist == Integer.MAX_VALUE) {
                    break; // or (i+1) would overflow
                }
        }
        return bothSmooth;
    }

    private int nextSmoothPosDist(int pos, int dist) {
        // the smooth numbers are defined over a bigger factor base
        int nestPos = smoothNumber.nextSetBit(pos + dist + 1);
        return nestPos - pos;
    }

    private void initSmoothNumbers() {
        long start = System.currentTimeMillis();
        smoothBits = 27;
        long n = 1l << smoothBits;
        // 23 bits is ~ 2 MB
        // 34 bits is ~ 4 GB -> in best case for 90 Bit numbers
//        final int smoothBits = 28;
//        smoothBits = 15;
//        if (bitsN/2 + 4 < smoothBits)
//            smoothBits = bitsN/2 + 2;
//        smoothBits = 21;
        smoothBound = 1 << smoothBits;
//        smoothBound = Integer.MAX_VALUE - 1;
        // this is only temporary we might do it in different runs.
//        byte[] numberLength = new byte [smoothBound];
//        byte[] numberLength = new byte [smoothBound];
        int factorBaseSize = (int) (factorBaseSize(n));
//        int factorBase []= Arrays.copyOf(primes, factorBaseSize);
        // a power of 2, a little bit bigger then the smoothBound

        boolean load = true;
        if (load) {
            // TODO just serialize the BitSet
    //        File smoothNumbersFile = new File(SMOOTH_FILE_NAME + smoothBits + ".txt");
            smoothNumber = (BitSet) readJavaObject(SMOOTH_FILE_NAME + smoothBits + ".dat");
            maxSmoothFactorIndex = (int[]) readJavaObject("smoothFactors" + smoothBits + ".dat");

            if (smoothNumber != null){
                long end = System.currentTimeMillis();
                System.out.println("read in file with smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");

                smoothNegativeHash = new BitSet(smoothBound);
                smoothNumber.stream().forEach(i -> smoothNegativeHash.set(smoothBound - i));
                return;
            }
        }

//        if (smoothNumbersFile.exists()){
//            try (BufferedReader reader = new BufferedReader(new FileReader(smoothNumbersFile))) {
//                String numberString;
//                while ((numberString= reader.readLine()) != null){
//                    final Integer number = Integer.valueOf(numberString);
//                    smoothNumber.set(number);
//                    smoothNegativeHash.set(smoothBound - number);
//                }
//            }
//            catch (IOException e){
//            }
//            System.out.println("found smooth number file for " + smoothBits + " bits. Will use it");
//            return;
//        }
//        try (PrintWriter pw = new PrintWriter(smoothNumbersFile)) {
//            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(smoothNumbersFile)), false)) {

//        for (int j=290; j > 0 ; j--){
        smoothNumber = new BitSet(smoothBound);
//        smoothPosHash = new BitSet(smoothBound);
        maxSmoothFactorIndex = new int[smoothBound];
        smoothNumber.set(1);
        Arrays.fill(maxSmoothFactorIndex, Integer.MAX_VALUE);
        // initialize each prime of the factor base with its index.
        int baseSizeMax = (int) factorBaseSize((long)smoothBound*smoothBound);
        for (int i=1; i<baseSizeMax; i++){
            maxSmoothFactorIndex[primes[i]] = (int) i;
        }
        maxSmoothFactorIndex[1] = 1;
            int smoothCount = 0;
            for (long j = 2; j < smoothBound; j++) {
                double baseSize = factorBaseSize(j * j);
                SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
//                boolean isSmooth = smallFactoriser.factor((int) j, baseSize, factors);
                // since we store the high factors for small numbers use the maximal factor base
                int factorIndex = smallFactoriser.findSingleFactorIndex((int) j, baseSizeMax) + 1;
                if (factorIndex > 0) {
                    int jDivFactor = (int) (j / primes[factorIndex]);
//                    if (isSmooth != smoothNumber.get(jDivFactor)) {
//                        System.out.println();
//                    }
                    int maxFactorIndex = maxFactorIndex(jDivFactor);
                    maxFactorIndex = max(factorIndex, maxFactorIndex);
                    maxSmoothFactorIndex[(int) j] = maxFactorIndex;
//                    if ((isSmooth && maxFactorIndex >= baseSize) || (!isSmooth && maxFactorIndex < baseSize) /*|| jDivFactor == 22801*/) {
//                        maxFactorIndex = maxFactorIndex(jDivFactor);
//                    }
//                    maxSmoothFactorIndex[(int) j] = maxFactorIndex;
                    if (/*isSmooth*/ maxFactorIndex < baseSize) {
                        smoothNumber.set((int) j);
//                        if (factors.size() == 0)
//                            System.out.println(" number " + j + " is smooth but has no factor");
                    } else {
//                        System.out.println(j + " is not smooth");
                    }
                }
            }
//        }
//        catch (IOException e) {
//            e.printStackTrace();
//        }
        writeJavaObject("smoothFactors" + smoothBits + ".dat", maxSmoothFactorIndex);
        writeJavaObject(SMOOTH_FILE_NAME + smoothBits + ".dat", smoothNumber);
        long end = System.currentTimeMillis();
        System.out.println();
        System.out.println("Searched for smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");
        System.out.println("smooth numbers : " + smoothNumber.cardinality());
        System.out.println("ratio : " + ((double)smoothBound)/ smoothNumber.cardinality());

        // we map blocks of size sqrt(x) to a range O(sqrt(x) / factorBaseSize(x))
        //
        // TODO we do not need max factor!?
        // around the same size as smooth numbers BitSet
        System.out.println("initialization ready ");
    }

    private int maxFactorIndex(int n) {
        int m = n;
        int maxFactorIndex = -1;
        int currentFactorIndex = -1;
        while(m > 1 && (currentFactorIndex = maxSmoothFactorIndex[m]) < Integer.MAX_VALUE){
            maxFactorIndex = max(maxFactorIndex, currentFactorIndex);
            m = m/primes[currentFactorIndex];
        }
        // either we found all factors, or there must be a factor outside the factor base
        return max(maxFactorIndex, currentFactorIndex);
    }

    private void writeJavaObject(String objectName, Object javaObject) {
        try (FileOutputStream fos = new FileOutputStream(objectName)) {
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(javaObject);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Object readJavaObject(String objectName) {
        try (FileInputStream fis = new FileInputStream(objectName)) {
            ObjectInputStream ois = new ObjectInputStream(fis);
            return ois.readObject();
        } catch (IOException | ClassNotFoundException e) {
            if (! (e instanceof  FileNotFoundException))
            e.printStackTrace();
        }
        return null;
    }

        /**
     * sqrt(exp(sqrt(log(n) * log(log(n))))) / 2
     * @param n
     * @return
     */
    private double factorBaseSize(long n) {
        final double logN = log(n);
        return .7 * sqrt(exp(sqrt(logN * log(logN)))) +3;
    }

}
