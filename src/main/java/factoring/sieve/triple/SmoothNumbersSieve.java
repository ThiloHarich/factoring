package factoring.sieve.triple;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.SortedMultiset_BottomUp;
import factoring.math.Row;
import factoring.math.SquareFinder3;
import factoring.trial.TDiv31Barrett;

import java.io.*;
import java.math.BigInteger;
import java.util.*;

import static java.lang.Math.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Starting with the square t^2 , t = ceil(sqrt(n)) + i we
 * have t^2 - n < 2*i*sqrt(n)<br>
 * we can build<br>
 * (t-c)*(t+c) - n = t^2 - c^2 -n<br>
 * with c = sqrt(t^2 - n) + j < sqrt(2*i*sqrt(n))+j < sqrt(2i) * n^1/4 + j we have :<br>
 * (t-c)*(t+c) - n = t^2 - c^2 -n < t^2 - n - (sqrt(t^2 - n) +j)^2<br>
 * < t^2 - n - (t^2 - n) + 2j*sqrt(t^2 - n) + j^2<br>
 * < 2j* sqrt(t^2 - n) + j^2<br>
 * < j*sqrt(8i) *n^1/4 + j^2<br>
 * With subpolynomial i and j we will have numbers of size O(n^1/4+epsilon) compared to O(n^1/2 + epsilon) of the
 * regular Quadratic Sieve.
 * For a factor base which size depends on n we need an efficient way to<br>
 * a) determine smooth numbers (t-c)*(t+c) for a fixed t<br>
 * b) get the factorization of a number of size sqrt(8i) *n^1/4 + j^2 < 4*n^1/2<br>
 * In this algorithms we give efficient algorithms for this problems by lookup tables.
 * For b we store the biggest factor of all numbers below a limit around n^1/2.
 * By retrieving the biggest factor we can decide fast if the number is smooth over the factor base.
 * The a) we have a bit set of smooth numbers in ascending and descending order.
 * We then calculate the and of smooth numbers above t and the smooth numbers below t (in descending order).
 * For b we need to store the index of the prime. For ~ 54 Bits this fits in a byte.
 * We use an int array -> storing solutions for b takes 8 times the memory then the solutions for a), if both
 * numbers have the same size.
 * We might use a hash here to save memory.
 *
 * In the matrix step we not only have to guarantee that the right side is a square,
 * we also have to make sure that the left side is a square as well.
 * If we just make sure that the product of the left and the right side is a square we can make
 * both sides of the relation a square by multiplying with a prime p if both sides have an odd exponent.
 * If we multiply a relation l * p == r * p mod n with p we get squares on both sides of the relation:
 * l * p^2 == r * p^2 mod n
 * With this extension, the matrix consists out of |factor base| columns.
 *
 * We choose j such that the numbers on the right are below n^1/2, because
 * the numbers on the left (t+c) side will be around n^1/2.
 * By introducing an aditional factor on the left we can bring both sides below n^(3/8) but the algorithm
 * will be more complicated.
 *
 *
 */

// TODO it does not work when extending FactorAlgorithm
public class SmoothNumbersSieve extends FactorAlgorithm {

    private static final String SMOOTH_NUMBERS = "smoothNumbers";
    public static final String FACTORISATION = "factorisations";
    public static final String SMOOTH_EXPONENTS = "exponentsMod2";
    public static final String MAX_FACTOR = "maxFactor";
    public static final String MAX_SMOOTH_FACTOR = "maxSmoothFactor";
    public static final Long NO_SMOOTH_FACTOR = Long.MAX_VALUE;

    public static int [] primes = {-1, 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
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
            1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,
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
			3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733};

    private static boolean print = true;
    private static boolean showTiming = false;
    private static boolean check = false;
    // for a range we hash the smooth numbers.
    private static BitSet smoothNegative;
    private static BitSet smoothNumber;
    // we only store factorisations up to smoothBound / factorisationFraction;
    public static int factorisationBound;
    public static long [] factorisation;


    // TODO only use the hash to reduce memory.
    // instead of 2* n^1/2 bits we only need n^1/2 / sqrt(exp(sqrt(log(n) * log(log(n)))))
    public static short[] maxFactorIndex;
    public static short[] maxSmoothFactorHashed;
    public static int factorisationWords;



    //    long[] smoothEponents;
    BitSet[] exponentsMod2;
    long[] exponentsMod2AsLong;
    private int smoothBits = 24;
    private int smoothBound = 1 << smoothBits;
    final int hashMod = smoothBound / 4;
    private int hashMask = hashMod - 1;
    final int exponentsStored = smoothBound / 64;
    final int maxFactorAlLLimit = smoothBound / 16;


    TDiv31Barrett smallFactoriser = new TDiv31Barrett();
    private static int steps = 0;
    private int allTries;
    private int tries;
    private int[] primeCountPos;

    Set<Integer> relations;
//    long[] relationsHash;
    private int[] factorBase;
    private int exponentsMod2Words;

    //    @Override
    public String getName() {
        return "Smooth numbers ";
    }

//    @Override
    public BigInteger findSingleFactor(BigInteger N) {
        // TODO we need a static initializer here.
        initSmoothNumbers();

        long start = System.nanoTime();
        // we only need to factorize over the factor base.
        long n = N.longValue();

        final int factorBaseSize = (int) Math.ceil(factorBaseSize((long) n));
        relations = new HashSet<>();
        // make the hash table a power of 2 and hits rate < 1/8
//        relationsHash = new long[8 << PrimeMath.ceil(factorBaseSize)];
//        int relationsBits = 31 - Integer.numberOfLeadingZeros(relationsHash.length);

        int maxPrime = primes[factorBaseSize - 1];
        if (print){
            System.out.println("n : " + n);
            System.out.println("factor base size : " + factorBaseSize);
            System.out.println("max factor : " + maxPrime);
//            int length = 0;
//            for (int i=1; i< factorBaseSize; i++){
//                final int bits = 32 - Integer.numberOfLeadingZeros((int) Math.ceil(maxPrime / primes[i]));
//                length += bits;
//            }
//            System.out.println(length + " bit needed for the exponents");
        }
//        Integer[] factorBaseArr = IntStream.of(primes).limit(factorBaseSize).boxed().toArray(Integer[]::new);
        factorBase = new int[factorBaseSize];
        System.arraycopy(primes, 0, factorBase, 0, factorBaseSize);
        SquareFinder3 finder = new SquareFinder3(factorBase, n);
        allTries = 0;
        tries = 0;
        double sqrtNBits = Math.log(n)/Math.log(2)/2;
        // for big numbers this has to be 1 (or lower) otherwise we get an IndexArrayOutOfBounds
        // with this factor we can control the size of the numbers on the right side.
        // for small numbers the numbers on the right have to be higher.
        long t = (long) Math.ceil(Math.sqrt((double) (n)));
        // squares have no serach range, so we will directly return
        if (t*t == n)
            return BigInteger.valueOf(t);
        long tDiff = 0;

        // if we want to ensure that we can factozize every number with prob 1, we need to add the bitlength to the
        // factor base size
        int nBits = N.bitLength();
        while (relations.size() < 1 *factorBaseSize) {
            // TODO we only need the factors up to a bound. Which?
            // c is the minimal c for |(t^2 - c^2) - n|
            int factor = sieve(n, factorBaseSize, finder, t, tDiff);
            tDiff++;
            t++;
            if (factor > 0)
                return BigInteger.valueOf(factor);
        }
        long endSieve = System.nanoTime();
        if (print || showTiming)
        {
            System.out.println("time sieve  : " + (endSieve - start));
        }
        if (print)
            System.out.println("max multiplier : " + (tDiff + 1));

        if (false)
            return BigInteger.valueOf(1);
        long factor = solveMatrix(start, n, finder, endSieve, factorBaseSize);
        if (factor > 0)
        {
            return BigInteger.valueOf(factor);
        }
        // For most of the numbers the first solution of the matrix leads to a factor. Unfortunately some numbers like
        // 8299 we can find 14 = 21-7 colums with only '0' (ven exponents), but none of them leads to a factor
        // it might be that solutions of the matrix can not lead to a factor for some unkown reason.
        // we might get even more relations by extending the search. More negative numbers on the right of the equation.
        // at the moment the algorithm just fails for these numbers.
        // Only a verry small portion of numbers ~2% ? reach this special handling were we try to get as many relations
        // as we can.
        if (false) {
            while (relations.size() < 3 * factorBaseSize && tDiff < 150) {
                // TODO we only need the factors up to a bound. Which?
                // c is the minimal c for |(t^2 - c^2) - n|
                factor = sieve(n, factorBaseSize, finder, t, tDiff);
                tDiff++;
                t++;
                if (factor > 0)
                    return BigInteger.valueOf(factor);
            }
            factor = solveMatrix(start, n, finder, endSieve, factorBaseSize);
            if (factor > 0) {
                return BigInteger.valueOf(factor);
            }
        }
        return BigInteger.valueOf(factor);
    }


    private int sieve(long n, int factorBaseSize, SquareFinder3 finder, long t, long tDiff) {
        int lastRelationSize = 0;
        double c = sqrt(t * t - n);
        // if we have an upperbound u of smooth numbers:
        // (t^2 - (c+l)^2) - n > -u >= -n^1/2
        // (t^2 - (c^2+2cl+l^2)) - n > -u
        //  2cl+l^2 < u
        //  since we ensure l < c it holds : l^2 < cl -> 2cl+l^2 <= 3cl
        //  -> l < u/3
//        double sieveIntervalHalf1 = Math.sqrt(smoothBound/2);
        // TODO to have bigger positive numbers we might shift c to the left
        double sieveIntervalHalf2 = smoothBound / (3.0 * c);
//        if (sieveIntervalHalf1 > sieveIntervalHalf2)
//            System.out.println("need both sieve intervals");
        // sieve interval can not be greather then c
//        double sieveIntervalHalf = min(smoothFactor * min(sieveIntervalHalf1, sieveIntervalHalf2), c);
        double sieveIntervalHalf = min(sieveIntervalHalf2, c);
//            List<Integer> smoothPos = smoothPosPre((int) t, (int)Math.round(c), (int)sieveIntervalHalf);
        // TODO dist should be a multiple of 64 to improve performance of BitSet.get
        // TODO clean up up int variables
        // dist must be much greater then the distance between smooth numbers
        int dist = (int) sieveIntervalHalf + 1;
        int cInt = (int) Math.round(c);
        int tInt = (int) t;
        long tN = tInt * tInt - n;
        int begin = cInt - dist;
        BitSet smoothPosDistHash = getBitSet(dist, cInt, tInt, begin);
        tries++;
        steps++;
        long tryDist;
        int factor = -1;
        int exponentsSize = (int) ceil(factorBaseSize / 64.0);
        // TODO instead of BitSet we might directly iterate over the long array
        for (int smoothDist = smoothPosDistHash.nextSetBit(0); smoothDist >= 0 && factor < 0 && tDiff < 100;
             smoothDist = smoothPosDistHash.nextSetBit(smoothDist + 1)) {
            factor = getFactor(n, factorBaseSize, finder, tInt, tN, begin, factor, exponentsSize, smoothDist);
        }
        final int hits = relations.size() - lastRelationSize;
        if (print && hits > 0) {
            double hitRate = hits == 0 ? Integer.MAX_VALUE : (tries + 0.0) / hits;
            SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
            smallFactoriser.factor((int) t, factorBaseSize, factors);
            System.out.println("Mulitplier : " + (tDiff + 1));
            System.out.println("best dist : " + c);
            System.out.println("sieve interval : " + 2 * dist);
            System.out.println("num hits  " + hits);
            System.out.println("rel  hit rate 1 to " + hitRate);
            System.out.println("work done : " + (100 * relations.size()) / factorBaseSize + "%");
            System.out.println();
        }
        return factor;
    }

    private BitSet getBitSet(int dist, int cInt, int tInt, int begin) {
        BitSet smoothPosDistHash = smoothNumber.get(tInt + begin, tInt + (cInt + dist));
        BitSet smoothNegDistHash = smoothNegative.get(smoothBound - (tInt - begin), smoothBound - (tInt - cInt - dist));
        smoothPosDistHash.and(smoothNegDistHash);
        return smoothPosDistHash;
    }

    private int getFactor(long n, int factorBaseSize, SquareFinder3 finder, int tInt, long tN, int begin, int factor, int exponentsSize, int smoothDist) {
        long tryDist;
        tryDist = begin + smoothDist;
        tries++;
        steps++;
        int right = (int) (tN - tryDist * tryDist);
//              // this might be the most time consuming step since we have (pseudo) random memory access
        // if we know first smoothPos we might access (smoothPos + i)^2 = smoothPos^2 + i *(smoothPos + i) fast!?
        final int rightPos = abs(right);
        // when only checking if the number is smooth we might miss small right numbers
        // when checking/finding the biggest smooth prime, we can divide it out and lookup the smaller factorization
//        boolean isRightSmooth = aFactorIndex(rightPos) < factorBaseSize;
        boolean isRightSmooth2 = smoothNumber.get(abs(right));
        if (false && right < Math.sqrt(n)){
            // TODO we get ~ 15% more relations with this. is it worse the efford?
            final long[] rightExponent = exponentsMod2(abs(right), exponentsSize, false);
            final int maxFactorIndex = getMaxFactorIndex(rightExponent);
            if (rightExponent.length != 0 && maxFactorIndex <= factorBaseSize && maxFactorIndex >= 0 && !isRightSmooth2)
                System.out.println("found additional small right " + right + " < " + Math.sqrt(n) + " with max factor "
                        + primes[getMaxFactorIndex(rightExponent)] + " and index "  + maxFactorIndex +
                        " in factor base of size " + factorBaseSize + " exponents mod 2 : " + exponentsToString(rightExponent));
        }
//            if (isRightSmooth2 ^ isRightSmooth) {
//                if (isRightSmooth2)
//                    System.out.println("rigth " + right + " is smooth but max factor " + primes[maxFactorIndex] + " outside the factor base");
//                else
//                    System.out.println("rigth " + right + " is not smooth but max factor " + primes[maxFactorIndex] + " is inside the factor base");
//            }
//            boolean rightCanBeADuplicat = containsHash(right, relationsHash, hashMask);
//            final boolean contains = relations.contains(right);
        // we first check if the hash of the rigth number machtes the hash of a already added right number
        if (isRightSmooth2 &&  !relations.contains(right))
//                if (isRightSmooth && (!containsHash(right, relationsHash, hashMask) || !relations.contains(right)))
//                if (isRightSmooth2 &&
//                        maxSmoothFactorIndex[rightPos] < factorBaseSize &&
//                        (!containsHash(right, relationsHash, hashMask) || !relations.contains(right)))
        {
            final long[] rightExponent = exponentsMod2(abs(right), exponentsSize, true);
            if (rightExponent.length == 0 || getMaxFactorIndex(rightExponent) > factorBaseSize) {
                if (print)
                    System.out.println("max exponent " + getMaxFactorIndex(rightExponent) + " for right " + right + " outside factor base of size " + factorBaseSize);
            }
            else {
                if (right < 0)
                    rightExponent[0] |= 1;
                checkExponentsMod2(rightExponent, right);
                final long higher = tInt + tryDist;
                final long[] higherExponent = exponentsMod2((int) higher, exponentsSize, true);
                checkExponentsMod2(higherExponent, higher);
                if (higherExponent.length == 0 || getMaxFactorIndex(higherExponent) > factorBaseSize)
                    System.out.println("max exponent " + getMaxFactorIndex(higherExponent) + " for hiher " + higher + " outside factor base of size " + factorBaseSize);
                else {

//            boolean isRightSmooth = tries % 10 == 0;
//            if (maxFactorInFactorBase) {
                    final long lower = tInt - tryDist;
                    // TODO we have to
//                    int right = (int) (lower * higher - n);
                    factor = addRelation(relations, finder, (int) lower, (int) higher, higherExponent, right, rightExponent,factorBaseSize, exponentsSize, n);
//                relations.add(right);
//            setHash(right, relationsHash, hashMask);
                }
            }
        }
        return factor;
    }

    /**
     * retuns a factor dividing number.
     * We store the maximal factor dividing the number. Since we might have collitions we might return a max factor for
     * a different (smaller) number, which is not the maximal factor of this number.
     * @param number
     * @return
     */
    private int aFactorIndex(int number) {
        if (number < maxFactorAlLLimit)
            return maxFactorIndex[number];
        int index = (hash(number)+ hashMod) & hashMask;
        int factorIndex = maxSmoothFactorHashed[index];
        while (factorIndex != 0){
            final int prime = primes[factorIndex];
            if ((number / prime) * prime == number)
                return factorIndex;
            index++;
            index = (index) & hashMask;
            factorIndex = maxSmoothFactorHashed[index];
        };

        return Integer.MAX_VALUE;
        // TODO if we find a factor outside of the factorbase we should store the factor, tN, tryDist, tInt
        // to be able to reconstruct the number if there is a second for this factor and reuse
        // this would give us some relations for free, since we have the max factor by design
    }



    private int addRelation(Set<Integer> relations, SquareFinder3 finder, int lower, int higher, long[] higherExponent, int right, long[] rightExponent, int factorBaseSize, int exponentsSize, long n) {
        // ensure that this relation is not a duplicate, and all numbers fit in the matrix
        // right should not be dividable by k
        // TODO use a int [] based map
//        if (relations.contains(right))
//            return -1;
        // for small numbers we do not have to wait that the whole matrix is build
        if (right == 0 && lower != 1) {
            int max = max(higher, lower);
            if (print) System.out.println("right side is 0. Directly found factor " + max);
            return max;
        }
        long[] exponents = getExponents(lower, higherExponent, rightExponent, exponentsSize);
//        if (right < 0)
//            rightExponents[0] |= 1;


//        if (maxSmoothFactorIndex[Math.abs(higher)] >= factorBaseSize) {
//            if (print) System.out.println("higher number 't+c' has a factor " + primes[maxSmoothFactorIndex[Math.abs(higher)]] + " outside the factor base.");
//            return -1;
//        }
//        long[] exponentsLeft = new long[longsPerNumber];
//        long[] exponentsRight = new long[longsPerNumber];
//        for (int i = 0; i < longsPerNumber; i++) {
//            exponentsLeft[i] = smoothEponents[lower * longsPerNumber + i];
//            exponentsLeft[i] += smoothEponents[higher * longsPerNumber + i];
//            exponentsRight[i] = smoothEponents[right * longsPerNumber + i];
//        }


        // TODO we might only use relationsHash
        relations.add(right);
//        Multiset<Integer> leftFactors = HashMultiset.create();
//        Multiset<Integer> rightFactors = HashMultiset.create();
//        addFactors(lower, leftFactors);
//        addFactors(higher, leftFactors);
//        addFactors(right, rightFactors);
//        long prodLeft = leftFactors.stream().map(v -> (long) v).reduce(1l, (a, b) -> a*b);
//        long prodRight = rightFactors.stream().map(v -> (long) v).reduce(1l, (a,b) -> a*b);
//        if (prodLeft % n != (prodRight +n) %n)
//            System.out.println(prodLeft + " % " + n + " = " + prodRight + " but should be " +  prodLeft % n);
//        finder.addFactors(leftFactors, rightFactors);
//        long[] exponents = {leftExponents ^ rightExponents};
//        System.out.println("try at dist " + (higher - t));
        SquareFinder3.Relation rel = finder.addFactors(lower, higher, right, exponents, factorBaseSize);
        if (print) {
            final long[] rightExponents = exponentsMod2(abs(right), exponentsSize, true);
            if (right < 0)
                rightExponents[0] |= 1;
            final long[] lowerExponents = exponentsMod2(lower, exponentsSize, true);
            final long[] higherExponents = exponentsMod2(higher, exponentsSize, true);

             // can  happen if k is high or higher value is to high -> the generated numbers might be bigger then sqrt(n)
//            if (getMaxFactorIndex(leftExponents) >= factorBaseSize) {
//                System.out.println("left numbers " + lower + " * " + higher + "  have a factor " + primes[getMaxFactorIndex(leftExponents)] + " outside the factor base.");
//            }
            System.out.println("tries " + tries);
            System.out.println("found smooth lower  " + lower + " odd factors : " + exponentsToString(lowerExponents));
            System.out.println(" and  smooth higher " + higher + " odd factors : " + exponentsToString(higherExponents));
//            System.out.println(" odd factors left                   : " + exponentsToString(leftExponents) +
//                    " : " + SquareFinder3.printExponents(rel.getLeftFactrorisation(), factorBase));
            //                        SortedMultiset<BigInteger> factorsRight = smallFactoriser.factor(BigInteger.valueOf(Math.abs(right)));
            System.out.println(" and  smooth right  " + right + " odd factors : " + exponentsToString(rightExponents) );
//                    " : " + SquareFinder3.printExponents(rel.getRightFactrorisation(), factorBase));
            System.out.println();

            System.out.println("------ new relation ------------");
        }
        return -1;
    }

    private long[] getExponents(int lower, long[] higherExponent, long[] rightExponent, int exponentsSize) {
        long [] exponents = new long[exponentsSize];

        // TODO we might chnage order of array access
//        for (int i = 0; i < exponentsSize; i++) {
//            lowerExponents[i] = exponentsMod2AsLong[lower * exponentsMod2Words + i];
//            higherExponents[i] = exponentsMod2AsLong[higher * exponentsMod2Words + i];
        final long[] lowerExponent = exponentsMod2(lower, exponentsSize, true);
        checkExponentsMod2(lowerExponent, lower);
//            rightExponents[i] = rightExponent;
//            final long leftExponent = exponentsMod2AsLong[lower * exponentsMod2Words + i] ^ exponentsMod2AsLong[higher * exponentsMod2Words + i];
//            leftExponents[i] = leftExponent;
            // TODO if everything is running only this is needed
        for (int i = 0; i < exponentsSize; i++) {
            exponents[i] = rightExponent[i] ^ lowerExponent[i] ^ higherExponent[i];
        }
//        checkExponentsMod2(exponents, ((long)lower) * ((long)higher) * ((long)right));
        return exponents;
    }

    private void checkExponentsMod2(long[] exponentsMod2, long number) {
        if (!check)
            return;
        long n = number;
        for (int i = 1; i < 64 * exponentsMod2.length; i++) {
            int exp = 0;
            while ((number / primes[i]) * primes[i] == number){
                number /= primes[i];
                exp++;
            }
            if ((exponentsMod2[i >> 6] & (1l << i)) == (1l<< i)){
                 if (exp % 2 == 0)
                    System.out.println(" wrong exponent " + exponentsToString(exponentsMod2) + " for number : " + n);
            }
            else{
                if (exp % 2 != 0)
                    System.out.println(" wrong exponent " + exponentsToString(exponentsMod2) + " for number : " + n);
            }
        }

    }

    private long[] exponentsMod2(int number, int exponentsSize, boolean checkSmooth) {
        if (number < exponentsStored) {
            long[] exponents = new long[exponentsSize];
            if ( exponentsMod2AsLong[number*exponentsMod2Words] == NO_SMOOTH_FACTOR)
                return new long [] {};
            for (int i = 0; i < exponentsSize; i++) {
                exponents[i] = exponentsMod2AsLong[number*exponentsMod2Words + i];
            }
            return exponents;
        }
        // TODO create a method which directy returns number / maxFactor
        int maxFactorIndex = aFactorIndex(number);
        if (maxFactorIndex == Short.MAX_VALUE || maxFactorIndex >= exponentsSize * 64) {
            // TODO handle case when we did not find the maximal number before
//            if (checkSmooth)
//                System.out.println("number not smooth :" + number + " max factor index : " + maxFactorIndex);
            return new long [] {};
        }
        long[] exponents = exponentsMod2(number / primes[maxFactorIndex], exponentsSize, checkSmooth);
        if (exponents.length == 0)
            return exponents;
        exponents [maxFactorIndex >> 6] ^= (1l << maxFactorIndex);
        return exponents;
    }

    public static String exponentsToString(long [] exponents) {
        String s = "";
        for (int i = 0; i < 64 * exponents.length; i++) {
            if ((exponents[i >> 6] & (1l << i)) == (1l<< i)){
                s += primes[i] + " ";
            }
        }
        return s;
    }



//    private void addFactors(final int number, final Multiset<Integer> factors) {
//        if (number < 0)
//            factors.add(-1);
//        int x = Math.abs(number);
//        while (x != 1){
//            final int factor = primes[maxSmoothFactorIndex[x]];
//            factors.add(Integer.valueOf(factor));
//            // TODO replace division
//            x /= factor;
//        }
//    }

    private long solveMatrix(long start, long n, SquareFinder3 finder, long endSieve, int factorBaseSize) {
        List<Row> smoothMatrix = finder.matrix;
//        int [] colCounts = new int [factorBaseSize];

//        List<Column> reducedMatrix = finder.reduceMatrix(smoothMatrix, colCounts);

        long factor = finder.doGaussElimination(smoothMatrix);

        long endMatrix = System.nanoTime();
        if (print || showTiming)
            System.out.println("time matrix : " + (endMatrix - endSieve));
        if (print) {
            System.out.println("Number of relations considered overall : " + steps);
            System.out.println("Runtime calculated exponent : " + Math.log(steps) / Math.log(n));
        }

        return factor;
    }

    private void initSmoothNumbers() {
        final long start = System.currentTimeMillis();
        long n = 1L << smoothBits;
        // 23 bits is ~ 2 MB
        // 34 bits is ~ 4 GB -> in best case for 90 Bit numbers
//        int maxFactorIndex = (int) factorBaseSize(smoothBound);
        int baseSizeMax = (int) factorBaseSize((long)smoothBound*smoothBound) - 1;
        int maxFactor = primes[baseSizeMax];
//        int[] primeCountPos = getPrimeCountPos(baseSizeMax, maxFactor);
//        factorisationWords = (int) ceil(baseSizeMax / 8.0);
//        exponentsMod2Words = (int) ceil(baseSizeMax / 64.0);
        exponentsMod2Words = (int) ceil(factorBaseSize((long)smoothBound/64l *smoothBound) / 64);
        factorisationBound = smoothBound / 32;

        if (smoothNumber != null /*&& maxSmoothFactorIndex != null*/)
            return;
        final boolean load = false;
        if (load) {
            // TODO just serialize the BitSet
    //        File smoothNumbersFile = new File(SMOOTH_FILE_NAME + smoothBits + ".txt");
            smoothNumber = (BitSet) readJavaObject(SMOOTH_NUMBERS + smoothBits + ".dat");
            maxSmoothFactorHashed = (short[]) readJavaObject(MAX_SMOOTH_FACTOR + smoothBits + ".dat");
            maxFactorIndex = (short[]) readJavaObject(MAX_FACTOR + smoothBits + ".dat");
            exponentsMod2AsLong = (long[]) readJavaObject(SMOOTH_EXPONENTS  + smoothBits + ".dat");
//            factorisation = (long[]) readJavaObject(FACTORISATION + smoothBits + ".dat");

//            exponentsMod2 = (BitSet[]) readJavaObject(SMOOTH_EXPONENTS  + smoothBits + ".dat");
//            smoothEponents = (long[]) readJavaObject(SMOOTH_EXPONENTS  + smoothBits + ".dat");
            if (smoothNumber != null /*&& maxSmoothFactorIndex != null*/ && exponentsMod2AsLong != null){
                long end = System.currentTimeMillis();
                System.out.println("read in file with smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");

                // we still have to determine smoothNegative
                smoothNegative = new BitSet(smoothBound);
                smoothNumber.stream().forEach(i -> smoothNegative.set(smoothBound - i));
                return;
            }
        }

        // no file for the smoothBits found -> calculate the numbers
        smoothNumber = new BitSet(smoothBound);
//        factorisationIndex = new int[smoothBound];
        maxFactorIndex = new short[maxFactorAlLLimit];
        maxSmoothFactorHashed = new short [hashMod];
//        exponentsMod2 = new BitSet[smoothBound];
        exponentsMod2AsLong = new long[exponentsStored * exponentsMod2Words];
//        factorisation = new long[factorisationBound * factorisationWords];
        smoothNumber.set(1);
//        Arrays.fill(maxSmoothFactorIndex, Integer.MAX_VALUE);
        // initialize each prime of the factor base with its index.
        for (int i=1; i<baseSizeMax; i++){
            maxFactorIndex[primes[i]] = (byte) i;
        }

        maxFactorIndex[0] = Short.MIN_VALUE;
        maxFactorIndex[1] = 1;

        // We need a fast initialization of the matrix step
        // for solving the matrix we just need the factors mod 2 -> BitSet does it.
        // For getting the solution we need to mulitply all factors of the combinded solution.
        // we just multiply the numbers together and apply Shanks tonelli.
        // for each factor we reserve a section which can represent a value of Size max Factor.
        // -> we need log_2(maxFactor / factor) bits to repesent
        // for the i-th prime we start at loc Sum_j<i log_2(maxFactor/i*logi)
        // = Sum_j<i (log_2(maxFactor) - log(i*log(i)))
        // = i* log_2(maxFactor) - Sum_j<i log(i*log(i))
        // ~ i * (log_2(maxFactor)-1)
        // for maxFactor = 512 -> log_2(maxFactor) = 9
        // i=1 , p=2 -> 8 bits
        // i=2 , p=3 -> 8 bits
        // i=3 , p=5 -> 7 bits
        // i=4 , p=7 -> 7 bits
        // i=5 , p=11 -> 6 bits
        // i=6 , p=17 -> 5 bits
        // i=? , p=maxFactor/2 -> 1 bits
//        smoothEponents = new long[longsPerNumber * smoothBound];

        // TODO sieve here for performance
//        int smoothCount = 0;
        int factIndex = 0;
        double maxFactorSmoothCount = 0;
        int collition = 0;
//        for (long j = 5087; j < smoothBound; j++) {
        for (long j = 2; j < smoothBound; j++) {
//            int baseSize = (int) Math.ceil(factorBaseSize(j * j));
//            int factorizationWords = (int) Math.ceil((baseSize+1) * .125);
            SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
            // since we store the high factors for small numbers use the maximal factor base
//            int factorIndex = smallFactoriser.findSingleFactorIndex((int) j, baseSizeMax) + 1;
            // we take 8 bits per count, a word has 64 bits
            long [] factorCounts = new long[factorisationWords];
            if (j == 25348)
                System.out.println();
            // we store
            // the max fastor Index of smooth numbers above smoothBound / 8 as byte
            // the max fastor Index of all numbers below smoothBound / 8 as byte
            // the exponents mod 2 of all numbers below smoothBound / 64*i as long []
            TDiv31Barrett.ExponentsMaxFactor result = smallFactoriser.fillFactorCounts((int) j, baseSizeMax);
            // if we define it like this, there can only be smoothBound bits used
            if (result.maxFactorIndex < Short.MIN_VALUE || result.maxFactorIndex > Short.MAX_VALUE)
                System.out.println("XXXXXXXXX");
            if (j < maxFactorAlLLimit) {
                maxFactorIndex[(int) j] = (short) result.maxFactorIndex;
            }
            else if (result.isSmooth){
                maxFactorSmoothCount++;
                int index = (hash((int) j)+ hashMod) & hashMask;
                while (maxSmoothFactorHashed[index] != (short) result.maxFactorIndex && maxSmoothFactorHashed[index] != 0){
                    index++;
                    index = (index) & hashMask;
                    collition++;
                }
                maxSmoothFactorHashed[index] = (short) result.maxFactorIndex;
            }
//            else{
//                int maxFactorIndex = getMaxFactorIndex(factorCounts)
//                final int prime = primes[maxFactorIndex];
//                int jDivFactor = (int) (j / prime);
//                int maxfactorIndexSub = maxSmoothFactorIndex[jDivFactor];
//                maxFactorIndex = max(maxfactorIndexSub, maxFactorIndex);
//                maxSmoothFactorIndex[(int)j] = maxFactorIndex;
//                // this is for the maxtrix only
//                final long[] longs = exponentsMod2.toLongArray();
//                for (int i = 0; i < longs.length; i++) {
            // if we define it like this there can only be smoothBound bits used
            if (j < exponentsStored){
                if (!result.exponents.isEmpty() || !result.isSmooth) {
                    int length = result.isSmooth ? result.exponents.toLongArray().length : exponentsMod2Words;
                    for (int i = 0; i < Math.min(exponentsMod2Words, length); i++) {
                        exponentsMod2AsLong[(int) (j * exponentsMod2Words) + i] = result.isSmooth ? result.exponents.toLongArray()[i]
                                : NO_SMOOTH_FACTOR;
                    }
                }
            }
//                maxSmoothFactorIndex[(int) j] = maxFactorIndex;
//                if (maxFactorIndex < Integer.MAX_VALUE) {
//                    if (maxFactorIndex > baseSizeMax)
//                        System.out.println();
//                    int number = (int) j;
//                    do {
//                        addPrimeIndex(primeCountPos[maxFactorIndex], longsPerNumber, number);
//                        number = number / primes[maxFactorIndex];
//                    } while (number > 1);
//                }

                // we get an bit representation of all smooth factors (but no info on how often it divides the number).
                // This a compromise between low memory (no count),
//                if (exponentsMod2 == null){
//                     If we have nt found a factor we set all factors to true -> factors outside considered factor base are set
//                    exponentsAsLong[(int)j] = Long.MAX_VALUE;
//                }
//                else{
                int baseSize = (int) factorBaseSize((long)j*j);
                // since smooth is for the left factors only, which are of nearly the same size n^1/2 (at least for
                // bigger numbers) we assume they have max Factors of the same size within the factor base
                if (result.maxFactorIndex <= baseSize) {
                    smoothNumber.set((int) j);
                }
//                smoothCount++;
//                factorisationIndex[(int)j] = factIndex;
                // we store the full factorization of a number. But we use it only for numbers divided by the biggest
                // prime factor only -> we just have to store them for a smaller portion of the numbers
//                if ((j+1) < factorisationBound) {
//                    for (int i = 0; i < factorisationWords; i++) {
//                        final int index = (int) (factorisationWords * j + i);
//                        factorisation[index] = factorCounts[i];
//                    }
//                factorCounts[factIndex] = factIndex;
//                }
//                }
//                int maxFactorFound = getMaxFactorIndex(exponentsAsLong[(int) j]);
//                if (/*isSmooth*/ maxFactorIndex < baseSize) {
//                    System.out.println(exponentsToString(exponentsMod2AsLong[(int)j]));
//                    boolean isSmooth2 = maxFactorFound < baseSize;
//                    assertTrue(isSmooth2);
//                    smoothNumber.set((int) j);
//                }else{
//                    assertTrue(maxFactor >= baseSize);
//                }
                // BitSet has really bad performance in serialization it is slow and memory wasting -> we use the internal long representation
//                exponentsMod2[(int) j] = factorIndexMod2;
//            }
        }
        double maxFactorSmoothRate = maxFactorSmoothCount / (smoothBound - maxFactorAlLLimit);

//        writeJavaObject(FACTORISATION + smoothBits + ".dat", factorisation);
        writeJavaObject(SMOOTH_NUMBERS + smoothBits + ".dat", smoothNumber);
        writeJavaObject(MAX_SMOOTH_FACTOR + smoothBits + ".dat", maxSmoothFactorHashed);
        writeJavaObject(MAX_FACTOR + smoothBits + ".dat", maxFactorIndex);
        writeJavaObject(SMOOTH_EXPONENTS + smoothBits + ".dat", exponentsMod2AsLong);
        long end = System.currentTimeMillis();
        System.out.println();
        System.out.println("Searched for smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");
        System.out.println("smooth numbers : " + smoothNumber.cardinality());
        System.out.println("ratio : " + ((double)smoothBound)/ smoothNumber.cardinality());
        System.out.println("initialization ready ");

        // we still have to determine smoothNegative
        smoothNegative = new BitSet(smoothBound);
        smoothNumber.stream().forEach(i -> smoothNegative.set(smoothBound - i));
    }
    int hash(int x) {
        x = ((x >>> 16) ^ x) * 0x45d9f3b;
        x = ((x >>> 16) ^ x) * 0x45d9f3b;
        x = (x >>> 16) ^ x;
        return x;
    }

    private int getMaxFactorIndex(long exponents) {
        return 63 - (Long.numberOfLeadingZeros(exponents));
    }
    private int getMaxFactorIndex(long[] factorCounts) {
        int factorIndex = -1;
        int wordIndex = factorCounts.length;
        while(factorIndex < 0 && wordIndex >= 1) {
            wordIndex--;
            long word = factorCounts[wordIndex];
            factorIndex = 63 - (Long.numberOfLeadingZeros(word));
        }
//        return (factorIndex >> 3) + (wordIndex << 3);
        return factorIndex + (wordIndex << 6);
    }

//    static boolean containsHash (int index, long[] hashSet, int hashMask){
//        int hashCode = index & hashMask;
//        final int wordIndex = hashCode >> 3;
////        final int indexInWord = (hashCode - (wordIndex << 3)) << 3;
//        return (hashSet[wordIndex] & 1l << hashCode) != 0;
//    }

//    static void setHash (int index, long[] hashSet, int hashMask){
//        // mod hashSet.length * 64 to fir in hashSet
//        int hashCode = index & hashMask;
//        final int wordIndex = hashCode >> 3;
////        final int indexInWord = (hashCode - (wordIndex << 3)) << 3;
//        hashSet[wordIndex] |= 1l << hashCode;
//    }


    private int[] getPrimeCountPos(int baseSizeMax, int maxFactor) {
        primeCountPos = new int[baseSizeMax +1];
        primeCountPos[0]=0;
        int length = 0;
        for (int i = 1; i<= baseSizeMax; i++){
            final int bits = 32 - Integer.numberOfLeadingZeros((int) Math.ceil(maxFactor / primes[i]));
            length += bits;
            primeCountPos[i] = length;
        }
        return primeCountPos;
    }

//    private void addPrimeIndex(int primeCountPos, int longsPerNumber, int number) {
//        final int primeBlock = primeCountPos / 64;
//        final int primeIndex = primeCountPos - (64 * primeBlock);
//        final int numberIndex = (int) (number * longsPerNumber) + primeBlock;
//        smoothEponents[numberIndex] += 1 << primeIndex;
//    }
//
//    private int maxFactorIndex(int n) {
//        int m = n;
//        int maxFactorIndex = -1;
//        int currentFactorIndex = -1;
//        while(m > 1 && (currentFactorIndex = maxSmoothFactorIndex[m]) < Integer.MAX_VALUE){
//            maxFactorIndex = max(maxFactorIndex, currentFactorIndex);
//            m = m/primes[currentFactorIndex];
//        }
//        // either we found all factors, or there must be a factor outside the factor base
//        return max(maxFactorIndex, currentFactorIndex);
//    }

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
        return .8 * sqrt(exp(sqrt(logN * log(logN))));
    }

}
