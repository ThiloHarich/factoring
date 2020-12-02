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
import factoring.math.SquareFinder2;
import factoring.trial.TDiv31Barrett;

import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static org.junit.Assert.assertEquals;

/**
 * Starting with the square t^2 , t = ceil(sqrt(n)) + i we
 * have t^2 - n < 2*i*sqrt(n)
 * can build
 * (t-c)*(t+c) - n = t^2 - c^2 -n
 * with c = sqrt(t^2 - n) + j < sqrt(2*i*sqrt(n))+j < sqrt(2i) * n^1/4 + j we have :
 * (t-c)*(t+c) - n = t^2 - c^2 -n < t^2 - n - (sqrt(t^2 - n) +j)^2
 * < t^2 - n - (t^2 - n) + 2j*sqrt(t^2 - n) + j^2
 * < 2j* sqrt(t^2 - n) + j^2
 * < j*sqrt(8i) *n^1/4 + j^2
 * With subpolynomial i and j we will have numbers of size O(n^1/4+epsilon) compared to O(n^1/2 + epsilon) of the
 * regular Quadratic Sieve.
 * For a factor base which size depends on n we need an efficient way to
 * a) determine smooth numbers (t-c)*(t+c) for a fixed t
 * b) get the factorization of a number of size sqrt(8i) *n^1/4 + j^2 < 4*n^1/2
 * In this algorithms we give efficient algorithms for this problems by lookup tables.
 * For b we store the biggest factor of all numbers below a limit around n^1/2.
 * By retrieving the biggest factor we can decide fast if the number is smooth over the factor base.
 * The a) we have a bit set of smooth numbers in ascending and descending order.
 * We then calculate the and of smooth numbers above t and the smooth numbers below t (in descending order).
 * In the matrix step we not only have to guarantee that the right side is a square,
 * we also have to make sure that the left side is a square as well.
 * If we just make sure that the product of the left and the right side is a square we can make
 * both sides of the relation a square by multiplying with a prime p if both sides have an odd exponent.
 * If we multiply a relation l * p == r * p mod n with p we get squares on both sides of the relation:
 * l * p^2 == r * p^2 mod n
 * With this extension, the matrix consists out of |factor base| columns.
 */

public class SmoothNumbersSieve extends FactorAlgorithm {

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
            1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063};

    // for a range we hash the smooth numbers.
    private BitSet smoothNegative;
    private BitSet smoothNumber;

    // TODO only use the hash to reduce memory.
    // instead of 2* n^1/2 bits we only need n^1/2 / sqrt(exp(sqrt(log(n) * log(log(n)))))
    int[] maxSmoothFactorIndex;
    private int smoothBits;
    private int smoothBound;


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

        Set<Integer> relations = new HashSet<>();
        final int factorBaseSize = (int) factorBaseSize((long) n);
        System.out.println("factor base size : " + factorBaseSize);
        int maxPrime = primes[factorBaseSize -1];
        System.out.println("max factor : " + maxPrime);
        int lastRelationSize = 0;
        Integer[] factorBaseArr = IntStream.of(primes).limit(factorBaseSize).boxed().toArray(Integer[]::new);
        SquareFinder2 finder = new SquareFinder2(Arrays.asList(factorBaseArr), n);
        int allTries = 0;
        int tries = 0;
        long t = (long) Math.ceil(Math.sqrt((double) (n)));
        long tDiff = 0;
        while (relations.size() <factorBaseSize + 3) {
             tDiff++;
             t++;
            // TODO we only need the factors up to a bound. Which?
            // c is the minimal c for |(t^2 - c^2) - n|
            double c = sqrt(t * t - n);
            // (t^2 - (c+l)^2) - n > -n^1/2
            // (t^2 - (c^2+2cl+l^2)) - n > -n^1/2
            //  2cl+l^2 < n^1/2
            //  l < n^1/4 and
            //  l < n^1/2 / 2c
            double sieveIntervalHalf1 = Math.pow(n, .25);
            double sieveIntervalHalf2 = Math.sqrt(n) / (2*c);
            double smoothFactor = 2;
//            double smoothFactor = 20;
            double sieveIntervalHalf = min(smoothFactor * min(sieveIntervalHalf1, sieveIntervalHalf2), c);
            List<Integer> smoothPos = smoothPosPre((int) t, (int)Math.round(c), (int)sieveIntervalHalf);
            long tN = t*t - n;
            for (int tryDist : smoothPos) {
                tries++;
                steps++;
                // TODO per calculate t^2 - n and only calc tryDist * tryDist
                int right2 = (int) (tN - tryDist*tryDist);
//                assertEquals(right, right2);
                // this might be the most time consuming step since we have (pseudo) random memory access
                // if we know first smoothPos we might access (smoothPos + i)^2 = smoothPos^2 + i *(smoothPos + i) fast!?
                boolean isRightSmooth2 = maxSmoothFactorIndex[Math.abs(right2)] < factorBaseSize;
                if (isRightSmooth2) {
                    final long lower = t - tryDist;
                    final long higher = t + tryDist;
                    int right = (int) (lower * higher - n);
                    addRelation(relations, finder, (int) lower, (int) higher, right, factorBaseSize, n);
                }
            }
            final int hits = relations.size() - lastRelationSize;
            int hitRate = hits == 0 ? Integer.MAX_VALUE : (tries / hits);
            SortedMultiset<BigInteger> factors = new SortedMultiset_BottomUp<>();
            smallFactoriser.factor((int) t, factorBaseSize, factors);
            System.out.println("t diff : " + tDiff + " : " + factors);
            System.out.println("best dist : " + c);
            System.out.println("sieve interval : " + 2 * sieveIntervalHalf);
            System.out.println("rel  " + hits + " hit rate 1 : " + hitRate);
            System.out.println("work done : " + (100 *relations.size()) / (factorBaseSize) + "%");
            System.out.println();
        }
        System.out.println("relations considered : " + steps);
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


    private void addRelation(Set<Integer> relations, SquareFinder2 finder, int lower, int higher, int right, int factorBaseSize, long n) {
        // ensure that this relation is not a duplicate, and all numbers fit in the matrix
        // right should not be dividable by k
        if (relations.contains(right))
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
        addFactors(lower, leftFactors);
        addFactors(higher, leftFactors);
        addFactors(right, rightFactors);
        int prodLeft = leftFactors.stream().reduce(1, (a,b) -> a*b);
        int prodRight = rightFactors.stream().reduce(1, (a,b) -> a*b);
        if (prodLeft % n != (prodRight +n) %n)
            System.out.println(prodLeft + " % " + n + " = " + prodRight + " but should be " +  prodLeft % n);
        finder.addFactors(leftFactors, rightFactors);
//        System.out.println("try at dist " + (higher - t));
//        System.out.println("tries " + tries);
        System.out.println("found smooth lower " + lower );
        System.out.println(" and  smooth higher" + higher);
        System.out.println(" factors           " + leftFactors);
        //                        SortedMultiset<BigInteger> factorsRight = smallFactoriser.factor(BigInteger.valueOf(Math.abs(right)));
        System.out.println(" and  smooth right " + right + " = " + rightFactors);
        System.out.println();
        System.out.println("------ new relation ------------");
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
     * Gives back the c+s  with t + c + s is smooth and t - c - s, -dist < s < dist is smooth.
     * We take the Binary representation of smooth numbers pos + i and
     * the Binary representation of smooth numbers pos - i.
     * Do an and, and iterate over the results. Since getting the next set bit is fast
     * when the dist is less then 64 (long bit length).
     * @param t
     * @param c
     * @param dist
     * @return
     */
    public List<Integer> smoothPosPre(int t, int c, int dist){
//        final double smoothDist = factorBaseSize(pos) / 8;
        // TODO dist should be a multiple of 64 to improve performance of BitSet.get
        BitSet smoothPosDistHash = smoothNumber.get(t + c - dist, t + (c + dist));
        BitSet smoothNegDistHash = smoothNegative.get(smoothBound-(t-c+dist), smoothBound - (t- c-dist));
        smoothPosDistHash.and(smoothNegDistHash);
        List<Integer> bothSmooth = new ArrayList<>();
//        while ((tryPos = smoothLowerHash.nextSetBit((int) pos)) <= dist) {
        // TODO directly use the iterator; do not use the List
        for (int smoothDist = smoothPosDistHash.nextSetBit(0); smoothDist >= 0; smoothDist = smoothPosDistHash.nextSetBit(smoothDist+1)) {
            int higherNumber = t + c - dist + smoothDist;
                // since we do not use a hash at the moment this is always true.
                if (smoothNumber.get(higherNumber) && smoothNumber.get(2 * t - higherNumber)) {
                    bothSmooth.add(higherNumber-t);
                }
                // TODO remove
                if (smoothDist == Integer.MAX_VALUE) {
                    break; // or (i+1) would overflow
                }
        }
        return bothSmooth;
    }


    private void initSmoothNumbers() {
        long start = System.currentTimeMillis();
        smoothBits = 27;
        long n = 1l << smoothBits;
        // 23 bits is ~ 2 MB
        // 34 bits is ~ 4 GB -> in best case for 90 Bit numbers
        smoothBound = 1 << smoothBits;

        boolean load = true;
        if (load) {
            // TODO just serialize the BitSet
    //        File smoothNumbersFile = new File(SMOOTH_FILE_NAME + smoothBits + ".txt");
            smoothNumber = (BitSet) readJavaObject(SMOOTH_FILE_NAME + smoothBits + ".dat");
            maxSmoothFactorIndex = (int[]) readJavaObject("smoothFactors" + smoothBits + ".dat");

            if (smoothNumber != null){
                long end = System.currentTimeMillis();
                System.out.println("read in file with smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");

                smoothNegative = new BitSet(smoothBound);
                smoothNumber.stream().forEach(i -> smoothNegative.set(smoothBound - i));
                return;
            }
        }

        smoothNumber = new BitSet(smoothBound);
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
                // since we store the high factors for small numbers use the maximal factor base
                int factorIndex = smallFactoriser.findSingleFactorIndex((int) j, baseSizeMax) + 1;
                if (factorIndex > 0) {
                    int jDivFactor = (int) (j / primes[factorIndex]);
                    int maxFactorIndex = maxFactorIndex(jDivFactor);
                    maxFactorIndex = max(factorIndex, maxFactorIndex);
                    maxSmoothFactorIndex[(int) j] = maxFactorIndex;
                    if (/*isSmooth*/ maxFactorIndex < baseSize) {
                        smoothNumber.set((int) j);
                    }
                }
            }
        writeJavaObject("smoothFactors" + smoothBits + ".dat", maxSmoothFactorIndex);
        writeJavaObject(SMOOTH_FILE_NAME + smoothBits + ".dat", smoothNumber);
        long end = System.currentTimeMillis();
        System.out.println();
        System.out.println("Searched for smooth numbers below " + smoothBound + " in : " + (0.0 + end - start) / 1000 + " seconds.");
        System.out.println("smooth numbers : " + smoothNumber.cardinality());
        System.out.println("ratio : " + ((double)smoothBound)/ smoothNumber.cardinality());
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
        return .7 * sqrt(exp(sqrt(logN * log(logN)))) + 5;
    }

}
