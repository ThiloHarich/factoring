package factoring.fermat.residue;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.math.SquaresMod;

import java.util.Collection;

import static factoring.math.PrimeMath.mod;

/**
 * Factorizes the number n by investigate the fermat equation for some residues.
 * This is x^2 - n = y^2 mod p. Due to the fact (-i)^2 = i^2, this equation
 * has at most (p+1)/2 solutions. If p is a prime, there are only (p+1)/2 values x^2 mod p.
 * So there can not be more solutions for a prime. We are putting together solutions for the small
 * primes. For p = 360 = 2^3 * 3^2 * 5 there are at most 24 (18.75 in average) solutions if n is
 * not dividable by 2,3 or 5. For each prime >5 multiplied to the Mod p the number of solutions is halved.
 * Unfortunately this only results in a speedup of O(exp( ln(n')/ld(ld(n')))) which is sub-polynomial.
 * To get from one solution for a residue which is the product of the first primes to the solution
 * it uses
 * s(360) = 360/24
 * s(360 * 7) = s(p)*2 , where q is a prime.
 * s(360 * first k primes greater then 7) = s(p)*2^k
 *
 * the product of the first k primes greather then p:
 * prod_i<k (i*log(i)) = 2^ ( sum_i<k ld(i*ld(i)))  , ld 0 logatithmus dualis
 * = 2^ ( sum_i<k ld(i)+ ld(log(i))))
 * = 2^ ( k *(ld(k)+ ld(ld(k))))
 * = 2^ ( k * ld(k)) * 2^( k * ld(ld(k))
 * = 2^ ( k * ld(k)) * ld(k)^k)
 * if  n' = 2^ ( k * ld(k)) * ld(k)^k) -> speedup 2^k
 * n' / ld(k)^k)  = 2^ ( k * ld(k))
 * ld (n') - ld(ld(k)^k) = k * ld(k)
 * ld (n')/ ld(k) - ld(ld(k)^k)/ ld(k) = k
 * -> k ~ ln(n')/ld(ld(n'))
 * 2^k = 2^( ln(n')/ld(ld(n')))
 * 2^k = 2^( ln(n')* (1/ld(ld(n')))
 *
 * -> there is only a subploimial speedup.
 *
 *
 * primorial n# is the product of all primes lower n
 *
 * Created by Thilo Harich on 29.11.2017.
 */
public class FermatResiduesRec extends FermatFact {

    private static final int MAX_FACTOR = 15;
    long minFactor = 3;
    protected int[] residueClasses;
    protected int [] residuesArray = {2,3,5,7,11,13,17,19}; // 25, 27, 49

    protected boolean[][] squares;
    protected int[][] xArray;
//    protected int[] xMask;
int[] modOldInvert;
    protected int[] modProd;

    protected int lastLevel = 0;
    private long xEnd;

    public FermatResiduesRec(int minFactor) {
        this.minFactor = minFactor;
    }

    public FermatResiduesRec() {
    }

    public void init (long n) {
        xEnd = (n / minFactor + minFactor) / 2;
        int [] residueClasses = calculateResidueClasses(xEnd);
        init(n, residueClasses);
    }

    public void initSquares(int level, int mod) {
        // we might also use (-i)^2 = i^2, and stop the loop in the middle, but does not help much
        squares[level] = new boolean[mod];
        int square = 0;
        for (int i = 1; i < mod+1; i += 2) {
            squares[level][square]= true;
            // we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
            // we might move out the last call
            square += i;
            square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
        }
    }
    public void initX(int level, int mod, long n)
    {
        xArray[level] = new int[mod + 1];
        int nMod = PrimeMath.mod(n, mod);
        int resIndex= 0;
        int squareMinN = nMod == 0 ? 0 : mod-nMod;

        for (int i = 1; i <2*mod; i += 2) {
            if (squares[level][squareMinN]) {
                xArray[level][resIndex++] = i / 2;

//                xMask[level] |= (1L <<  i / 2);
            }
            // we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
            // we might move out the last call
            squareMinN += i;
            squareMinN -= squareMinN >= mod ? (squareMinN >= 2*mod ? 2*mod : mod) : 0;
        }
        xArray[level][resIndex] = -1;
    }

//    public void xMask(int n)
//    {
//        xMask = 0;
//        int resIndex= 0;
//        int nMod = PrimeMath.mod(n, mod);
//        int squarePlusN = nMod;
//        if (((1L << squarePlusN) & squares) != 0) {
//            xMask |= 1L;
//        }
//        squarePlusN += 1;
//        squarePlusN -= squarePlusN >= mod ?  mod : 0;
////			for (int j = 1; j < 2*mod+1; j += 2) {
//        for (int j = 3; j < mod+1; j += 2) {
//            if (((1L << squarePlusN) & squares) != 0) {
//                xMask |= (1L << j/2);
//                xMask |= (1L << (mod - j/2));
//            }
//            squarePlusN += j;
//            squarePlusN -= squarePlusN >= mod ? mod : 0;
//        }
//        // since all primes different from 2 (which call xArray) are odd this will never be executed
//        if ((mod & 1) == 0 && ((1L << squarePlusN) & squares) != 0) {
//            xMask |= (1L << mod/2);
//        }
//    }

    public void init( long n, int[] residueClasses) {
        int prod = 1;
        modOldInvert = new int [residueClasses.length];
        modProd = new int [residueClasses.length];
        squares = new boolean [lastLevel][];
        xArray = new int [lastLevel][];
//        xMask = new int [lastLevel];
        for (int level=0; level < residueClasses.length-1 && residueClasses[level] > 0; level++) {
            int mod = residueClasses[level];
            initSquares(level, mod);
            initX(level, mod, n);
            if (mod > 0) {
                modProd[level] = prod;
                if (prod > 1) {
                    modOldInvert[level] = (int) PrimeMath.invert(prod, mod);
                }
                prod *= mod;
            }
        }
    }

    /**
     * Determine the residue classes for a given search interval, such that sieving with the
     * solutions has maximal speed.
     * This is we have a minimal number of solutions for the product of the generated
     * residueClasses.
     * if all P_i^(e_i) are equal it looks like a good solution.
     *
     * We get good results for products of the following numbers:
     * - 2 and small powers -> 8 (maximal : 4, avg : 2,3)
     * - 3 and small powers -> 9 (maximal : 2, avg : 2)
     * - any prime q (maximal (q+1)/2, avg )
     *
     *
     * @param searchInterval
     * @return
     */
    public int [] calculateResidueClasses(long searchInterval) {
        residueClasses = new int [residuesArray.length+1];
        int modProd = 1;
        int i = 0;
        for (; i < residuesArray.length && modProd * 4 *residuesArray[i] <= searchInterval && i<residuesArray.length; i++)
        {
            int prime = residuesArray[i];
            int modPerPrime = prime;
            modProd *= prime;
            while(modPerPrime * prime <= MAX_FACTOR)
            {
                modPerPrime *= prime;
                modProd *= prime;
            }
            residueClasses [i] = modPerPrime;
            lastLevel++;
        }
        // be sure we have at least one thing to search
        if (lastLevel == 0){
            lastLevel++;
            residueClasses [0] = 4;
        }

        return residueClasses;
    }

    public long findFactors(long n, Collection<Long> factors, int x, int level) {
        if (residueClasses[level] == 0) {
            // adjust sqrtN to be a multiple of modOld
            long xL = x;
//			if (modulus * mod  >  xRange/SEARCH_COUNT) {
            while (xL <= xEnd) {
                long right = xL * xL - n;
                if (SquaresMod.isSquare(right)) {
                    long y = PrimeMath.sqrt(right);
                    long factorHigh = xL + y;
                    long factorLow = xL - y;
                    factors.add(factorHigh);
                    return factorLow;
                }
                xL+= modProd[level-1];
            }
            return n;
        }
        long factor = findFactorsByMerge(n, factors, x, level);
        if (factor != n)
            return factor;
        return n;
    }

    /**
     * Here the solutions x modulo modProd , ... will be merged with the solutions
     * Modulo residueClasses[level] = modNew.
     * This is done by checking if the value x + i* modOld = xSolNew modulo modNew .
     * i* modProd = xSolNew - x modulo modNew
     * i = (xSolNew - x) * modProd^-1 modulo modNew
     *
     * This only works if the number n and the residueClasses have no common divisor
     * 0,4   mod 8
     * 0,3,6 mod 9
     *
     */
    public long findFactorsByMerge(long n, Collection<Long> factors, int x, int level) {
        for (int k = 0; xArray[level][k] >= 0; k++ ) {
            int i = mod((xArray[level][k] - x)*modOldInvert[level], residueClasses[level]);
            int xMerge = x + i * modProd[level];
            long factor = findFactors(n, factors, xMerge, level+1);
            if (factor != n)
                return factor;
        }
        return n;
    }

//    public long findFactorsByMerge(long n, Collection<Long> factors, int x, int level) {
//        int modOld = residueClasses[level - 1];
//        int mod = residueClasses[level];
//
//        final int step = (modOld / mod + 1) * mod;
//            int leftMod = x;
//            while (leftMod >= mod)
//                leftMod -= mod;
//            for (int left = x; left < modOld * mod; left += mod) {
//                if (((1L << leftMod) & xMask[level]) != 0) {
//                    long factor = findPrimeFactors(n, factors, left, level + 1);
//                    if (factor != n)
//                        return factor;
//                }
//                leftMod += mod;
//                // We avoid calculating expensive % operation
//                leftMod -= leftMod >= step ? step : (step - mod);
//            }
//        return n;
//    }

    public long findFactors(long n, Collection<Long> factors, long nOrig) {
        // to be sure we do not have some bad effects with the merging
        if (minFactor < 19)
            for (int prime : residuesArray) {
                if (n % prime == 0) {
                    factors.add((long) prime);
                    return n / prime;
                }
            }
        init(n);

        for (int i=0; xArray[0][i] >= 0; i++) {
            long factor = findFactors(n, factors, xArray[0][i],1);
            if (factor != n)
                return factor;
        }
        return n;
    }
}
