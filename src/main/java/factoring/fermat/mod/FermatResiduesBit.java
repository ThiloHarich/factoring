package factoring.fermat.mod;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;
import factoring.math.SquaresMod;

import java.util.Collection;

import static factoring.math.PrimeMath.mod;

/**
 * It factorizes the number n, by
 *
 * For x^2 + 1 = y^2 mod 24 there is only two solution for x (0, 12). 24 = 2^3 * 3
 * We need a n' = n * m = -1 mod 24. m = -n mod 24 does this since:
 * 5*-5   =  -25 = -1 mod 24
 * 7*-7   =  -49 = -1 mod 24
 * 11*-11 = -121 = -1 mod 24
 *
 * Created by Thilo Harich on 29.11.2017.
 */
public class FermatResiduesBit extends FermatFact {

    private static final int MAX_FACTOR = 15;
    private final long minFactor = 3;
    private int[] residueClasses;
    protected int [] residuesArray = {2,3,5,7,11,13,17,19}; // 25, 27, 49
    private long[] squares;
    private int[][] xArray;
    private int[] modOldInvert;
    private long modProd = 1;
    private int lastLevel = 0;
    private long xEnd;


    public FermatResiduesBit() {
    }

    public void init (long n) {
        xEnd = (n / minFactor + minFactor) / 2;
        int [] residueClasses = calculateResidueClasses(xEnd);
        init(n, residueClasses);
    }

    public void initSquares(int level, int mod) {
        // TODO use (-i)^2 = i^2
        int square = 0;
        // in squares the bit k is set if k is a square modulo mod
        for (int i = 1; i < mod+1; i += 2) {
            squares[level] |= (1L << square);
            // we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
            // we might move out the last call
            square += i;
            square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
        }
    }
    public void initX(int level, int mod, long n)
    {
        // TODO put the x with x^2 - n = 0 mod mod to the front
        // i.e. check if squareMinN == 0 then exchange j/2 with the front
        int exchangePos = 0;
        xArray[level] = new int[mod + 1];
        int nMod = PrimeMath.mod(n, mod);

        int resIndex= 0;
        int squareMinN = nMod == 0 ? 0 : mod-nMod;

        for (int i = 1; i <2*mod; i += 2) {
            if (((1L << squareMinN) & squares[level]) != 0) {
                xArray[level][resIndex++] = i / 2;
            }
            // we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
            // we might move out the last call
            squareMinN += i;
            squareMinN -= squareMinN >= mod ? (squareMinN >= 2*mod ? 2*mod : mod) : 0;
        }
        xArray[level][resIndex] = -1;
    }

    public void init( long n, int[] residueClasses) {
        long modProd = 1;
        modOldInvert = new int [residueClasses.length];
        xArray = new int [lastLevel][];
        squares = new long [lastLevel];
        for (int level=0; level < residueClasses.length-1 && residueClasses[level] > 0; level++) {
            int mod = residueClasses[level];
            initSquares(level, mod);
            initX(level, mod, n);
            if (mod > 0) {
                if (modProd > 1)
                modOldInvert[level] = (int) PrimeMath.invert(modProd, mod);
                modProd *= mod;
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
        // be sure we have at least one thing to search
        residueClasses [0] = 4;
        modProd = 1;
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
                xL+= modProd;
            }
            return n;
        }
        // TODO put the x with x^2 - n = 0 mod mod to the front
        long factor = findFactorsByMerge(n, factors, x, level);
        if (factor != n)
            return factor;
        return n;
    }

    public long findFactorsByMerge(long n, Collection<Long> factors, int x, int level) {
        for (int k = 0; xArray[level][k] >= 0; k++ ) {
            int i = mod((xArray[level][k] - x)*modOldInvert[level], residueClasses[level]);
            int xMerge = x + i * residueClasses[level-1];
            long factor = findFactors(n, factors, xMerge, level+1);
            if (factor != n)
                return factor;
        }
        return n;
    }

    public long findFactors(long n, Collection<Long> factors) {
        if(n < 200)
            return super.findFactors(n, factors);

        init(n);

        for (int i=0; xArray[0][i] >= 0; i++) {
            long factor = findFactors(n, factors, xArray[0][i],1);
            if (factor != n)
                return factor;
        }
        return n;
    }
}
