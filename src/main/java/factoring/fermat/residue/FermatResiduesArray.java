package factoring.fermat.residue;

import factoring.math.PrimeMath;
import factoring.math.SquaresMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 03.01.2018.
 */
public class FermatResiduesArray extends FermatResiduesRec {

    int mod; // = 24
    int repeat = 3;
    public int xRange = mod * repeat;
    public int xRangeM1 = mod * (repeat-1);
    public int x4Mod[][];
    boolean[] xMask;

    public FermatResiduesArray(int mod) {
        this.mod = mod;
        residueClasses = new int[]{mod};
        x4Mod = new int [mod][];
        lastLevel = 1;
        // TODO we need only one level
        squares = new boolean [lastLevel][];
        for (int level=0; level < residueClasses.length; level++) {
            int residue = residueClasses[level];
            initSquares(level, residue);
        }
    }

    public void init(long n, int[] residueClasses) {
        int nMod = PrimeMath.mod(n, mod);
        if (x4Mod[nMod] != null)
            return;
        x4Mod[nMod] = new int[mod+1];
        int prod = 1;
//        modProd = new int [residueClasses.length];
//        xArray = new int [lastLevel][];
//        xMask = new boolean [lastLevel];
        for (int level=0; level < residueClasses.length; level++) {
            int mod = residueClasses[level];
//            if (level < residueClasses.length-1 || level == 0)
                initX(level, mod, n);
//            else
//                initXMask(level, mod, n);
//            if (mod > 0) {
//                modProd[level] = prod;
//                prod *= mod;
//            }
        }
    }
    public void initX(int level, int mod, long n)
    {
        int nMod = PrimeMath.mod(n, mod);
        int resIndex= 0;
        int squareMinN = nMod == 0 ? 0 : mod-nMod;

        for (int i = 1; i <2*mod; i += 2) {
            if (squares[level][squareMinN]) {
                x4Mod[nMod][resIndex++] = i / 2;
            }
            // we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
            // we might move out the last call
            squareMinN += i;
            squareMinN -= squareMinN >= mod ? (squareMinN >= 2*mod ? 2*mod : mod) : 0;
        }
        x4Mod[nMod][resIndex] = -1;
    }

//    public void initX(int level, int mod, long n)
//    {
//        int resIndex= 0;
//        int nMod = PrimeMath.mod(n, mod);
//        int squareMinN = nMod == 0 ? 0 : mod - nMod;
//        if (squares[level][squareMinN]) {
//            x4Mod[nMod][resIndex++] = 0;
//        }
//        squareMinN += 1;
//        squareMinN -= squareMinN >= mod ?  mod : 0;
////			for (int j = 1; j < 2*mod+1; j += 2) {
//        // since we only want to iterate over the first half of the possible values,
//        // buy using (-xArray)^2 = xArray^2, we have to handle the special cases xArray=0 and xArray=mod/2
//        for (int j = 3; j < mod+1; j += 2) {
//            if (squares[level][squareMinN]) {
//                x4Mod[nMod][resIndex++] = j/2;
//                x4Mod[nMod][resIndex++] = mod - j/2;
//            }
//            squareMinN += j;
//            squareMinN -= squareMinN >= mod ? mod : 0;
//        }
//        // since all primes different from 2 (which call xArray) are odd this will never be executed
//        if ((mod & 1) == 0 && squares[level][squareMinN]) {
//            x4Mod[nMod][resIndex++] = mod/2;
//        }
//        x4Mod[nMod][resIndex] = -1;
//   }


//    public void initXMask(int level, int mod, long n)
//    {
//        xMask = new boolean[mod];
//        int nMod = PrimeMath.mod(n, mod);
//        int squareMinN = nMod == 0 ? 0 : mod - nMod;
//        for (int j = 1; j < 2* mod; j += 2) {
//            if (squares[level][squareMinN]) {
//                xMask[j/2] = true;
//            }
//            squareMinN += j;
//            squareMinN -= squareMinN >= mod ? mod : 0;
//        }
//    }

    public long findFactors(long n, Collection<Long> factors, long nOrig) {
        // to be sure we do not have some bad effects with the merging
        init(n, residueClasses);
        int nMod = PrimeMath.mod(n, mod);
        // we reuse the results
//        if (x4Mod[nMod] == null)
//            merge(0, nMod);
        long sqrtN = PrimeMath.sqrt(n);
        long xBegin = sqrtN;
        xBegin = ((xBegin) / mod) * mod;
        long xEnd = xBegin + xRangeM1;
        for (int i=0; x4Mod[nMod][i] >= 0; i++) {
            long x = xEnd + x4Mod[nMod][i];
            while (x >= sqrtN) {
                long right = x * x - n;
                if (SquaresMod.isSquare(right)) {
                    long y = PrimeMath.sqrt(right);
                    long factor = PrimeMath.gcd(nOrig, x - y);
                    if (factor != 1) {
                        factors.add(factor);
                        return nOrig / factor;
                    }
                }
                x -= mod;
            }
        }
        return n;
    }


    public void merge(int level, int nMod) {
        x4Mod[nMod] = new int[modProd[level+1] * residueClasses[level+1] + 1];
        int resIndex = 0;
//        for (int l = 0; xArray[level][l] >= 0; l++ )
//        {
//            int x = xArray[level][l];
//            for (int k = 0; xArray[level+1][k] >= 0; k++ ) {
//                int i = mod((xArray[level+1][k] - x)*modOldInvert[level+1], residueClasses[level+1]);
//                int xMerge = x + i * modProd[level+1];
//                x4Mod[nMod][resIndex++] = xMerge;
//            }
//        }
        int modOld = residueClasses[level];
        int modNew = residueClasses[level+1];
        final int step = ( modOld / modNew + 1) * modNew;

        for (int l = 0; xArray[level][l] >= 0; l++ ) {
            int x = xArray[level][l];
            // avoid % here, maybe we can keep old xMod
            int xMod = x;
            while (xMod >= modNew)
                xMod -= modNew;
            while (x < modOld * modNew) {
                if (xMask[xMod]) {
                x4Mod[nMod][resIndex++] = x;
                }
                x += modOld;
                xMod += modOld;
                // We avoid calculating expensive % operation
                xMod -= xMod >= step ? step : (step - modNew);
            }
        }
        x4Mod[nMod][resIndex] = -1;
    }


}
