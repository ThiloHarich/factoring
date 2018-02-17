package factoring.math;

import java.util.Collection;

/**
 * Created by Thilo Harich on 26.12.2017.
 */
public class FermatResidueIter  extends FermatResiduesRecursive{

    private static FermatResidueIter[] factorizers;
    int modOldInvert;

    public static FermatResidueIter create (int[] residueClasses2, long n) {
        factorizers = new FermatResidueIter[residueClasses2.length];
        residueClasses = residueClasses2;
        long modProd = 1;
        int i = 0;
        for (; residueClasses2[i] > 0; i++) {
            int mod = residueClasses2[i];
            factorizers[i] = new FermatResidueIter(mod, modProd, n);
            if (mod > 0)
                modProd *= mod;
        }
        factorizers[i] = new FermatResidueIter(0, modProd, n);

        return factorizers[0];
    }

    public FermatResidueIter(int mod, long modProd, long n) {
        this.mod = mod;
        this.previousResidueProduct = modProd;
        if (mod > 0) {
            initSquares(mod);
            if (mod % 2 == 1)
            initX(n);
        }
        if (mod > 0 && modProd > 1)
            modOldInvert = (int) PrimeMath.invert(modProd, mod);

    }


    public void initX(long n)
    {
        // TODO put the x with x^2 - n = 0 mod2Pow mod2Pow to the front
        // i.e. check if squareMinN == 0 then exchange j/2 with the front
        int exchangePos = 0;
        xArray = new int[mod + 1];
        int nMod = PrimeMath.mod(n, mod);

        resIndex= 0;
        int squareMinN = nMod == 0 ? 0 : mod-nMod;
        if (squares[squareMinN]) {
            xArray[resIndex++] = 0;
        }
        squareMinN += 1;
        squareMinN -= squareMinN >= mod ?  mod : 0;
        for (int j = 3; j < mod+1; j += 2) {
            if (squares[squareMinN]) {
                xArray[resIndex++] = j/2;
                xArray[resIndex++] = mod - j/2;
//				if (squareMinN == 0 && resIndex >= exchangePos + 3)
//				{
//					int tmp = xArray[exchangePos];
//					xArray[exchangePos++] = xArray[resIndex-2];
//					xArray[resIndex-2] = tmp;
//					tmp = xArray[exchangePos];
//					xArray[exchangePos++] = xArray[resIndex-1];
//					xArray[resIndex-1] = tmp;
//				}
            }
            squareMinN += j;
            squareMinN -= squareMinN >= mod ? mod : 0;
        }
        if ((mod & 1) == 0 && squares[squareMinN]) {
            xArray[resIndex++] = mod/2;
        }
        xArray[resIndex] = -1;
    }

    public static FermatResidueIter create(long n) {
        int [] residueClasses = calculateResidueClasses(n/6 + 2);
        return FermatResidueIter.create(residueClasses, n);
    }

    public long findFactors(long n, Collection<Long> factors) {
        // TODO merge the xArray call in here. Avoid creating the array
//        int [] xNewArr = factorizer.xArray(mod2Pow(n,modNew));

        int nMod = PrimeMath.mod(n, mod);
        long factor;
        int squareMinN = mod - nMod;
        // TODO move this to the end since nMod != 0 -> y^2 does not divide mod2Pow
        FermatResidueIter factorizer = factorizers[1];
        if (squares[squareMinN]) {
            factor = factorizer.findFactors(n, factors, 0, 1,mod);
            if (factor != n)
                return factor;
        }
        // TODO first do a round but just check if (squareMinN == 0)
        squareMinN += 1;
        squareMinN -= squareMinN >= mod ?  mod : 0;
        for (int j = 3; j < mod+1; j += 2) {
            if (squares[squareMinN]) {
                // TODO check if /2 is slow
                factor = factorizer.findFactors(n, factors, j / 2, 1, mod);
                if (factor != n)
                    return factor;
                factor = factorizer.findFactors(n, factors, mod - j/2, 1, mod);
                if (factor != n)
                    return factor;
            }
            squareMinN += j;
            squareMinN -= squareMinN>= mod ? mod : 0;
        }
        if ((mod & 1) == 0 && squares[squareMinN]) {
            factor = factorizer.findFactors(n, factors, mod/2, 1, mod);
            if (factor != n)
                return factor;
        }
        return n;
    }

    /**
     * generates all possible x values for the given mod2Pow and calls {@link #findFactors(long, Collection, int, int, int)} for
     * x. if a real factor is returned, it will be returned, otherwise the next x is going to be generated.
     * Precondtition: squares are beeing calculated
     */
    public long findFactors(long n, Collection<Long> factors, int x, int residueIndex, int modOld) {
        if (mod == 0) {
//            int modOld = (int) previousResidueProduct;
            long xEnd = (n / 3 + 3) / 2;
            // adjust sqrtN to be a multiple of modOld
//			int xBegin = (sqrtN / modOld) * modOld;
            long xL = x;
//			if (modulus * mod2Pow  >  xRange/SEARCH_COUNT) {
            while (xL <= xEnd) {
                long right = xL * xL - n;
                if (SquaresMod.isSquare(right)) {
                    long y = PrimeMath.sqrt(right);
                    long factorHigh = xL + y;
                    long factorLow = xL - y;
                    factors.add(factorHigh);
                    return factorLow;
                }
                xL+= previousResidueProduct;
            }
            return n;
        }
        // TODO put the x with x^2 - n = 0 mod2Pow mod2Pow to the front
        long factor = findFactorsByMerge(n, factors, x, residueIndex+1, modOld);
        if (factor != n)
            return factor;
        return n;
    }

    /**
     * Here the solutions x modulo modOld , ... will be merged with the solutions
     * Modulo modNew.
     * This is done by checking if the value xSolOld + i* modOld = xSolNew modulo modNew .
     * i* modOld = xSolNew - xSolOld modulo modNew
     * i = (xSolNew - xSolOld) * modOld^-1 modulo modNew
     *
     * complexity Invert: #oldSoltuions * #newSolutions
     * complexity Merge : #oldSolutions * modNew
     *
     */
    public long findFactorsByMerge(long n, Collection<Long> factors, int x, int residueIndex, int modOld) {
        // TODO merge the xArray call in here. Avoid creating the array
//        int [] xNewArr = factorizer.xArray(mod2Pow(n,modNew));

        FermatResidueIter factorizer = factorizers[residueIndex];
//        FermatResidueIter factorizerOld = factorizers[residueIndex-1];
//        int modOld = factorizerOld.mod2Pow;
        int nMod = PrimeMath.mod(n, mod);
        long factor;
        int squareMinN = PrimeMath.mod(mod - nMod, mod);
        // TODO move this to the end since nMod != 0 -> y^2 does not divide mod2Pow
        if (squares[squareMinN]) {
            int xMerge = mergeX(x, 0, modOld);
            factor = factorizer.findFactors(n, factors, xMerge, residueIndex,mod);
            if (factor != n)
                return factor;
        }
        // TODO first do a round but just check if (squareMinN == 0)
        squareMinN += 1;
        squareMinN -= squareMinN >= mod ?  mod : 0;
        for (int j = 3; j < mod+1; j += 2) {
            if (squares[squareMinN]) {
                // TODO check if /2 is slow
                int xMerge = mergeX(x, j / 2, modOld);
                factor = factorizer.findFactors(n, factors, xMerge, residueIndex, mod);
                if (factor != n)
                    return factor;
                xMerge = mergeX(x, mod - j/2, modOld);
                factor = factorizer.findFactors(n, factors, xMerge, residueIndex, mod);
                if (factor != n)
                    return factor;
            }
            squareMinN += j;
            squareMinN -= squareMinN>= mod ? mod : 0;
        }
        if ((mod & 1) == 0 && squares[squareMinN]) {
            int xMerge = mergeX(x,  mod/2, modOld);
            factor = factorizer.findFactors(n, factors, xMerge, residueIndex, mod);
            if (factor != n)
                return factor;
        }
        return n;
    }

    public int mergeX(int x, int xNew, int modOld) {
        int i = PrimeMath.mod((xNew - x)* modOldInvert, mod);
        return x + i * modOld;
    }

}
