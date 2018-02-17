package factoring.math;

import java.util.Collection;

/**
 * TODO
 * - analyze how many solutions of the fermat residue equations exist
 * - use solutions x^2 - n = 0 mod2Pow p first for factor > 3
 * - first do Trial Division
 * - check against competitiors
 * - use bigger precision can reuse this solution
 * - use the error shift
 * - use one line factorizer
 * Created by Thilo Harich on 10.12.2017.
 */
public class FermatResidueMergeByInversion extends FermatResiduesRecursive{

    int modOldInvert;



    public static FermatResiduesRecursive create( long n, int[] residueClasses) {
        FermatResiduesRecursive fermat = new FermatResidueMergeByInversion(residueClasses[0], 0, n);
        FermatResiduesRecursive root = fermat;
        long modProd = residueClasses[0];
        int i=1;
        for (; i < residueClasses.length-1 && residueClasses[i] > 0; i++) {
            int mod = residueClasses[i];
            fermat.nextFermatResidue = new FermatResidueMergeByInversion(mod, modProd, n);
            fermat = fermat.nextFermatResidue;
            if (mod > 0)
                modProd *= mod;
        }
        fermat.nextFermatResidue = new FermatResidueMergeByInversion(modProd, n);

        return root;
    }

    public FermatResidueMergeByInversion(int mod, long modProd, long n) {
        this.mod = mod;
        this.previousResidueProduct = modProd;
        initSquares(mod);
        initX(n);
        if (mod > 0 && modProd > 0)
            modOldInvert = (int) PrimeMath.invert(modProd, mod);
    }
    public FermatResidueMergeByInversion(int mod) {
        this.mod = mod;
        initSquares(mod);
    }

    public FermatResidueMergeByInversion(long modProd, long n) {
        this.mod = 0;
        this.previousResidueProduct = modProd;
        this.sqrtN = (int) PrimeMath.sqrt(n);
    }

    public static FermatResiduesRecursive create(long n) {
        int [] residueClasses = calculateResidueClasses(n/6 + 2);
//        int [] residueClasses = {64,9,5,7};
        return FermatResidueMergeByInversion.create(n, residueClasses);
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
     * 0,4   mod2Pow 8
     * 0,3,6 mod2Pow 9
     *
     */
    public long findFactorsByMerge(long n, Collection<Long> factors, int x, int modOld) {
        // we always use the precomputed solutions for x
        for (int k = 0; xArray[k] >= 0; k++ ) {
            int i = PrimeMath.mod((xArray[k] - x)*modOldInvert, mod);
            int xMerge = x + i * modOld;
            long factor = nextFermatResidue.findFactors(n, factors, xMerge, mod);
            if (factor != n)
                return factor;
        }
        return n;
    }

    /**
     * Returns an array of ints xArray with xArray^2 - n = y^2 mod2Pow mod2Pow.
     * xArray=-1 ends the array.
     * We take calculate xArray^2 - n.
     * These numbers must be in the set of the squares.
     * We use the fact (xArray+1)^2 -n - xArray^2 - n = 2x+1 -> in each interation we add j=2x.
     * using the fact that (Mod-xArray)^2 = xArray^2 modolus mod2Pow we can stop the loop at (mod2Pow+1)/2.
     * This means j=2x <= mod2Pow so we cam subtract mod2Pow instead of using expensive % operations
     * for reducing the numbers modulus mod2Pow.
     * @param n
     * @return
     */
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

}
