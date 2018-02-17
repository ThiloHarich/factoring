package factoring.math;

import java.util.Iterator;

/**
 *
 * Looks for solutions xArray for xArray^2 - n = y^2 mod2Pow m
 *
 * we calculate f(y) = y^2 mod2Pow m and store the values 2^f(xArray) in a boolean Array. This set is represented in a long value
 * squaresMask.
 * Then we iterate over xArray and test if there is a f(y) = xArray^2 + n. This is done by testing the specific boolean value.
 *
 * since (-xArray)^2 = xArray^2 we can only have to check the first m/2 +1 different values for xArray^2 + n mod2Pow m
 *
 * @author thiloharich
 *
 */
public class FermatResiduesBooleanArray{

	/**
	 * This defines the number of iterations we need in the sieve to be faster then
	 * the creation of the sieve values. Depends on the impl to find the sieve values.
	 */
	private static final int SEARCH_COUNT = 64;
	private static int MAX_FACTOR = 10;
	protected static int[] residueClasses;


	protected boolean[] squares;
	public int mod;
	protected int resIndex;
	public int[] xArray;
	protected boolean[] xMask;



	/**
	 * The list of primes to be used as residues
	 */
//	private static int [] residuesArray = {2,3,5,7,11,13,17,23, 29, 31}; // 25, 27, 49
	public static int [] residuesArray = {2,3,5,7,11,13,17,19}; // 25, 27, 49
	public static int multiplier = -1;

	/**
	 * Compute the squares xArray^2 modulo mod2Pow. We avoid expensive % calculations by
	 * using (xArray+1)^2 - xArray^2 = 2x + 1. Since (-xArray)^2 = xArray^2 we only have to iterate until mod2Pow/2 +1
	 * @param mod
	 */
	public FermatResiduesBooleanArray(int mod)
	{
		this.mod = mod;
		initSquares(mod);
	}

	public void initSquares(int mod) {
		// TODO use (-i)^2 = i^2
		squares = new boolean[mod];
		int square = 0;
		// in squares the bit k is set if k is a square modulo mod2Pow
		for (int i = 1; i < mod+1; i += 2) {
			squares[square]= true;
			// we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
			// we might move out the last call
			square += i;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
		}
	}

	public FermatResiduesBooleanArray(int ... residueClasses2)
	{
		residueClasses = residueClasses2;
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
	public static int [] calculateResidueClasses(long searchInterval) {
		residueClasses = new int [residuesArray.length+1];
		long mod = 1;
		// be sure we have at least one thing to search
		residueClasses [0] = 4;
		int i = 0;
		for (; i < residuesArray.length && mod * 4 *residuesArray[i] <= searchInterval && i<residuesArray.length; i++)
		{
			int prime = residuesArray[i];
			int modPerPrime = prime;
			mod *= prime;
			while(modPerPrime * prime <= MAX_FACTOR)
			{
				modPerPrime *= prime;
				mod *= prime;
			}
			residueClasses [i] = modPerPrime;
		}
//		int missing = (int) (searchInterval/mod2Pow);
////		int missingFactor2 = Integer.highestOneBit(missing)/16;
//
//		if (missing > 4)
//		residueClasses[0] *= 2;
		return residueClasses;
	}


	public int getMod() {
		return mod;
	}






	/**
	 * takes the solutions xArray for the modulus mod2Pow and calculates the solutions for
	 * modulus mod2Pow*modMultiplier. This is done by testing all candidates xArray+ i*mod2Pow.
	 * @param modMultiplier
	 * @return
	 */
	public int [] xArray4ModMultiplyer(int n,int modMultiplier)
	{
//		int nMod = PrimeMath.mod2Pow(n, mod2Pow*modMultiplier);
		int[]  newX = new int [xArray.length*modMultiplier];
		for (int k = 0; xArray[k] >= 0; k++ ) {
			for (int i = 0; i < modMultiplier; i++) {
				int x2 = xArray[k]+i*mod;
//				if (x2*x2-n )
			}
		}
		return  newX;
	}

	/**
	 * Calculates the possible solutions for x and stores the possible xArray values in a
	 * bit Mask instead of an array. This mask will be used in  {@link #merge(int[], int, int)}.
	 * @param n
	 * @return
	 */
	public void initX(long n)
	{
		xMask = new boolean[mod];
		resIndex= 0;
		int nMod = PrimeMath.mod(n, mod);
		int squareMinN = nMod == 0 ? 0 : mod - nMod;
		if (squares[squareMinN]) {
			xMask[0] = true;
		}
		squareMinN += 1;
		squareMinN -= squareMinN >= mod ?  mod : 0;
//			for (int j = 1; j < 2*mod2Pow+1; j += 2) {
		// since we only want to iterate over the first half of the possible values,
		// buy using (-xArray)^2 = xArray^2, we have to handle the special cases xArray=0 and xArray=mod2Pow/2
		for (int j = 3; j < mod+1; j += 2) {
			if (squares[squareMinN]) {
				xMask[j/2] = true;
				xMask[mod - j/2] = true;
			}
			squareMinN += j;
			squareMinN -= squareMinN >= mod ? mod : 0;
		}
		// since all primes different from 2 (which call xArray) are odd this will never be executed
		if ((mod & 1) == 0 && squares[squareMinN]) {
			xMask[mod/2] = true;
		}
	}

//	/**
//	 * Returns an array of ints xArray with xArray^2 + n = y^2 modulo the product of the numbers stored in {@link #residueClasses}.
//	 * xArray=-1 ends the array.
//	 * The first xArray is calculated by calling {@link #xArray(int)} then {@link #merge(int[], int, int)} is called
//	 * @param n
//	 * @return
//	 */
//	public int[] x4ResidueClasses(int n) {
//		int modProd = residueClasses[0];
//		mod2Pow = residueClasses[0];
//		FermatResiduesBooleanArray squareMods1 = new FermatResiduesBooleanArray(residueClasses[0]);
//		int nMod = PrimeMath.mod2Pow(n, residueClasses[0]);
//		int [] x1 = squareMods1.xArray(nMod);
//		for (int i = 1; i < residueClasses.length && residueClasses[i] > 1; i++) {
//			FermatResiduesBooleanArray squareMods2 = new FermatResiduesBooleanArray(residueClasses[i]);
//			nMod = PrimeMath.mod2Pow(n, residueClasses[i]);
//			squareMods2.initX(nMod);
////			int [] x2 = squareMods2.merge(x1, squareMods1.xLength(), modProd);
//			int [] x2 = squareMods2.merge(x1, squareMods1.xLength(), mod2Pow);
//			squareMods1 = squareMods2;
////			modProd *= residueClasses[i];
//			mod2Pow *= residueClasses[i];
//			x1 = x2;
//		}
//		return x1;
//	}

	public int xLength()
	{
		return resIndex;
	}

	MergeIterator mergeIterator (int [] xOld, int xLength, int mod)
	{
		return new MergeIterator (xOld, xLength, mod);
	}

	/**
	 * combines the solutions (i*mod2Pow + xArray) and (j* this.mod2Pow + this.xArray).
	 * We calculate l=(i*mod2Pow + xArray) mod2Pow this.mod2Pow. then we check if initX[l] == true.
	 * initX is a boolean array which reflects all possible values of (j* this.mod2Pow + this.xArray)
	 *
	 * A differnt approach:
	 *
	 *
	 *
	 * @param xOld
	 * @param xLength
	 * @param mod
	 * @return
	 */
	public int [] merge(int [] xOld, int xLength, int mod)
	{
		xArray = new int[xLength  * this.mod + 1];
		resIndex = 0;
		final int step = (mod / this.mod + 1) * this.mod;
		for (int k=0; xOld[k] >= 0; k++ )
		{
			// avoid % here, maybe we can keep old xMod
			int leftMod = xOld[k];
			while (leftMod >= this.mod)
				leftMod -= this.mod;
			for (int left = xOld[k]; left <  mod*this.mod; left += mod)
			{
				if (xMask[leftMod])
					xArray[resIndex++] = left;
				leftMod += mod;
				// We avoid calculating expensive % operation
				leftMod -= leftMod >= step ? step : (step - this.mod);
			}
		}
		xArray[resIndex++] = -1;
		return xArray;
	}



	/**
	 * To avoid storing the merged candidates in a big array, we want to create them on the fly,
	 * this is done by this Iterator
	 */
	public class MergeIterator implements Iterator<Integer> {
		int [] x2;
		int x2Length;
		int mod2;
		int x;
		int xMod;
		final int step;
		int k=0;
		int modNew;
		private int result;


		public MergeIterator(int [] x2, int x2Length, int mod2) {
			this.x2 = x2;
			this.x2Length = x2Length;
			this.mod2 = mod2;
			step = (mod2 / mod + 1) * mod;
			k=0;
			x = x2[k];
			xMod = x;
			modNew = mod2 * mod;
		}

		@Override
		public boolean hasNext() {
			while(x < modNew || (x >= modNew && x2[k+1] > -1)) {
				if (xMask[xMod]) {
					result = x;
					doStep();
					return true;
				}
				doStep();
			}
			return false;
		}

		@Override
		public Integer next() {
			return result;
		}

		private void doStep() {
			if (x + mod2 < modNew) {
				x += mod2;
				xMod += mod2;
				xMod -= xMod >= step ? step : (step - mod);
			}
            else {
                k++;
                x = x2[k];
                xMod = x % mod;
            }
		}
	}
}
