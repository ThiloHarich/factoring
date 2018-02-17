package factoring.math;

import com.google.common.math.IntMath;

/**
 * Looks for solutions xArray for xArray^2 + n = y^2 mod2Pow m
 *
 * we calculate f(y) = y^2 mod2Pow m and store the values 2^f(xArray) in a Set. This set is represented in a long value
 * squaresMask.
 * Then we iterate over xArray and test if there is a f(y) = xArray^2 + n. This is done by testing the specific bit
 * in the mask.
 *
 * since (-xArray)^2 = xArray^2 there can only be m/2 +1 different values for xArray^2 + n mod2Pow m
 *
 * @author thiloharich
 *
 */
public class QuadraticDiophantineModBit {

	private int[] residueClasses;
	private long squares;
	protected int mod;
	protected int resIndex;
	protected int[] x;
	protected int xMask;

	/**
	 * The list of primes to be used as residues
	 */
	public static int [] residuesArray = {3,5,7,11,13,17,19}; // 25, 27, 49


	/**
	 * Compute the squares xArray^2 modulo mod2Pow. We avoid expensive % calculations by
	 * using (xArray+1)^2 - xArray^2 = 2x + 1. Since (-xArray)^2 = xArray^2 we only have to iterate until mod2Pow/2 +1
	 * @param mod
	 */
	public QuadraticDiophantineModBit(int mod)
	{
		this.mod = mod;
		int square = 0;
		// in squares the bit k is set if k is a square modulo mod2Pow
		for (int i = 1; i < mod+1; i += 2) {
			squares |= (1L << square);
			// we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
			// we might move out the last call
			square += i;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
		}
	}

	public QuadraticDiophantineModBit(int ... residueClasses)
	{
		this.residueClasses = residueClasses;
	}

	/**
	 * Determine the residue classes for a given search interval, such that sieving with the
	 * solutions {@link #x4ResidueClasses(int)} has maximal speed.
	 * This is we have a minimal number of solutions for the product of the generated
	 * residueClasses and still have enough sieving (~ 40) for each residue.
	 * Minimal number of solutions are created by using using prime residue with a maximal exponent, fit
	 * in one long (64 bits)
	 * Good candidates are powers of 2,3,5,7
	 * 2 -> 4,8,16,32,64
	 * 3 -> 9,27
	 * 5 -> 25
	 * 7 -> 49
	 *
	 *
	 * @param searchInterval
	 * @return
	 */
	public static int [] getResidueClasses(long searchInterval) {
		int [] residues = new int [residuesArray.length+1];
		int mod = 64; // sieving is appreciated, the highest prime factor will we represented by 2^i*3^j
		for (int i = 0; i < residuesArray.length && mod * residuesArray[i] * residuesArray[i] <= searchInterval; i++)
		{
			residues [i+1] = residuesArray[i];
			mod *= residuesArray[i];
		}
		// adjust the rest with powers of 2 and 3
		int log = PrimeMath.log2(searchInterval/mod);
		mod *= IntMath.pow(3, (log)/3);
		mod /= residues[1] > 0 ? residues[1] : 1;
		residues[1] = IntMath.pow(3, (log)/3);
		log = PrimeMath.log2(searchInterval/mod);
		residues[0] = 1 << Math.max(2, log); // be sure we have 2^i >= 4

		mod = 1;
		for (int i = 0; residues[i] > 0; i++)
		{
			mod *= residues[i];
		}
		long loop = searchInterval / mod;
//		System.out.println(loop);

		return residues;
	}


	public int getMod() {
		return mod;
	}



	/**
	 * Returns an array of ints xArray with xArray^2 + n = y^2 mod2Pow mod2Pow.
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
	public int [] xArray(int n)
	{
		x = new int[mod + 1];
		int nMod = PrimeMath.mod(n, mod);

		resIndex= 0;
		int squarePlusN = nMod;
		if (((1L << squarePlusN) & squares) != 0) {
			x[resIndex++] = 0;
		}
		squarePlusN += 1;
		squarePlusN -= squarePlusN >= mod ?  mod : 0;
//		for (int j = 3; j < 2*mod2Pow+1; j += 2) {
		for (int j = 3; j < mod+1; j += 2) {
			if (((1L << squarePlusN) & squares) != 0) {
				x[resIndex++] = j/2;
				x[resIndex++] = mod - j/2;
			}
			squarePlusN += j;
			squarePlusN -= squarePlusN>=mod ? mod : 0;
		}
		if ((mod & 1) == 0 && ((1L << squarePlusN) & squares) != 0) {
			x[resIndex++] = mod/2;
		}
		x[resIndex] = -1;
		return x;
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
		int[]  newX = new int [x.length*modMultiplier];
		for (int k=0; x[k] >= 0; k++ ) {
			for (int i = 0; i < modMultiplier; i++) {
				int x2 = x[k]+i*mod;
//				if (x2*x2-n )
			}
		}
		return  newX;
	}

	/**
	 * Does the same as {@link #xArray(int)} but stores the possible xArray values in a
	 * bit Mask instead of an array. This mask will be used in  {@link #merge(int[], int, int)}.
	 * @param n
	 * @return
	 */
	public void xMask(int n)
	{
		xMask = 0;
		resIndex= 0;
		int nMod = PrimeMath.mod(n, mod);
		int squarePlusN = nMod;
		if (((1L << squarePlusN) & squares) != 0) {
			xMask |= 1L;
		}
		squarePlusN += 1;
		squarePlusN -= squarePlusN >= mod ?  mod : 0;
//			for (int j = 1; j < 2*mod2Pow+1; j += 2) {
		for (int j = 3; j < mod+1; j += 2) {
			if (((1L << squarePlusN) & squares) != 0) {
				xMask |= (1L << j/2);
				xMask |= (1L << (mod - j/2));
			}
			squarePlusN += j;
			squarePlusN -= squarePlusN >= mod ? mod : 0;
		}
		// since all primes different from 2 (which call xArray) are odd this will never be executed
		if ((mod & 1) == 0 && ((1L << squarePlusN) & squares) != 0) {
			xMask |= (1L << mod/2);
		}
	}

	/**
	 * Returns an array of ints xArray with xArray^2 + n = y^2 modulo the product of the numbers stored in {@link #residueClasses}.
	 * xArray=-1 ends the array.
	 * The first xArray is calculated by calling {@link #xArray(int)} then {@link #merge(int[], int, int)} is called
	 * @param n
	 * @return
	 */
	public int[] x4ResidueClasses(int n) {
		int modProd = residueClasses[0];
		mod = residueClasses[0];
		QuadraticDiophantineModBit squareMods1 = new QuadraticDiophantineModBit(residueClasses[0]);
		int nMod = PrimeMath.mod(n, residueClasses[0]);
		int [] x1 = squareMods1.xArray(nMod);
		for (int i = 1; i < residueClasses.length && residueClasses[i] > 1; i++) {
			QuadraticDiophantineModBit squareMods2 = new QuadraticDiophantineModBit(residueClasses[i]);
			nMod = PrimeMath.mod(n, residueClasses[i]);
			squareMods2.xMask(nMod);
//			int [] x2 = squareMods2.merge(x1, squareMods1.xLength(), modProd);
			int [] x2 = squareMods2.merge(x1, squareMods1.xLength(), mod);
			squareMods1 = squareMods2;
//			modProd *= residueClasses[i];
			mod *= residueClasses[i];
			x1 = x2;
		}
		return x1;
	}

	public int xLength()
	{
		return resIndex;
	}

	/**
	 * combines the solutions i*mod2Pow + xArray and j* this.mod2Pow + this.xArray by
	 * calculating if i*mod2Pow + xArray mod2Pow this.mod2Pow \in initX.
	 * calculate all i*m mod2Pow this.mod2Pow and see if there are values xArray,xArray' with
	 * xArray + i*m = xArray' mod2Pow this.mod2Pow
	 *
	 * 3,6,9,12,15 ->
	 * 3,1,4,2,0
	 *
	 *  xArray + (xArray' + jm')m
	 *  xArray + (xArray' + jm')m  = xArray mod2Pow m
	 *
	 *  xArray + (xArray' + jm')m  = mod2Pow m'
	 *  xArray + xArray'm + jm'm  = mod2Pow m'
	 *  xArray + xArray'm   = mod2Pow m'
	 *
	 *
	 * @param xOld
	 * @param xLength
	 * @param mod
	 * @return
	 */
	public int [] merge(int [] xOld, int xLength, int mod)
	{
		x = new int[xLength  * this.mod + 1];
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
				if (((1L << leftMod) & xMask) != 0)
					x[resIndex++] = left;
				leftMod += mod;
				// We avoid calculating expensive % operation
				leftMod -= leftMod >= step ? step : (step - this.mod);
			}
		}
		x[resIndex++] = -1;
		return x;
	}
}
