package factoring.math;

import static org.junit.Assert.assertEquals;

/**
 * Looks for solutions xArray for xArray^2 - n = y^2 mod2Pow m
 *
 * we calculate f(xArray) = xArray^2 mod2Pow m and store the values 2^f(xArray) in an long value squaresMask
 * we calculate g(xArray) = xArray^2-n mod2Pow m and test if there is a f(xArray)
 * by testing 2^g(xArray) & squaresMask != 0
 * we
 **
 * @author thiloharich
 *
 */
public class QuadraticDiophantineModBig extends QuadraticDiophantineModBit {

	private long[] squares;
	private static int bits = 64;
	private static int bitsLog2 = 6;

	public static QuadraticDiophantineModBit create (int mod)
	{
		if (mod <= bits)
			return new QuadraticDiophantineModBit(mod);
		else
			return new QuadraticDiophantineModBig(mod);
	}

	public static int getMod(long xBound) {
//		int mod2Pow = 1 << (PrimeMath.log2(xBound / 100) / 2);
//		if (mod2Pow < 4)
			return 16;
//		return mod2Pow;
	}

	/**
	 * Initialize the set of squares xArray^2 mod2Pow mod2Pow.
	 * Such we can easily calculate xArray^2 - n, and can check if they are in the set of squares.
	 * TODO if we have mod2Pow = p^i we can use the squares modulu p^(i-1) (without mod2Pow) to calculate them
	 * @param mod
	 */
	public QuadraticDiophantineModBig(int mod)
	{
		// TODO cache the squares for all mods -> two dim array
		squares = new long[Math.max(1, (2*mod-1) >> bitsLog2)];
		this.mod = mod;
		int square = 0;
		// in squares the bit k is set if k is a square mod2Pow mod2Pow
		for (int i = 1; i < mod+1; i += 2) {
			int index = square >> bitsLog2;
			squares[index] |= (1L << square);
			// we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
			square += i;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
		}
	}
	/**
	 * Returns an array of ints xArray with xArray^2 + n = y^2 mod2Pow mod2Pow.
	 * xArray=-1 ends the array.
	 * @param n
	 * @return
	 */
	public int [] xArray(int n)
	{
		x = new int[mod + 1];
		int nMod = PrimeMath.mod(n, mod);
//		long [] squaresMinN = new long [squares.length];
//		int nWord = nMod / bits;
//		int nIndex = nMod - nWord*bits;
//		assertEquals(nMod, nWord*bits+nIndex);
//		// shifting the squares by n is like subtracting n from the squares
//		for (int sWord = 0; sWord < squares.length; sWord++)
//		{
////			upperIndex = sWord+bits-1+nMod
////			lowerIndex = sWord+nMod
//			long lowerBits = squares[sWord] << nIndex;
//			int lowerWord = PrimeMath.mod2Pow(sWord + nWord, squares.length);
//			squaresMinN[lowerWord] |= lowerBits;
//			if (nIndex != 0) {
//				long upperBits = squares[sWord] >> (bits - nIndex);
//				int upperWord = PrimeMath.mod2Pow(lowerWord + 1, squares.length);
//				squaresMinN[upperWord] |= upperBits;
//			}
//		}
		// for debugging the squares - n
//		for (int i = 0; i < squares.length; i++) {
//			for(int j=0; j<bits; j++)
//				if((squaresMinN[i] & (1L<< j)) > 0)
//					System.out.print(i*bits + j + ",");
//		}
//		System.out.println();

		// g(xArray) = xArray^2 + n
		// g(xArray-1) = (xArray-1)^2 + n  = xArray^2 - 2x - 1 + n
		resIndex= 0;
		int squarePlusN = nMod;
		for (int j = 1; j < 2*mod+1; j += 2) {
			int word = squarePlusN >> bitsLog2;
//			int squareIndex = squarePlusN % bits;
			if (((1L << squarePlusN) & squares[word])!= 0) {
				x[resIndex++] = j/2;
			}
			// square = i^2 = (j/2)^2; since i^2 - (i-1)^2 = 2i+1 = j
			// we do not need to calculate i^2 mod2Pow mod2Pow when adding j <= 2*mod2Pow
			squarePlusN += j;
			squarePlusN -= squarePlusN>=mod ? (squarePlusN>=2*mod ? 2*mod : mod) : 0;
		}
		x[resIndex] = -1;
		return x;
	}
}
