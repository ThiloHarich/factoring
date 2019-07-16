package factoring.math;

import java.util.BitSet;

/**
 * Looks for solutions x for x^2 - n = y^2 mod m
 *
 * we calculate f(x) = x^2 mod m and store the values 2^f(x) in an bit set squares
 * we calculate g(x) = x^2-n mod m and test if there is a f(g(x))
 * by testing 2^g(xArray) & squaresMask != 0
 * we
 **
 * @author thiloharich
 *
 */
public class SquaresModBitSet {

	private int[] residueClasses;
	public final BitSet squares = new BitSet();
	private final int mod;
	private int resIndex;
	private int[] x;
	private int xMask;
	int [][] nextADiff;

	public static int getMod(long xBound) {
		final int mod = 1 << (PrimeMath.log2(xBound / 100) / 2);

		return Math.max(mod, 4);
	}

	public SquaresModBitSet(int mod, boolean odd)
	{
		this.mod = mod;
		nextADiff = new int[mod][];
		int square = 0;
		// in squares the bit k is set if k is a square modulo mod2Pow
		for (int i = 1; i < mod+1; i += 2) {
			squares.set(square);
			// we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
			square += i;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
		}
		for (int i = 1; i < mod; i++)
			calcNextADiff(i, odd);
	}

	public int getMod() {
		return mod;
	}

	/**
	 * Returns an array of ints xArray with x^2 + n = y^2 modulo mod.
	 * We take the set of all squares y^2 and calculate y^2 - n.
	 * These numbers must be in the set of the squares.
	 * The subtraction by n is done by shifting by n bits.
	 * xArray=-1 ends the array.
	 * @param n
	 * @return
	 */
	public int [] xArray(int n, boolean odd)
	{
		x = new int[mod + 1];
		final BitSet squaresMinN = new BitSet();
		squares.stream().map(i -> PrimeMath.mod(i + n, mod)).forEach(j -> squaresMinN.set(j));
		//		squaresMinN.stream().forEach(i -> System.out.print(i + ","));
		//		System.out.println();
		resIndex= 0;
		int square = 0;
		for (int j = 1; j < 2*mod+1; j += 2) {
			if (squaresMinN.get(square)&& (!odd  || (odd && ((j & 2) == 2)))) {
				x[resIndex++] = j/2;
			}
			// odd number theorem
			square += j;
			square -= square>=mod ? (square>=2*mod ? 2*mod : mod) : 0;
		}
		x[resIndex] = -1;
		return x;
	}

	public int [] nextADiff (int nMod)
	{
		return nextADiff[nMod];
	}

	private void calcNextADiff(int nMod, boolean odd) {
		int nextSol;
		int solIndex = 0;
		final int modOdd = odd ? 2*mod : mod;
		nextADiff[nMod] = new int[modOdd];
		final int[] sol = xArray(nMod, odd);
		for (int index=0; index < modOdd; index++)
		{
			nextSol = sol[solIndex];
			if (nextSol == -1)
				nextSol = sol[0] + modOdd;
			if (nextSol == index)
				solIndex++;
			nextADiff[nMod][index] = nextSol - index;
		}
	}
	public int xLength()
	{
		return resIndex;
	}

}
