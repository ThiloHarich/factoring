package factoring.math;

import com.google.common.math.IntMath;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.StreamSupport;

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
public class SquaresModBitSet {

	private int[] residueClasses;
	private BitSet squares = new BitSet();
	private int mod;
	private int resIndex;
	private int[] x;
	private int xMask;

	public static int getMod(long xBound) {
		int mod = 1 << (PrimeMath.log2(xBound / 100) / 2);

		return Math.max(mod, 4);
	}

	public SquaresModBitSet(int mod)
	{
		this.mod = mod;
		int square = 0;
		// in squares the bit k is set if k is a square modulo mod2Pow
		for (int i = 1; i < mod+1; i += 2) {
			squares.set(square);
			// we want to avoid the % and since (s+1)^2 = s^2 + 2s + 1, we always add i=2s+1
			square += i;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
		}
	}

	public int getMod() {
		return mod;
	}

	/**
	 * Returns an array of ints xArray with xArray^2 + n = y^2 mod2Pow mod2Pow.
	 * We take the set of all squares y^2 and calculate y^2 - n.
	 * These numbers must be in the set of the squares.
	 * The subtraction by n is done by shifting by n bits.
	 * xArray=-1 ends the array.
	 * @param n
	 * @return
	 */
	public int [] xArray(int n)
	{
		x = new int[mod + 1];
//		xArray^2 = y^2 - n mod2Pow mod2Pow
		// shifting the squares by n is like subtracting n from the squares
		BitSet squaresMinN = new BitSet();
		squares.stream().map(i -> PrimeMath.mod(i - n, mod)).forEach(j -> squaresMinN.set(j));
//		squaresMinN.stream().forEach(i -> System.out.print(i + ","));
//		System.out.println();
		resIndex= 0;
		int square = 0;
		for (int j = 1; j < 2*mod+1; j += 2) {
			if (squaresMinN.get(square)) {
				x[resIndex++] = j/2;
			}
			// square = i^2 = (j/2)^2; since i^2 - (i-1)^2 = 2i+1 = j
			// we do not need to calculate i^2 mod2Pow mod2Pow when adding j <= 2*mod2Pow
			square += j;
			square -= square>=mod ? (square>=2*mod ? 2*mod : mod) : 0;
		}
		x[resIndex] = -1;
		return x;
	}

	public int xLength()
	{
		return resIndex;
	}

}
