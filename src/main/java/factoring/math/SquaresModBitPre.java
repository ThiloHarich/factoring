package factoring.math;

/**
 * Looks for solutions xArray for xArray^2 - n = y^2 mod2Pow m
 *
 * we calculate f(xArray) = xArray^2 mod2Pow m and store the values 2^f(xArray) in an long value squaresMask
 * we calculate g(xArray) = xArray^2-n mod2Pow m and test if there is a f(xArray)
 * by testing 2^g(xArray) & squaresMask != 0
 * we
 *
 * for f(xArray) = xArray^2 - n mod2Pow m we
 * define a Map_n xArray^2 - n  -> xArray
 *
 * @deprecated use {@link QuadraticDiophantineModBit}
 *
 * @author thiloharich
 *
 */
public class SquaresModBitPre {

	public long squaresMask;
	public int[] squares;
	int mod;

	public SquaresModBitPre(int mod)
	{
		this.mod = mod;
		squares = new int [mod];
		int square = 0;
		for (int i = 0; i < mod; i++) {
			squaresMask |= 1 << square;
			square += 2*i+1;
			square = square>=mod ? (square>=2*mod ? square-2*mod : square-mod) : square;
			squares [i] = square;
		}
	}


    /**
     * Returns an array of ints xArray with xArray^2 - n = y^2 mod2Pow n.
	 * xArray=-1 ends the array.
	 * @param n
     * @return
     */
	public int [] x(int n)
	{
		int[] x = new int[mod + 1];
		int resIndex= 0;
		for (int j = 0; j < mod; j++) {
			int squareMinN = squares[j] + n;
			squareMinN = squareMinN >= mod ? squareMinN - mod : squareMinN;
			if ((1 << squareMinN & squaresMask) != 0)
				x[resIndex++] = j;
		}
		x[resIndex] = -1;
		return x;
	}

}
