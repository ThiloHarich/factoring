package factoring.math;

/**
 * Looks for solutions xArray for  xArray^2 - n = y^2 mod2Pow m
 *
 * for f(xArray) = xArray^2 - n mod2Pow m we
 * define a Map_n xArray^2 - n  -> xArray
 *
 * @deprecated use {@link QuadraticDiophantineModBit}
 *
 * @author thiloharich
 *
 */
public class SquaresModArray {

	public static final int UNDEFINED = -1;

	public int[][] values;
	int mod;

	public SquaresModArray(int mod)
	{
		values = new int [mod][2];
		this.mod = mod;
	}

	public static SquaresModArray mod (int mod)
	{
		SquaresModArray t = new SquaresModArray(mod);
		for (int i = 0; i < mod; i++) {
			t.values[i][0] = -1;
			t.values[i][1] = -1;
		}
		for (int i = 0; i < mod; i++) {
			t.values[mod(i*i, mod)][0] = i;
			t.values[mod(i*i, mod)][1] = mod - i;
		}
		return t;
	}

	public SquaresModArray plusRight(long n)
	{
		SquaresModArray squares = new SquaresModArray(mod);
		int[][] valuesPlus = new int [mod][2];
		int iPn = (int) (n);
		for (int i = 0; i < mod; i++, iPn ++) {
			if (iPn >= mod)
				iPn -= mod;
			valuesPlus[iPn][0] = values[i][0];
			valuesPlus[iPn][1] = values[i][1];
		}
		squares.values = valuesPlus;
		return squares;
	}

	public void disjuctionRight (SquaresModArray b)
	{
		for (int i = 0; i < mod; i++) {
			if (values[i][0] >= 0 && b.values[i][0] < 0) {
				values[i][0] = -1;
				values[i][1] = -1;
			}
		}
	}


	public static int mod (long x, int mod)
	{
		int xMod = (int) (x % mod);
		if (xMod < 0)
			xMod += mod;
		return xMod;
	}

}
