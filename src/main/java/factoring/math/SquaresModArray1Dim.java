package factoring.math;

/**
 * Looks for solutions xArray for  xArray^2 - n = y^2 mod2Pow m
 *
 * for f(xArray) = xArray^2 - n mod2Pow m we
 * define a Map_n xArray^2 - n  -> xArray
 *
 * @deprecated use {@link QuadraticDiophantineModBit}
 * @author thiloharich
 *
 */
public class SquaresModArray1Dim {

	public static final int UNDEFINED = -1;

	public int[] values;
	int mod;

	public SquaresModArray1Dim(int mod)
	{
		values = new int [2*mod];
		this.mod = mod;
	}

	public static SquaresModArray1Dim mod (int mod)
	{
		SquaresModArray1Dim t = new SquaresModArray1Dim(mod);
		for (int i = 0; i < 2*mod; i++) {
			t.values[i] = -1;
		}

		int i = 0;
		while (i < mod/2) {
			t.values[2*mod(i*i, mod)] = i;
			t.values[2*mod(i*i, mod)+1] = mod - i;
			i++;
		}
		if ((mod & 1) == 0)
		{
			t.values[2*mod(i*i, mod)+1] = mod - i;
		}
		return t;
	}

	public SquaresModArray1Dim plusRight(long n)
	{
		SquaresModArray1Dim squares = new SquaresModArray1Dim(mod);
		int[] valuesPlus = new int [2*mod];
		int iPn = (int) (n);
		for (int i = 0; i < mod; i++, iPn ++) {
			if (iPn >= mod)
				iPn -= mod;
			valuesPlus[2*iPn] = values[2*i];
			valuesPlus[2*iPn + 1] = values[2*i+1];
		}
		squares.values = valuesPlus;
		return squares;
	}

	public void disjuctionRight (SquaresModArray1Dim b)
	{
		for (int i = 0; i < mod; i++) {
			if (values[2*i] >= 0 && b.values[2*i] < 0) {
				values[2*i] = -1;
				values[2*i+1] = -1;
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
