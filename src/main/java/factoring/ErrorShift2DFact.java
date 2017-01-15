package factoring;

public class ErrorShift2DFact {

	private final int maxL = 100;
	long operations = 0;
	private double exp;
	int x;
	int y;



	public long findFactor(long q, long p)
	{
		final long n = p*q;
		final long xSol = (p+q)/2;
		// factor lower factors with trial division

		float  a=1;
		final double sqrt = Math.sqrt(n);
		final int sqrtN = (int) Math.ceil(sqrt);

		exp = 2;
		int step = 100;
		int yTraget = 0;
		for (int i=1; yTraget < sqrtN; i++)
		{
			System.out.println("XX next i " + i);
			yTraget = i*step;
			int tRange = (int) Math.pow(i*step, 1/exp);
//			int tRange = 3;

			final Integer xShifted = findFactor(n, xSol, sqrtN, step, i, tRange);
			if (xShifted != null) return xShifted;
		}
		// no factor found
		return -1;
	}

	private Integer findFactor(long n, long xSol, int sqrtN, int step, int i, int tRange) {
		float a;
		for (int t = -tRange; t <= tRange; t++)
        {
            System.out.println("t:\t" + t + " \ts old \t" + s(t,i,step) + " \ts max \t" + s(t,i+1,step));
            for (int s = s(t,i,step); s <= s(t,i+1,step)-1; s++)
            {
                x = sqrtN + s;
                final long right = x * x - n;
                y = (int) Math.floor(Math.sqrt(right)) + t;
                operations++;

                    int error = (int) Math.abs(right - y * y);

                    if (error == 0) {
                        System.out.println("found with fermat");
                        return y + x;
                    }
                    a = 1;
                    if (error % 2 == 0) {
                        error = error / 2;
                        a = .5f;
                    }

                    final int sign = (int) Math.signum(right - y * y);
                    final int xShifted = x + error;
                    printSolution(xSol, a, x, t, error, sign, sqrtN, xShifted, s);
                    final long right2 = xShifted * xShifted - n;

//					System.out.format("s: \t%d t: \t%d : \t%d\n", s, t, right2);

                    operations++;
                    if (MyMath.isSquare(right2))
                    {
                        final int sqrtR = (int) Math.sqrt(right2);
                        if (sqrtR * sqrtR == right2) {
                            printSolution(xSol, a, x, t, error, sign, sqrtN, xShifted, s);
                            return (sqrtR + xShifted);
                        }
                    }
                }
        }
		return null;
	}

	private int s(int t, int i, int step) {
		int t1 = Math.abs(2 * t - 1);
		if (t1 < 8)
			return i*step/t1;
		else {
			double range = Math.pow((i) * step, 1 / exp);
			double diff = range / Math.abs(t);
			return diff==i ? 0 : 1;
		}
//		return (int) (i*step/Math.pow(t1, exp));

	}

	protected void printSolution(final long xSol, float a, int x, long t, final int error, int sign, int sqrtN, long xShifted, int s) {
		double div = 4;
		if ((0.0 + xSol-x)*div/error == Math.ceil((xSol-x)*div/error))
		{
			double rightA = (0.0 + xSol - x) / error;
			if (a==.5f)
				rightA /= 2;
//			if (rightA == 1 && a==.5f)
//				rightA = a;
			double speedup = (0.0 + xShifted - sqrtN) / operations;
			if (speedup > rightA)
			System.out.println("t:\t" + t + " \ts:\t"+ s + " right 'a' " + rightA + " \tspeedup " + speedup + " \tsearchInterval " + (x-sqrtN) + " \tx = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}
}
