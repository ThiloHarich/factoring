package factoring;

public class ErrorShiftFact {

	private final int maxL = 100;
	long operations = 0;

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		final long xSol = (p+q)/2;
		// factor lower factors with trial division

		float  a=1;
		//			if(l>0)
		//				a=7;
		System.out.println("next a = " + a + " for new interval");
		final double sqrt = Math.sqrt(n);
		final int tMax = (int) Math.floor(Math.sqrt(sqrt));
		final int sqrtN = (int) Math.ceil(sqrt);
		final double searchInterval = sqrtN;
		for (int x=sqrtN; x < sqrtN + searchInterval; x++){
			final long right = x*x - n;
			final int yFloor = (int) Math.floor(Math.sqrt(right));
			operations++;

			int t = 2;
			//									int t = tMax / (x-sqrtN +1) / 7;
			t = Math.max(2, t);
			for (int y=yFloor-(t-1); y<=yFloor+t; y++)
			{
				int error = (int) Math.abs(right - y*y);

				if(error == 0)
					return y+x;
				if (error % 2 == 0) {
					error = error / 2;
					a = .5f;
				}
				else
					a = 1;

				final int sign = (int) Math.signum(right - y*y);
				final long xShifted = x + error;
				printSolution(xSol, a, x, y, error, sign, sqrtN, xShifted);
				final long right2 = xShifted*xShifted - n;

				operations++;
				if (MyMath.isSquare(right2))
				{
					final long sqrtR = (long) Math.sqrt(right2);
					if (sqrtR*sqrtR == right2)
					{
						printSolution(xSol, a, x, y, error, sign, sqrtN, xShifted);
						return (sqrtR + xShifted);
					}
				}
			}
		}
		// no factor found
		return -1;
	}

	protected void printSolution(final long xSol, float a, int x, long y, final int error, int sign, int sqrtN, long xShifted) {
		if ((0.0 + xSol-x)*4/error == Math.ceil((xSol-x)*4/error))
		{
			double rightA = (0.0 + xSol - x) / error;
			if (a==.5f)
				rightA /= 2;
			if (rightA == 1 & a==.5f)
				rightA = a;
			System.out.println("right 'a' " + rightA + " \tspeedup " + (0.0 + xShifted-sqrtN)/operations + " \tsearchInterval " + (x-sqrtN) + " \tx = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}
}
