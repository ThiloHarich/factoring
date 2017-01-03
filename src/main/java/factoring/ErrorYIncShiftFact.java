package factoring;

public class ErrorYIncShiftFact {

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
//		final int tMax = (int) Math.floor(Math.sqrt(sqrt));
		final int tMax = 200;
		final int sqrtN = (int) Math.ceil(sqrt);
		final double searchInterval = 50;
		for (int x=sqrtN; x < sqrtN + searchInterval; x++){
			final long right = x*x - n;
			final int yFloor = (int) Math.floor(Math.sqrt(right));
			operations++;

//			int t = 4;
			int s = x - sqrtN + 1;
			int t = (int) (tMax / Math.sqrt(s) + 5);
			for (int yShift=0; yShift<=t && yShift < yFloor; yShift++)
			{
				long factor = factor(n, xSol, sqrtN, x, right, yShift, yFloor + yShift);
				if (factor != -1l) return factor;
				factor = factor(n, xSol, sqrtN, x, right, yShift, yFloor - yShift);
				if (factor != -1l) return factor;
			}
		}
		// no factor found
		return -1;
	}

	private long factor(long n, long xSol, int sqrtN, int x, long right, int yShift, int y2) {
		float a;
		int error = (int) Math.abs(right - y2 *y2);

		if(error == 0) {
            System.out.println("found with fermat");
            return y2 + x;
        }
		a = 1;
		if (error % 2 == 0) {
            error = error / 2;
            a = .5f;
        }

		final int sign = (int) Math.signum(right - y2*y2);
		final long xShifted = x + error;
		printSolution(xSol, a, x, yShift, error, sign, sqrtN, xShifted);
		final long right2 = xShifted*xShifted - n;

		operations++;
		if (MyMath.isSquare(right2))
        {
            final long sqrtR = (long) Math.sqrt(right2);
            if (sqrtR*sqrtR == right2)
            {
                printSolution(xSol, a, x, yShift, error, sign, sqrtN, xShifted);
                return (sqrtR + xShifted);
            }
        }
		return -1l;
	}

	protected void printSolution(final long xSol, float a, int x, long t, final int error, int sign, int sqrtN, long xShifted) {
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
			System.out.println("right 'a' " + rightA + " \tspeedup " + speedup + " \tsearchInterval " + (x-sqrtN) + "\t t " + t + " \tx = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}
}
