package factoring;

public class ErrorShiftFact {

	private final int maxL = 100;
	long operations = 0;
	int sqrtN;
	long xSol;

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		xSol = (p+q)/2;
		// factor lower factors with trial division

		float  a;
		final double sqrt = Math.sqrt(n);
		int shift = 2;
		sqrtN = (int) Math.ceil(sqrt);
		sqrtN = (int) (Math.ceil(sqrtN / shift) * shift);

		for (int x = sqrtN; x < 2*sqrtN ; x+= shift) {
			final long right = x * x - n;
			double sqrtRight = Math.sqrt(right);
			int sqrtRightInt = (int) Math.ceil(sqrtRight);
//			int sqrtRightInt = (int) Math.floor(sqrtRight);
//			sqrtRightInt += sqrtRight -sqrtRightInt <.5 ? 1 : 0;
			for (int y = sqrtRightInt; y <= sqrtRightInt+1; y++) {
				//			if (Math.abs(Math.round(sqrtRight)-sqrtRight) <.1) {
				//				y =
				//			}
				operations++;

				int error = (int) Math.abs(right - y * y);

				if (error == 0) {
					System.out.println("found with fermat");
					double speedup = (0.0 + xSol - sqrtN) / operations;
					System.out.println("XXXXXXXXX Speedup : " + speedup);
					return y + x;
				}
				a = 1;
				if (error % 2 == 0) {
					error = error / 2;
					a = .5f;
				}

				final int sign = (int) Math.signum(right - y * y);
				final long xShifted = x + error;
//				final long xShifted = x - error;
				printSolution(xSol, a, x, error, sign, sqrtN, xShifted);
				final long right2 = xShifted * xShifted - n;

				operations++;
				if (MyMath.isSquare(right2)) {
					final long sqrtR = (long) Math.sqrt(right2);
					operations++;
					if (sqrtR * sqrtR == right2) {
						double speedup = (0.0 + xSol - sqrtN) / operations;
						System.out.println("XXXXXXXXX Speedup : " + speedup);
						//					printSolution(xSol, a, x, error, sign, sqrtN, xShifted);
						return (sqrtR + xShifted);
					}
				}
			}
		}
		return -1;
	}

	protected void printSolution(final long xSol, float a, long x, final int error, int sign, int sqrtN, long xShifted) {
		double div = 2*2*2;
		if ((0.0 + xSol-x)*div/error == Math.ceil((xSol-x)*div/error))
		{
			double rightA = (0.0 + xSol - x) / error;
			if (a==.5f)
				rightA /= 2;
//			if (rightA == 1 && a==.5f)
//				rightA = a;
			double speedup = (0.0 + xSol - sqrtN) / operations;
//			if (rightA < 1 || rightA == Math.round(rightA))
//				System.out.println(" right 'a' " + rightA+ " \tspeedup " + speedup + " \tsearchInterval " + (x-sqrtN) + " \tx = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}

	public long getOperationsFermat (){
		return xSol - sqrtN;
	}

}
