package factoring;

public class ErrorShiftFact {

	private final int maxL = 100;
	long operations = 0;

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		final long ySol = (q- p)/2;
		// factor lower factors with trial division

		for (int l = 0; l < maxL ; l++) {
			int a=1;
			if(l>0)
				a=7;
			System.out.println("next a = " + a + " for new interval");
			final double sqrt = Math.sqrt(n);
			final int sqrtN = (int) Math.ceil(sqrt);
			final double log = Math.log(n);
			final double searchInterval = 2*sqrtN/(log*log);
			for (int x=sqrtN; x < sqrtN + searchInterval; x++){
				long right = x*x - n;
				final int yFloor = (int) Math.floor(Math.sqrt(right));
				operations++;

				if (yFloor == 0)
					return x+yFloor;

				for (int y=yFloor; y<yFloor+1; y++)
				{
					final int error = (int) Math.abs(right - y*y);
					final int sign = (int) Math.signum(right - y*y);
					final long xShifted = x + a*error;
					printSolution(ySol, a, x, y, error, sign, sqrtN, xShifted);
					right = xShifted*xShifted - n;

					operations++;
					if (MyMath.isSquare(right))
					{
						printSolution(ySol, a, x, y, error, sign, sqrtN, xShifted);
						return (long) (Math.sqrt(right) + xShifted);
					}
				}
			}
		}
		// no factor found
		return -1;
	}

	protected void printSolution(final long ySol, int a, int x, int y, final int error, int sign, int sqrtN, long xShifted) {
		if ((0.0 + ySol-y)/error == Math.ceil((ySol-y)/error))
		{
			System.out.println("right 'a' " + (ySol-y)/error + " \tspeedup " + (0.0 + xShifted-sqrtN)/operations + " \tsearchInterval " + (x-sqrtN) + " \tx = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}
}
