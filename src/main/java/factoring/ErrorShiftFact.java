package factoring;

public class ErrorShiftFact {

	private final int maxL = 100;

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		final long ySol = (p+q)/2;
		// factor lower factors with trial division

		for (int l = 0; l < maxL ; l++) {
			int a=1;
			if(l>0)
				a=7;
			System.out.println("next a= " + a + " for new interval");
			final double sqrt = Math.sqrt(n);
			final int sqrtN = (int) Math.ceil(sqrt);
			final double log = Math.log(n);
			final double searchInterval = 2*sqrtN/(log*log);
			for (int x=sqrtN; x < sqrtN + searchInterval; x++){
				long right = x*x - n;
				final int yFloor = (int) Math.floor(Math.sqrt(right));

				if (yFloor == 0)
					return x+yFloor;

				for (int y=yFloor; y<yFloor+1; y++)
				{
					final int error = (int) Math.abs(right - y*y);
					if ((0.0 + ySol-y)/error == Math.ceil((ySol-y)/error))
						System.out.println("x=" + x + " right 'a' " + (ySol-y)/error);
					final long xShifted = x + a*error;
					right = xShifted*xShifted - n;

					if (MyMath.isSquare(right))
					{
						return (long) (Math.sqrt(right) + xShifted);
					}
				}
			}
		}
		// no factor found
		return -1;
	}
}
