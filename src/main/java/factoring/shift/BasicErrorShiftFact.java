package factoring.shift;

public class BasicErrorShiftFact {

	//	long operations = 0;
	//	long xSol;
	//	long sqrtN ;

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		//		xSol = (p+q)/2;
		// factor lower factors with trial division

		final double sqrt = Math.sqrt(n);
		//		int sqrtN = (int) Math.ceil(sqrt);
		long x = (int) Math.ceil(sqrt);

		while (true) {
			final long right = x * x - n;
			final long y = (int) Math.ceil(Math.sqrt(right));
			//			operations++;

			final long error = y * y - right;
			final long x2 = x + error;
			final long y2 = x2 * x2 - n;

			//			operations++;
			final double sqrtY2=Math.sqrt(y2);
			if (sqrtY2 == (long) sqrtY2) {
				//				final long sqrtY2 = (long) Math.sqrt(y2);
				//				operations++;
				//				if (sqrtY2 * sqrtY2 == y2) {
				//					final double speedup = (0.0 + xSol - sqrtN) / operations;
				//					System.out.println("XXXXXXXXX Speedup : " + speedup);
				return (long) (sqrtY2 + x2);
				//				}
			}
			x+=1;
		}
	}



}
