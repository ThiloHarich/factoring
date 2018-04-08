package factoring.shift;

import factoring.math.PrimeMath;

public class ErrorShift2DFact {

	private final int maxL = 100;
	long operations = 0;
	private double exp;
	long x;
	int y;
	long n;
	long xSol;
	int sqrtN;
	int tRange;


	public long findFactor(long q, long p)
	{
		n = p*q;
		xSol = (p+q)/2;
		// factor lower factors with trial division

		final double sqrt = Math.sqrt(n);
		sqrtN = (int) Math.ceil(sqrt);

		int step = 5;
		for (int i=0; true; i++)
		{
//			System.out.println("XX next i " + i);

			tRange = 6;   // 6
			Long factor;

			for (int t = -tRange; t <= tRange+1; t++) {
//				System.out.println("t: \t" + t + " s begin: \t" + s(i*step, t) + " e end: \t" + (s((i+1)*step, t)-1));
				for (int s = s(i*step, t); s < s((i+1)*step, t); s++) {
					if ((factor = findFactor(s, t)) != -1) return factor;
				}
			}
			for (int s = 0; s <= 5; s++) {   // 2
				for (int t = Math.max(tRange+1, i*step/8); t < (i+1)*step/8; t++) {
					if ((factor = findFactor(s, t+1)) != -1) return factor;
					if ((factor = findFactor(s, -t)) != -1) return factor;
				}
			}
		}
	}

	private int s(int sRange, int t) {
		return sRange;
//		return (int) ((sRange+0.0)/(Math.abs(2*t-1))*tRange);
	}

	private long findFactor(int s, long t) {
//		System.out.println("t:\t" + t + " \ts \t" + s);
		float a;
		x = sqrtN + s;
		final long right = x * x - n;
		y = (int) (Math.floor(Math.sqrt(right)) + t);
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
		long xShifted = x + error;
		long factor = findFactor(s, t, a, error, sign, xShifted);
		if (factor != -1) return factor;
		if (s<=3)
		{
			error *= 2;
			xShifted = x + error;
			factor = findFactor(s, t, 2, error, sign, xShifted);
			if (factor != -1) return factor;
		}
		return -1;
	}

	private long findFactor(int s, long t, float a, int error, int sign, long xShifted) {
		printSolution(xSol, a, x, t, error, sign, sqrtN, xShifted, s);
		final long right2 = xShifted * xShifted - n;
		operations++;
		if (PrimeMath.isSquare(right2))
        {
            final int sqrtR = (int) Math.sqrt(right2);
            if (sqrtR * sqrtR == right2) {
				double speedup = (0.0 + xSol - sqrtN) / operations;
				System.out.println("XXXXXXXXX Speedup : " + speedup);
				printSolution(xSol, a, x, t, error, sign, sqrtN, xShifted, s);
                return (sqrtR + xShifted);
            }
        }
		return -1;
	}

	protected void printSolution(final long xSol, float a, long x, long t, final int error, int sign, int sqrtN, long xShifted, int s) {
		double div = 2*2*2*2* 3*3 * 5*5 * 7*7 * 11 * 13;
		if ((0.0 + xSol-x)*div/error == Math.ceil((xSol-x)*div/error))
		{
			double rightA = (0.0 + xSol - x) / error;
			if (a==.5f)
				rightA /= 2;
//			if (rightA == 1 && a==.5f)
//				rightA = a;
			double speedup = (0.0 + xSol - sqrtN) / operations;
			if (rightA == 2.)
			System.out.println("t:\t" + t + " \ts:\t"+ s + " right 'a' " + rightA+ " right '1/a' " + 1/rightA + " \tspeedup " + speedup + " \tsearchInterval " + (x-sqrtN) + " \txArray = " + x + "\te = " + error + " \ta(2x + a|e|) - sign(e) = " + a*(2*x + error - sign));
		}
	}

	public long getOperationsFermat (){
		return xSol - sqrtN;
	}
	public long getOperations (){
		return operations;
	}
}
