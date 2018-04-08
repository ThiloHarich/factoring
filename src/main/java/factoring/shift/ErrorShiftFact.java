package factoring.shift;

import factoring.fermat.FermatFact;
import factoring.math.PrimeMath;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ErrorShiftFact extends FermatFact {

	private final int maxL = 100;
	long operations = 0;
	int sqrtN;
	long xSol;
	boolean print = false;

	public ErrorShiftFact(boolean print) {
		this.print = print;
	}

	public long findFactor(long q, long p)
	{
		final long n = p*q;
		xSol = (p+q)/2;
		// factor lower factors with trial division


		List<Long> factors = new ArrayList<>();
		findFactors(n, factors);
		return factors.get(0);
	}

	public long findFactors(long n, Collection<Long> factors) {
		float  a;
		final double sqrt = Math.sqrt(n);
		int shift = 1;
		sqrtN = (int) Math.ceil(sqrt);
		sqrtN = (int) (Math.ceil(sqrtN / shift) * shift);

		for (int x = sqrtN; x<n ; x+= shift) {
			final long right = x * x - n;
			double sqrtRight = Math.sqrt(right);
			int sqrtRightInt = (int) Math.ceil(sqrtRight);
//			int sqrtRightInt = (int) Math.floor(sqrtRight);
//			sqrtRightInt += sqrtRight -sqrtRightInt <.5 ? 1 : 0;
			for (int y = sqrtRightInt; y <= sqrtRightInt; y++) {
				if (print) operations++;

				int error = (int) Math.abs(right - y * y);

				if (error == 0) {
					if (print) System.out.println("--- found with fermat");
					double speedup = (0.0 + xSol - sqrtN) / operations;
					if (print)
//						System.out.println("XXXXXXXXX Speedup : " + speedup + "xArray mod 4 " + xArray % 4);
					printSolution(xSol, 1, x, error, 1, sqrtN, x);
					factors.add((long) (y + x));
					return n/(x+y);
				}
				a = 1;
				if (error % 2 == 0) {
//					error += 2*y + 1;
					error = error / 2;
					a *= .5f;
				}
				else
				{
//					error += 2*y + 1;
//					error = error / 2;
//					a = .5f;
				}
				final int sign = (int) Math.signum(right - y * y);

				final long xShifted = x + Math.abs(error);
//				final long xShifted = xArray - error;
				final long right2 = xShifted * xShifted - n;

				if (print)
					operations++;
//				printSolution(xSol, a, xArray, error, sign, sqrtN, xShifted);
				if (PrimeMath.isSquare(right2)) {
					final long sqrtR = (long) Math.sqrt(right2);
					if (print) operations++;
					if (sqrtR * sqrtR == right2) {
						double speedup = (0.0 + xSol - sqrtN) / operations;
						if (print) {
							printSolution(xSol, a, x, error, sign, sqrtN, xShifted);
//							System.out.println("XXXXXXXXX Speedup : " + speedup + speedup + "\txArray mod 4 : " + xArray % 4 + "\txArray shifted mod 4 : " + xShifted % 4);
						}
						factors.add((sqrtR + xShifted));
						return n/(x+y);
					}
				}
			}
		}
		return n;
	}

	protected void printSolution(final long xSol, float a, long x, final int error, int sign, int sqrtN, long xShifted) {
		double div = 2*2*2 * 3 * 5;
//		if ((0.0 + xSol-xArray)*div/error == Math.ceil((xSol-xArray)*div/error))
		{
			double rightA = (0.0 + xSol - x) / error;
			if (a==.5f)
				rightA /= 2;
//			if (rightA == 1 && a==.5f)
//				rightA = a;
			double speedup = (0.0 + xSol - sqrtN) / operations;
			ErrorShiftFact silentFact = new ErrorShiftFact(false);
			if (print)
				System.out.println(" right 'a' " + rightA+ " \tspeedup " + speedup + " \tsearchInterval " + (x-sqrtN) + " \txArray = " + x +
						"\terror = " + silentFact.findAllFactors(error) + " \t xArray = " + silentFact.findAllFactors(x) +
						" \t 2x+1 = " + silentFact.findAllFactors(2*x+1) +
			" \t 2x+1+e = " + silentFact.findAllFactors(2*x+1+error) );
		}
	}

	public long getOperationsFermat (){
		return xSol - sqrtN;
	}

}
