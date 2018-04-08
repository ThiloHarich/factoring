package factoring.fermat.lehman.playground;

import java.util.Collection;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 * This is a version of the lehman factorizationByFactors, which is a variant of the fermat
 * factorizationByFactors.
 * It is twice as fats as the java version of the yafu lehman factorizationByFactors.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 * 		xBegin = ceil ( sqrt( 4kn ) )
 * xEnd   = floor( sqrt( 4kn + 1/4 * n^1/6 / sqrt(k) )
 * k = a * n^1/3 , a <= 1
 * xEnd   = floor( sqrt( 4kn + 1/4 * n^1/6 / (sqrt( a) * n^1/6) )
 * xEnd   = floor( sqrt( 4kn + 1/(4 * sqrt( a) )
 * -> fractional part ( yEnd) = < 1/(4 * sqrt( a))
 * -> fractional part ( yEnd) = < 1/(4 * sqrt( a)) < 1/2 , a=1/4
 * -> fractional part ( yEnd) = < 1/(4 * sqrt( a)) < 1   , a=1/16
 * we will run the regular lehman only until k <= n^1/3 / 16 then we
 * try to find all x' = 4k'n= 4mkn  where k' = k * m , 1 <= m <= 16 which fulfill the above condition
 * so we miss all the numbers which do not have a factor 2,3,5,7,11,13, which is seldom -> can we live without it, by increasing the search range for k?
 * for each x with a k > 1/16^2 calculate
 * fractional part (x) and assign then to buckets b_i = {x|i/m <= x <= (i+1)/m} for a given m
 * for each bucket we have a list of the numbers m_i with
 * fractional part (m_i * b_i) <= 1/4 * sqrt(16/m_i) = 1/sqrt(m_i)
 * fractional part (m_i * b_i) <= 1/sqrt(m_i) >= 1/4
 * ~ 1/3 of the numbers provide this.
 * if the number does not provide this flag it in a boolean array
 * f(x) = 1/4 * n^/1/6 * x^-1/2, lower limit 1 upper limit n^1/3
 * F(x) = 1/2 * n^/1/6 * x^1/2
 * F(n^1/3) =  1/2 * n^/1/6 * n^(1/3*1/2)
 * F(n^1/3) = 1/2 * n^/1/6 * n^(1/6)
 * F(n^1/3) = 1/2 * n^/1/3
 * F(1) = 1/2 * n^/1/6 *1^(*1/2)
 * F(1) = 1/2 * n^/1/6
 * if we just run to n^1/3 / 16
 * we just have effort 1/8 * n^1/3.
 * with the old approch we then have 15/16 n^1/3 operations -> 17/16 operations.
 *		with f(x)= 1/x^-1/2 von 0 bis 16
 *		F(x) = 2 * x^1/2
 *		F(16) = 2*1/4 = 1/2
 * -> we will strike out half of the candidates
 */
public class LehmanSmallRangeFact implements FactorizationOfLongs {

	static double ONE_THIRD = 1.0/3;

	//	float balanceTrial = 1.0f;
	//	float balanceTrialCube;

	// This is a constant that is below 1 for rounding up double values to long
	protected static final double ROUND_UP_DOUBLE = 0.9999999665;
	// TODO see if it works for values above the precision/scale of the algorithm as well
	// If not check BigInteger
	// to initialize the square roots we have to determine the maximal k

	double [] sqrt;
	double [] sqrtInv;

	final TrialInvFact smallFactoriser;
	int maxTrialFactor;

	private boolean[] inRange;

	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 48 bits are not allowed since the precision of long does not allow more here
	 * @param bits
	 */
	public LehmanSmallRangeFact(int bits, float balanceTrial) {
		//		maxTrialFactor = (int) Math.ceil(balanceTrial * Math.pow(1L << bits, ONE_THIRD));
		maxTrialFactor = (int) Math.ceil(Math.pow(1L << bits, ONE_THIRD));
		// using the trial division algorithm more doe not help
		//		this.balanceTrial = balanceTrial;
		//		balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
		//        smallFactoriser = new TrialPrimesDynamicFact(maxTrialFactor);
		smallFactoriser = new TrialInvFact(maxTrialFactor);

		initSquares();
	}

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
		//		float balanceTrialCube = balanceTrial * balanceTrial * balanceTrial;
		//		int kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
		final int kMax = (int) (Math.ceil(maxTrialFactor));

		sqrt = new double[kMax + 10];
		sqrtInv = new double[kMax + 10];
		inRange = new boolean[kMax + 10];
		for(int i = 1; i< sqrt.length; i++)
		{
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the  root and the inverse
			sqrtInv[i] = 1.0d / sqrtI;
		}
		System.out.printf(" sqrt table[0..%d] built: ", sqrt.length);
		for(int i=0; i<Math.min(5, maxTrialFactor); i++){
			System.out.printf("%f ", sqrt[i]);
		}
		System.out.printf(" ... %f %f\n", sqrt[sqrt.length-2],sqrt[sqrt.length-1]);
	}

	@Override
	public long findFactors(long nOrig, Collection<Long> primeFactors) {
		// with this implementation the lehman part is not slower then the trial division
		// we do not have to use a multiplier for the maximal factor were we apply the
		// trial division phase
		//		double maxTrialFactor =  Math.ceil(balanceTrial * Math.pow(nOrig, ONE_THIRD));
		double maxTrialFactor =  Math.ceil(Math.pow(nOrig, ONE_THIRD));
		smallFactoriser.setMaxFactor((int) maxTrialFactor);
		// factor out all small factors
		final long n = smallFactoriser.findFactors(nOrig, primeFactors);

		if (n<maxTrialFactor)
			return n;

		if (PrimeMath.isSquare(n)){
			final long x = PrimeMath.sqrt(n);
			if (x*x == n) {
				primeFactors.add(x);
				return x;
			}
		}
		// re-adjust the maximal factor we have to search for. If factors were found, which is quite
		// often the case for arbitrary numbers, this cuts down the runtime dramatically.
		//		maxTrialFactor =  balanceTrial * Math.pow(n, ONE_THIRD);
		maxTrialFactor =  Math.pow(n, ONE_THIRD);
		final int kMax = (int) (Math.ceil(maxTrialFactor));
		//		int kMax = (int) (Math.ceil(maxTrialFactor / balanceTrialCube));
		final int multiplier = 4;
		final long n4 = n * multiplier;
		final int multiplierSqrt = 2;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		final double nPow2Third = maxTrialFactor * maxTrialFactor;
		// TODO is division by 4 more efficient if we use int?
		final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
		final double sqrt4N = multiplierSqrt * sqrtN;
		final int kFirstLoop = kMax/16;
		//		final int kFirstLoop = kMax > 2000 ? kMax/16 : kMax;


		for (int k = 1; k <= kFirstLoop; k++) {
			final double sqrt4kn = sqrt[k] * sqrt4N;
			// adding a small constant to avoid rounding issues and rounding up is much slower then
			// using the downcast and adding a constant close to 1. Took the constant from the yafu code
			int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
			// use only multiplications instead of division here
			// TODO use a float here
			// since it is much bigger as SQRT (Long.MaxValue) we have to take a long
			// for k > kMax / 16 xRange is less then 1 unfortunately i can not use this fact to
			// speed up the runtime. Since the range for x is very small applying mod arguments to the x values
			// makes not much sense.
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify. But reducing the numbers by a factor
			// 2 or 4 the runtime is only reduced by  something like 10%. This might be due to the Hotspot
			// Behavior of java?
			final double xRange = nPow1Sixth * sqrtInv[k];
			final int xEnd = (int) (sqrt4kn + xRange);

			//			final int xStep = 1;
			int xStep = 2;
			if (k % 2 == 0) {
				xBegin |= 1;
			}
			//			else{
			//				xBegin &= (Integer.MAX_VALUE -1);
			//			}
			else{
				//				if (k*16 < kMax) {
				// TODO use the k%8 cases as well
				xStep = 4;
				xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
				//				}
			}
			for(long x = xBegin; x <= xEnd; x+= xStep) {
				// in java the trick to replace the multiplication with an addition does not help
				final long x2 = x * x;
				final long right = x2 -  k * n4;
				// TODO use a less restrictive is square check and apply the error shift
				// instead of taking the square root (which is a very expensive operation)
				// and squaring it, we do some mod arguments to filter out non squares
				if (PrimeMath.isSquare(right)) {
					final long y = (long) Math.sqrt(right);
					final long factor = PrimeMath.gcd(n, x - y);
					if (factor != 1) {
						primeFactors.add(factor);
						return n / factor;
						// we know that the remaining factor has to be a prime factor
						// but this gives no speedup for ramdom numbers
						//						if (n != factor)
						//							factors.add(n / factor);
						//						return 1;
					}
				}
			}
		}
		// now we handle small ranges individually. There will be no more loop.
		// since the range is smaller then one we can skip the x if x and k
		// have different residues mod 2
		for (int k = kFirstLoop+1; k <= kMax; k++) {
			final double sqrt4kn = sqrt[k] * sqrt4N;
			final int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
			final int kXMod2 = (k + xBegin) % 2;
			if (kXMod2 == 1)
				//			if (k%2 == 0 && xBegin%2==1 || k%2 == 1 && xBegin%2==0)
			{
				final double xRange = nPow1Sixth * sqrtInv[k];
				final int xEnd = (int) (sqrt4kn + xRange);
				if (xBegin <= xEnd)
				{
					final long x = xBegin;
					final long x2 = x * x;
					final long right = x2 -  k * n4;
					if (PrimeMath.isSquare(right)) {
						final long y = (long) Math.sqrt(right);
						final long factor = PrimeMath.gcd(n, x - y);
						if (factor != 1) {
							primeFactors.add(factor);
							return n / factor;
						}
					}
				}
			}
		}
		// For k > kMax/16 we know that the range is less then one
		//		for k > kMax/4 we know that the fractional part has to be less then 1/2. We can check both conditions
		//		1/(4 *sqrt(kMax/k)) = 1/b
		//		(4 *sqrt(kMax/k)) = b
		//		 sqrt(kMax/k) = b/4
		//		 kMax/k = (b/4)^2
		// b = 1/sqrt(a)  =  1/sqrt()
		// only in 1 of 4 cases we have to calculate x^2-n
		//		for (int k = kMax/4+1; k <= kMax; k++) {
		//			final double sqrt4kn = sqrt4N * sqrt[k];
		//			final long xFloor = (long) sqrt4kn;
		//			// we miss xFraction == 0 here !?
		//			final long x = xFloor + 1;
		//			// TODO do the loop for even k and odd k in an extra loop
		//			if (k % 2 == 0) {
		//				if (x % 2 == 0)
		//					continue;
		//			}
		//			else{
		//				final int step = PrimeMath.mod(k + nMod4 - x, 4);
		//				if (step > 0)
		//					continue;
		//			}
		//			final double xFranction = sqrt4kn - xFloor;
		//			if (xFranction >= .5)
		//			{
		//
		//				final long x2 = x * x;
		//				final long right = x2 -  k * n4;
		//				if (PrimeMath.isSquare(right)) {
		//					final long y = (long) Math.sqrt(right);
		//					final long factor = PrimeMath.gcd(n, x - y);
		//					if (factor != 1) {
		//						factors.add(factor);
		//						return n / factor;
		//					}
		//				}
		//			}
		//		}
		return n;
	}
}
