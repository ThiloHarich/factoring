package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FindPrimeFact;
import factoring.math.PrimeMath;
import factoring.trial.TrialRangeFact;

/**
 * This is a version of the lehman factorization, which is a variant of the fermat
 * factorization.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster then the java version of the yafu lehman factorization.
 * By storing the square roots of the multiplier k the range can be done faster.
 * It also uses a version of trial division, where the multiple inverse of the primes are stored.
 * So instead of a division a multiplication is needed to find out if a number is dividable
 * by a prime.
 * In the lehman algorithm in most of the cases i.e. k > n^1/3 / 16 the upper bound for x is less the
 * The lower bound plus 1. In this case at most one value of x has to be considered, but
 * the calculation of the lower and upper rage has to be done all the time.
 * Here we ignore the fact that we always increase x by 2 or 4.
 * Since calculating the ranges of the inner loop requires at least one square root
 * and a division we try to reduce the cost for calculating this by precalculating the
 * square roots for the small multipliers k and the inversion of it.
 *
 * Like in the YAFU implementation we get no speed when using smooth multipliers (for k) first.
 * This is surprising since most of the implementations use small multipliers, since they should
 * increase the chance that a created number is a square.
 *
 * The Hart variant always just one x per multiplier k, this eliminates the determination of the
 * upper bound, but using it gives no extra speed. Again why?
 *
 * We need 2^42/3 = 2^14 = 16364 locations to store the squares and the inversive. This does not fit in the L1 cache,
 * but it is accessed in a sequential way -> still efficient
 *
 * Open questions, possible improvements :
 * - can we get rid of storing the square roots? how can we calculate them efficiently?
 * - In most of the cases (more then 93%) for k only one or none the value x^2 -n needs to be calculated.
 *
 * Since this is a java implementation we use just the basic operations "/" and "%"
 * and let the JVM do the optimization here. When adapting to other languages this should be done.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanReverseFact extends FindPrimeFact {

	static double ONE_THIRD = 1.0/3;
	// to be fast to decide if a number is a square we consider the
	// number modulo some values. First filer/modulus is a power of 2
	// since we can easily calculate the remainder by using bit and operation
	// here also low mod like 64 work fine as well
	static int mod_1024 = 1024;
	// using an additional second filter gives around 10% of speedup
	// Since we do not want to exceed the running time of  then n^1/3 ~ 20000 we
	// restrict the mod to the product of small primes below these limit
	static int mod_9_5_7_11 = 3*3 * 5 * 7 * 11; // 3465
	// we might extend this limit if we know we have to factorize more number we might extend
	//	static int mod_9_5_7_11 = 3*3 * 5 * 7 * 11 * 13; // 45045

	private static boolean[] isSquareMod_1024;
	// for the other mod arguments we have to do a real expensive "%" operation
	// so we might have time to do a bit more expensive but space saving representation
	// of the squares
	private static boolean[] isSquareMod_9_5_7_11;

	float maxFactorMultiplier = 1.0f;
	float maxFactorMultiplierCube;

	// This is a constant that is below 1 for rounding up double values to long
	protected static final double ROUND_UP_DOUBLE = 0.9999999665;

	double [] sqrt;
	float [] sqrtInv;
	// a fast way to do the trial division phase
	final TrialRangeFact smallFactoriser;
	int maxTrialFactor;

	static {
		isSquareMod_1024 = new boolean[mod_1024];
		//		isSquareMod9_5_7_11 = new BitSet(mod9_5_7_11);
		isSquareMod_9_5_7_11 = new boolean[mod_9_5_7_11];
		for (int i = 0; i < mod_9_5_7_11; i++) {
			if (i < mod_1024)
				//				isSquareMod_1024[(i*i) % mod_1024] = true;
				isSquareMod_1024[(i*i) & 1023] = true;
			//			isSquareMod_9_5_7_11[(i*i) % mod_9_5_7_11] = true;
			isSquareMod_9_5_7_11[(i*i) % 3456] = true;
		}
		System.out.println("sqrt mod table built                       bytes used : " + (isSquareMod_1024.length + isSquareMod_9_5_7_11.length));
	}



	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 42 bits are not allowed since the precision of long does not allow more here
	 * @param bits the maximum number of bits the algorithm is going to work.
	 * @param maxFactorMultiplierIn Defines the size of the factors the algorithm is first looking for.
	 * The algorithm is working best for number where the lowest factor is
	 * maxFactorMultiplier * n^1/3 or higher. If maxFactorMultiplier is 1 or lower we first do trial division
	 * for numbers below n^1/3 then do the lehman factorization above n^1/3. In this case the running time is
	 * bounded by max(maximal factor, c*n^1/3 / log(n)).
	 * For maxFactorMultiplier >= 1 we first inspect numbers above maxFactorMultiplier * n^1/3 with the lehman
	 * factorization.
	 * The running time in this stage is 1/maxFactorMultiplier * n^1/3. After completing this phase we do
	 * trial division for numbers below maxFactorMultiplier * n^1/3 in time maxFactorMultiplier * n^1/3  / log(n).
	 * This means if you know that there are
	 * only factor above maxFactorMultiplier * n^1/3 the runtime is bounded by 1/maxFactorMultiplier * n^1/3.<br><br>
	 * Use<br>
	 * maxFactorMultiplier = 1<br>
	 * - if you do not know anything about the numbers to factorize<br>
	 * - if you know the maximal factors can not exceed n^1/3<br>
	 * maxFactorMultiplier = 3 <br>
	 * - if you know for most of the numbers the maximal factors will exceed 3*n^1/3<br>
	 * In the last case {@link #findPrimeFactors(long, Collection)} might return a composite number
	 * TODO fix this behaviour
	 */
	public LehmanReverseFact(int bits, float maxFactorMultiplierIn) {
		if (bits > 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		maxFactorMultiplier = maxFactorMultiplierIn < 1 ? 1 : maxFactorMultiplierIn;
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(1L << bits, ONE_THIRD));
		maxFactorMultiplierCube = maxFactorMultiplier * maxFactorMultiplier * maxFactorMultiplier;
		smallFactoriser = new TrialRangeFact(maxTrialFactor);
		initSquares();
	}

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
		final int kMax = (int) (maxTrialFactor / maxFactorMultiplierCube);

		sqrt = new double[kMax + 10];
		sqrtInv = new float[kMax + 10];
		for(int i = 1; i< sqrt.length; i++)
		{
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the square root and the inverse
			sqrtInv[i] = (float) (1.0 / sqrtI);
		}
		System.out.println("sqrt tables for max trial factor "  + maxFactorMultiplier + " built bytes used : " + sqrt.length);
	}

	@Override
	public long findPrimeFactors(long n, Collection<Long> primeFactors) {
		if (n > 1l << 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		// with this implementation the lehman part is not slower then the trial division
		// we only apply the multiplier if we want to cut down the numbers to be tested
		maxTrialFactor =  (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
		//		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
		smallFactoriser.setMaxFactor(maxTrialFactor);
		// factor out the factors around n^1/3 with trial division first is this faster?
		double nPowOneThirdFact = 1 / maxFactorMultiplier;
		final long nAfterTrial = smallFactoriser.findPrimeFactors(n, primeFactors, nPowOneThirdFact, 1.);
		if (primeFactors == null && nAfterTrial != n)
			return nAfterTrial;

		n = nAfterTrial;

		// TODO this can not be
		if (n<maxTrialFactor)
			return n;

		if (PrimeMath.isSquare(n)){
			final long x = PrimeMath.sqrt(n);
			if (x*x == n) {
				if (primeFactors == null)
					return x;
				primeFactors.add(x);
				primeFactors.add(x);
				return 1;
			}
		}
		// re-adjust the maximal factor we have to search for. If factors were found, which is quite
		// often the case for arbitrary numbers, this cuts down the runtime dramatically.
		// TODO maybe we can do even better here?
		maxTrialFactor =  (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
		//		final int kMax = maxTrialFactor;
		final int kMax = (int) (Math.ceil(maxTrialFactor / maxFactorMultiplierCube));
		final long n4 = 4 * n;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		// we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
		final long nPow2Third = ((long)maxTrialFactor) * ((long)maxTrialFactor);
		final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
		// surprisingly it gives no speedup when using k's with many prime factors as lehman suggests
		// for k=2 we know that x has to be even and results in a factor more often
		final double sqrt4N = 2 * sqrtN;
		for (int k = 1; k <= kMax; k++) {
			final double sqrt4kn = sqrt[k] * sqrt4N;
			// adding a small constant to avoid rounding issues and rounding up is much slower then
			// using the downcast and adding a constant close to 1. Took the constant from the yafu code
			int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
			// use only multiplications instead of division here
			final float xRange = nPow1Sixth * sqrtInv[k];
			final int xEnd = (int) (sqrt4kn + xRange);
			// instead of using a step 1 here we use the mod argument from lehman
			// to reduce the possible numbers to verify.
			int xStep = 2;
			if (k % 2 == 0) {
				xBegin |= 1;
			}
			else{
				xStep = 4;
				xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
				//				if (xBegin % 2 == 1)
				//					xBegin++;
			}
			// since the upper limit is the lower limit plus a number lower then 1
			// in most of the cases checking the upper limit is very helpful
			for(long x = xBegin; x <= xEnd; x+= xStep) {
				// here we start using long values
				// in java the trick to replace the multiplication with an addition does not help
				final long x2 = x * x;
				final long right = x2 -  k * n4;
				// instead of taking the square root (which is a very expensive operation)
				// and squaring it, we do some mod arguments to filter out non squares
				//				if (PrimeMath.isSquare(right)) {
				if (isProbableSquare(right)) {
					final long y = (long) Math.sqrt(right);
					if (y*y == right){
						final long factor = PrimeMath.gcd(n, x - y);
						if (factor != 1) {
							if (primeFactors == null)
								return factor;
							// since maxTrialFactor >= n^1/3 we have done the trial division first -> the factor
							// has to be a prime factor, n/factor is of size < n^2/3 and can not be a composite factor
//							if (maxTrialFactor == 1) {
//								primeFactors.add(factor);
//								if (n != factor)
//									primeFactors.add(n / factor);
//								return 1;
//							}
						}
					}
				}
				//				if(LehmanYafuFact.issq1024[(int)(right & 1023)]){
				//					//					if((LehmanYafuFact.issq4199[(int)(right%3465)]&2) != 0){
				//					//						if((LehmanYafuFact.issq4199[(int)(right%4199)]&1) != 0){
				//					final int b = (int) Math.sqrt(right + 0.9);
				//					if(b*b==right){ //square found
				//						final long B2 = LehmanYafuFact.gcd64(x+b, n);
				//						assert(B2>1);
				//						//                                if(B2>=N){ System.out.printf("theorem failure: B2=%llu N=%llu\n", B2,N); }
				//						return(B2);
				//					}
				//					//						}
				//					//					}
				//				}
			}
		}
		// if we have not found a factor we still have to do the trial division phase
		if (maxFactorMultiplier > 1)
			n = smallFactoriser.findPrimeFactors(n, primeFactors, 0, nPowOneThirdFact);

		return n;
	}

	protected static boolean isProbableSquare(final long number) {
		if (isSquareMod_1024[(int) (number & 1023)]){
			// using mod9_5_7_11 instead of hard 3456 coded number causes double run time
			// the % is expensive, but saves ~ 10% of the time, since the sqrt takes even more time
			// another mod filter gives not gain, in the YAFU impl it is used
			if (isSquareMod_9_5_7_11[(int) (number % 3456)]){

				final long y = (long) Math.sqrt(number);
				return y*y == number;
			}
		}
		return false;
	}

	@Override
	public String toString() {
		return "LehmanNoSqrtFact{" +
				"maxFactorMultiplier=" + maxFactorMultiplier +
				'}';
	}
}
