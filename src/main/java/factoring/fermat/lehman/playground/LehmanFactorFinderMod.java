package factoring.fermat.lehman.playground;

import java.util.Collection;

import com.google.common.annotations.VisibleForTesting;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 * We analyze the lehman fermat function
 * x^2 - 4kn = y^2 mod 16
 * the solutions for x were stored in a array nextX such that the when indexing
 * nextX[kn%16][(x+1)%16]
 *
 * analyze x^2 - kn = y^2 mod 16    No 4 as multiplier seems to work as well.
 *
 * for k or kn = 2 mod 4 there seems to be no solution
 *
 * for odd k step can be much higher then 4, by using 2? alternating steps
 */
public class LehmanFactorFinderMod implements FactorizationOfLongs {

	static double ONE_THIRD = 1.0 / 3;
	// to be fast to decide if a number is a square we consider the
	// number modulo some values. First filer/modulus is a power of 2
	// since we can easily calculate the remainder by using bit and operation
	// here also low mod like 64 work fine as well
	static int mod_1024 = 1024;
	// using an additional second filter gives around 10% of speedup
	// Since we do not want to exceed the running time of  then n^1/3 ~ 20000 we
	// restrict the mod to the product of small primes below these limit
	static int mod_9_5_7_11 = 3 * 3 * 5 * 7 * 11; // 3465
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

	double[] sqrt;
	float[] sqrtInv;
	// a fast way to do the trial division phase
	final TrialInvFact smallFactoriser;
	int maxTrialFactor;

	static {
		isSquareMod_1024 = new boolean[mod_1024];
		//		isSquareMod9_5_7_11 = new BitSet(mod9_5_7_11);
		isSquareMod_9_5_7_11 = new boolean[mod_9_5_7_11];
		for (int i = 0; i < mod_9_5_7_11; i++) {
			if (i < mod_1024)
				//				isSquareMod_1024[(i*i) % mod_1024] = true;
				isSquareMod_1024[(i * i) & 1023] = true;
			//			isSquareMod_9_5_7_11[(i*i) % mod_9_5_7_11] = true;
			isSquareMod_9_5_7_11[(i * i) % 3465] = true;
		}
		System.out.println("sqrt mod table built                       bytes used : " + (isSquareMod_1024.length + isSquareMod_9_5_7_11.length));
	}

	//    byte[] getNextX() {
	byte[][] getNextX() {
		return nextX;
	}


	@VisibleForTesting
	public
	byte[][] nextX;
	//    byte[] nextX;

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
		final int kMax = maxTrialFactor;

		sqrt = new double[kMax + 10];
		sqrtInv = new float[kMax + 10];
		// TODO try other mods like 24
		nextX = new byte[32][32];
		//        nextX = new byte[32];
		for(int i = 1; i< sqrt.length; i++)
		{
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the square root and the inverse
			sqrtInv[i] = (float) (1.0 / sqrtI);
		}
		System.out.printf(" sqrt table[0..%d] built: ", sqrt.length);
		for(int i=0; i<Math.min(5, maxTrialFactor); i++){
			System.out.printf("%f ", sqrt[i]);
		}
		System.out.printf(" ... %f %f\n", sqrt[sqrt.length-2],sqrt[sqrt.length-1]);

		final int mod = 32;
		final boolean [] squares = new boolean[32];
		for (int i = 0; i< mod; i++)
			squares [(i*i) % 32] = true;

		for (int kn = 0; kn< mod; kn++)
		{
			int lastX = 32;
			for (int x = mod-1; x >= 0; x--){
				final int right = PrimeMath.mod(x*x - 4*kn, 32);
				nextX[kn][x] = (byte) (lastX - x);
				//                nextX[kn*32 + x] = (byte) (lastX - x);
				if (squares[right]){
					lastX = x;
				}
			}
			nextX[kn][mod-1] = (byte) (lastX + 1);
			//            nextX[kn * 32 + mod-1] = (byte) (lastX + 1);
		}
	}

	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 42 bits are not allowed since the precision of long does not allow more here
	 *
	 * @param bits                  the maximum number of bits the algorithm is going to work.
	 * @param maxFactorMultiplierIn Defines the size of the factors the algorithm is first looking for.
	 *                              The algorithm is working best for number where the lowest factor is
	 *                              maxFactorMultiplier * n^1/3 or higher. If maxFactorMultiplier is 1 or lower we first do trial division
	 *                              for numbers below n^1/3 then do the lehman factorizationByFactors above n^1/3. In this case the running time is
	 *                              bounded by max(maximal factor, c*n^1/3 / log(n)).
	 *                              For maxFactorMultiplier >= 1 we first inspect numbers above maxFactorMultiplier * n^1/3 with the lehman
	 *                              factorizationByFactors. After completing this phase we do
	 *                              trial division for numbers below maxFactorMultiplier * n^1/3<br>
	 *                              Use<br>
	 *                              maxFactorMultiplier <= 1<br>
	 *                              - if you do not know anything about the numbers to factorize<br>
	 *                              - if you know the maximal factors can not exceed n^maxFactorMultiplier < n^1/3<br>
	 *                              here all found factors are prime factors and will be stored in the primeFactors Collection.
	 *                              maxFactorMultiplier = 3 (this is the optimal value) <br>
	 *                              - if you know for most of the numbers the maximal factors will exceed 3*n^1/3<br>
	 *                              In the last case {@link #findFactors(long, Collection)} might return a composite number
	 */
	public LehmanFactorFinderMod(int bits, float maxFactorMultiplierIn) {
		if (bits > 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		maxFactorMultiplier = maxFactorMultiplierIn < 1 ? 1 : maxFactorMultiplierIn;
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(1L << bits, ONE_THIRD));
		maxFactorMultiplierCube = maxFactorMultiplier * maxFactorMultiplier * maxFactorMultiplier;
		smallFactoriser = new TrialInvFact(maxTrialFactor);
		initSquares();
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (n > 1l << 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		// with this implementation the lehman part is not slower then the trial division
		// we only apply the multiplier if we want to cut down the numbers to be tested
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
		//		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
		smallFactoriser.setMaxFactor(maxTrialFactor);
		// factor out all small factors if
		if (maxFactorMultiplier <= 1) {
			final long nAfterTrial = smallFactoriser.findFactors(n, primeFactors);
			if (primeFactors == null && nAfterTrial != n)
				return nAfterTrial;
			n = nAfterTrial;
		}

		// if number is already factorized return immediately without calling lehman
		if (n == 1)
			return n;

		// re-adjust the maximal factor we have to search for. If small factors were found by trial division,
		// which is quite often the case for arbitrary numbers, this cuts down the runtime dramatically.
		// TODO maybe we can do even better here?
		maxTrialFactor = (int) Math.ceil(maxFactorMultiplier * Math.pow(n, ONE_THIRD));
		// effectively kMax is reduced by maxFactorMultiplier^2 ->
		final int kMax = (int) (Math.ceil(maxTrialFactor / maxFactorMultiplierCube));
		final long n4 = 4 * n;
		final double sqrtN = Math.sqrt(n);
		final int nMod16 = (int) (n % 32);
		// we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
		final int nPow2Third = maxTrialFactor * maxTrialFactor;
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
			// TODO do a mod argument here.

			final int knMod16 = (k * nMod16) % 32;
			//            final int knMod16 = (k * nMod16) & 15;
			final int xBeginMod16 = xBegin % 32;
			//            int xMod16 = xBegin & 15;
			xBegin += nextX[knMod16][xBeginMod16];
			final int xBeginInitial = xBegin;
			//            xBegin += nextX[knMod16 * 32 + xMod16];

			while (xBegin <= xEnd && xBegin != xBeginInitial + 32) {
				long x = xBegin;
				while (x <= xEnd) {
					// here we start using long values
					// in java the trick to replace the multiplication with an addition does not help
					final long x2 = x * x;
					final long right = x2 - k * n4;
					// instead of taking the square root (which is a very expensive operation)
					// and squaring it, we do some mod arguments to filter out non squares
					//				if (PrimeMath.isSquare(right)) {
					if (isProbableSquare(right)) {
						final long y = (long) Math.sqrt(right);
						if (y * y == right) {
							final long factor = PrimeMath.gcd(n, x - y);
							if (factor != 1) {
								if (primeFactors == null)
									return factor;
								// when maxFactorMultiplier == 1 we have done the trial division first -> the factor
								// has to be a prime factor, n/factor is of size < n^2/3 and can not be a composite factor
								if (maxFactorMultiplier == 1) {
									primeFactors.add(factor);
									if (n != factor)
										primeFactors.add(n / factor);
									return 1;
								}
							}
						}
					}
					x += 32;
				}
				//                xMod16 = (int) (x & 15);
				xBegin += nextX[knMod16][ xBegin % 32];
				//                x += nextX[knMod16 * 32 + xMod16];
			}
		}
		// if we have not found a factor we still have to do the trial division phase
		if (maxFactorMultiplier > 1)
			n = smallFactoriser.findFactors(n, primeFactors);

		return n;
	}

	protected static boolean isProbableSquare(final long number) {
		if (isSquareMod_1024[(int) (number & 1023)]) {
			// using mod9_5_7_11 instead of hard 3465 coded number causes double run time
			// the % is expensive, but saves ~ 10% of the time, since the sqrt takes even more time
			// another mod filter gives not gain, in the YAFU impl it is used
			if (isSquareMod_9_5_7_11[(int) (number % 3465)]) {

				final long y = (long) Math.sqrt(number);
				return y * y == number;
			}
		}
		return false;
	}

	@Override
	public String toString() {
		return "LehmanFactorFinder{" +
				"maxFactorMultiplier=" + maxFactorMultiplier +
				'}';
	}
}
