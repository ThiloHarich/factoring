package factoring.fermat.lehman;

import java.util.Collection;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import factoring.trial.TrialInvFact;

/**
 * This is a version of the lehman factorization, which is a variant of the fermat
 * factorization.
 * It runs in O(n^1/3) and needs O(n^1/3) space.
 * It is about three times faster as the fastest (java) implementation I know, which is the yafu
 * lehman implementation found here https://github.com/DarkenCode/yafu/blob/master/factor/LehmanClean.c.
 * It is an extension of the {@link LehmanFactorFinder} which provides most of the performance improvements.
 * The performance related features different from LehmanFactorFinder are:
 * 1) Move the cross point when to switch from trial division to the fermat factoring towards trial division.
 *  This increases the numbers to be tested per loop and reduces the time for checking the bounds
 * 2) First use numbers for k (in x^2 - 4kn = y^2) with a factor 3, because they have more solutions then k != 3*k'.
 * Together with the arguments modulus 4 as used in LehmanFactorFinder this can be seen as an argument modulus 12, which gives the name.
 * It can also be seen as using a multiplier 12. And then going over the (more rare) solutions with multipliers different from 12
 * in a structured way from most likely to seldom multipliers.
 *
 * Both features work best if the number to be factorized have big factors.
 *
 * Why does it makes sense to do more trial divisions?
 * For a given "a" all primes below a*n^1/3 were pre calculated in O(n^1/3). After this there are only
 * a*n^1/3/log(a*n^1/3) = O(3* n^1/3 / log(n)) multiplications needed for finding all prime factors of n below a*n^1/3.
 * Since the Lehman phase takes O(n^1/3) time we can choose "a" bigger then 1 to balance both phases.
 * As "a" gets bigger, the upper bound for k goes down by a factor a^2.
 * The range for x is increased by a factor a^2 -> overall work is the same.
 * But since the range check for one k involves the expensive sqrt calculations for lookups, preventing it
 * - by reducing k max - is speeding up the algorithm.
 * And if the range (for x) gets bigger then 3 we can apply some more improvements mod 3.
 * Depending on the size of the factors "a" we can choose if we do the trial division or the lehman phase first.
 * The perfect value for "a" for bigger factors strongly depends on the implementation.
 * Here a=3.2 seems to be optimal. In this case the trial division phase will be done after the lehman phase.
 * For small numbers a is chosen low and the trial division phase is done first.
 *
 * The lehman equation x^2 - 4kn = y^2 is considered here for certain modulus.
 * Beside the arguments modulus 2 and 4 (lehman) and 32 (LehmanFactorFinder) for calculating the step for x
 * here we also use a fact modulus 3 (and 9).
 * While the modulo arguments due to powers of 2 will still keep the the sequence of k starting by 1 and then increase by one,
 * the argument for modulus = 3 will first consider numbers k = 3*k'.
 * All x might have a solution y since x^2 - 4kn = y^2 mod 3 <-> x^2 = y^2 mod 3
 * Then we investigate numbers with kn = 1 mod 3.
 * Here we know that x must be dividable by 3. There are 3 out of 9 x might have a solution y,
 * if we consider x^2 - 4kn = y^2 mod 9.
 * The numbers kn = 2 mod 3 have only 2 solutions modulo 9. I found no efficient way to iterate over them.
 * 0^2 = 0 mod 3
 * 1^2 = 1 mod 3
 * 2^2 = 1 mod 3
 * x^2 - 4kn = y^2 mod 3
 * x^2 = y^2 + kn mod 3
 * assume n != 0 mod 3 by checking if 3 divides n before applying the lehman algorithm :
 * kn = 0 mod 3 -> every x is possible -> k = 0 mod 3
 * kn = 1 mod 3 -> y^2 = 0, x^2 = 1 -> x = 1,2 mod 3 and k = n Mod 3
 * kn = 2 mod 3 -> y^2 = 1, x^2 = 0 -> x = 0 mod 3 and k = -n Mod 3
 *
 * Lehman suggests to first use numbers which have a small smooth factor like 24 = 2*2*2*3.
 * This has the following benefit:
 * (all?) the numbers fit in the range and are potentially might have a solution
 * -> we can remove the range check for those numbers
 * -> we do not have to determine the right start for the range
 * but in practice this is not performing better, why?
 *
 * We not only take numbers with a factor 3, we can also improve
 * the loop if numbers are not a multiple of 3.
 * We do this in a cascaded loop with adjusted parameters.
 *
 * Open questions, possible improvements :
 * - can we extend the mod 3 observations to mod 9 or 15? Since the the fermat formula uses k*n it seems to be
 * too complicated to use it.
 * Created by Thilo Harich on 28.06.2017.
 */
public class LehmanFactorFinderMod12 implements FactorizationOfLongs {

	static double ONE_THIRD = 1.0 / 3;
	// to be fast to decide if a number is a square we consider the
	// number modulo some values. First filer/modulus is a power of 2
	// since we can easily calculate the remainder by using bit and operation
	// here also lower modulus like 64 work fine as well
	static int mod_1024 = 1024;
	// using an additional second filter gives around 10% of speedup
	// Since we do not want to exceed the running time of  then n^1/3 ~ 20000 we
	// restrict the mod to the product of small primes below these limit
	static int mod_9_5_7_11 = 3 * 3 * 5 * 7 * 11; // 3465
	// we might extend this limit if we know we have to factorize more number we might extend
	//	static int mod_9_5_7_11 = 3*3 * 5 * 7 * 11 * 13; // 45045

	static int [] invMod9 = new int[9];


	private static boolean[] isSquareMod_1024;
	// for the other mod arguments we have to do a real expensive "%" operation
	// so we might have time to do a bit more expensive but space saving representation
	// of the squares
	private static boolean[] isSquareMod_9_5_7_11;

	float maxTrialFactorMultiplier = 1.0f;
	float maxTrialFactorMultiplierCube;
	boolean doTrialFirst;

	// This is a constant that is below 1 for rounding up double values to long
	protected static final double ROUND_UP_DOUBLE = 0.9999999665;

	double[] sqrt;
	float[] sqrtInv;
	// a fast way to do the trial division phase
	final FactorFinder smallFactoriser;
	int maxTrialFactor;

	static {
		for (int i=1; i<9; i++){
			invMod9[i] = (int) PrimeMath.invert(i, 9);
		}
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


	/**
	 * create a instance that can handle factorizations of numbers up to 2^bits
	 * This uses memory 2^(bits/3)
	 * more then 42 bits are not allowed since the precision of long does not allow more here
	 *  @param bits                  the maximum number of bits the algorithm is going to work.
	 * @param maxTrialFactorMultiplier Defines the size of the factors the algorithm is first looking for.
	 *                              The algorithm is working best for number where the second highest factor is
	 *                              maxTrialFactorMultiplier * n^1/3. If maxTrialFactorMultiplier is 1 or lower we first do trial division
	 *                              for numbers below n^1/3 then do the lehman factorizationByFactors above n^1/3.
	 *                              For maxTrialFactorMultiplier >= 1 we first inspect numbers above maxTrialFactorMultiplier * n^1/3 with the lehman
	 *                              factorization. After completing this phase we do
	 *                              trial division for numbers below maxTrialFactorMultiplier * n^1/3<br>
	 *                              Use<br>
	 *                              maxTrialFactorMultiplier <= 1<br>
	 *                              - if you do not know anything about the numbers to factorize<br>
	 *                              - if you know the second highest factor can not exceed n^1/3<br>
	 *                              here all found factors are prime factors and will be stored in the primeFactors Collection.
	 *                              maxTrialFactorMultiplier = 3.2 (this is the optimal value) <br>
	 *                              - if you know for most of the numbers the maximal or the second highest factors will exceed 3.2*n^1/3<br>
	 *                              In the last case {@link #findFactors(long, Collection)} might return a composite number
	 * @param doTrialFirst
	 */
	public LehmanFactorFinderMod12(int bits, float maxTrialFactorMultiplier, boolean doTrialFirst) {
		if (bits > 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		this.maxTrialFactorMultiplier = maxTrialFactorMultiplier < 1 ? 1 : maxTrialFactorMultiplier;
		maxTrialFactor = (int) Math.ceil(this.maxTrialFactorMultiplier * Math.pow(1L << bits, ONE_THIRD));
		maxTrialFactorMultiplierCube = this.maxTrialFactorMultiplier * this.maxTrialFactorMultiplier * this.maxTrialFactorMultiplier;
		smallFactoriser = new TrialInvFact(maxTrialFactor);
		this.doTrialFirst = doTrialFirst;
		initSquares();
	}

	protected void initSquares() {
		// precalculate the square of all possible multipliers. This takes at most n^1/3
		final int kMax = (int) (maxTrialFactor / maxTrialFactorMultiplierCube);

		sqrt = new double[kMax + 10];
		sqrtInv = new float[kMax + 10];
		for (int i = 1; i < sqrt.length; i++) {
			final double sqrtI = Math.sqrt(i);
			sqrt[i] = sqrtI;
			// Since a multiplication is
			// faster then the division we also calculate the square root and the inverse
			sqrtInv[i] = (float) (1.0 / sqrtI);
		}
		System.out.println("sqrt tables for max trial factor " + maxTrialFactorMultiplier + " built bytes used : " + sqrt.length);
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (n > 1l << 41)
			throw new IllegalArgumentException("numbers above 41 bits can not be factorized");
		// with this implementation the lehman part is not slower then the trial division
		// we only apply the multiplier if we want to cut down the numbers to be tested
		maxTrialFactor = (int) Math.ceil(maxTrialFactorMultiplier * Math.pow(n, ONE_THIRD));
		//		maxTrialFactor =  (int) Math.ceil(Math.pow(nOrig, ONE_THIRD));
		if (doTrialFirst) {
			smallFactoriser.setMaxFactor(maxTrialFactor);
		}
		else {
			// factor out the 3, needed for the mod 3 tricks
			// TODO check performance
			smallFactoriser.setMaxFactor(3);
		}
		// factor out all small factors if
		final long nAfterTrial = smallFactoriser.findFactors(n, primeFactors);
		if (primeFactors == null && nAfterTrial != n)
			return nAfterTrial;
		n = nAfterTrial;

		// if number is already factorized return immediately without calling lehman
		if (n == 1)
			return n;

		// re-adjust the maximal factor we have to search for. If small factors were found by trial division,
		// which is quite often the case for random numbers, this cuts down the runtime dramatically.
		// TODO maybe we can do even better here?
		maxTrialFactor = (int) Math.ceil(maxTrialFactorMultiplier * Math.pow(n, ONE_THIRD));
		// effectively kMax is reduced by maxFactorMultiplier^2
		final int kMax = (int) (Math.ceil(maxTrialFactor / maxTrialFactorMultiplierCube));
		final long n4 = 4 * n;
		final long n12 = 3*n4;
		final double sqrtN = Math.sqrt(n);
		final int nMod4 = (int) (n % 4);
		// we calculate n^1/6 = n^1/3*n^1/3 / n^1/2 -> we no not need to use Math
		// range is increased by maxFactorMultiplier^2
		final long nPow2Third = ((long) maxTrialFactor) * ((long) maxTrialFactor);
		final float nPow1Sixth = (float) ((nPow2Third / 4) / sqrtN);
		final double sqrt4N = 2 * sqrtN;
		final int nMod3 =  (int)  (n % 3);
		// we want kn = 4kn = x mod 3, x=2-> 2*1 = 2 , x=1 ->  2*2 = 1, 1*1 = 1,
		final int [] ksMod3 = {3,3-nMod3,nMod3};
		for (int i=0; i <3; i++)
		{
			final int kBegin = ksMod3[i];
			// for kn=2 mod 3 the step for x can be multiplied by 3
			final int xStepMultiplier = i == 1 ? 3 : 1;
			//            int xStepMultiplier = 1;
			long kn4 = kBegin * n4;
			for (int k = kBegin; k <= kMax; k+=3) {
				final double sqrt4kn = sqrt[k] * sqrt4N;
				// adding a small constant to avoid rounding issues and rounding up is much slower then
				// using the downcast and adding a constant close to 1. Took the constant from the yafu code
				int xBegin = (int) (sqrt4kn + ROUND_UP_DOUBLE);
				// use only multiplications instead of division here
				// xRange is increased by maxFactorMultiplier^2
				final float xRange = nPow1Sixth * sqrtInv[k];
				final int xEnd = (int) (sqrt4kn + xRange);
				// instead of using a step 1 here we use the mod argument from lehman
				// to reduce the possible numbers to verify.
				int xStep = 2;
				if (k % 2 == 0 ) {
					xBegin |= 1;
				} else {
					if (k * n % 4 == 3) {
						// this can be seen by analyzing mod 32
						xStep = 8;
						xBegin = xBegin + PrimeMath.mod(7 - k * n - xBegin, 8);
					} else
					{
						xStep = 4;
						xBegin = xBegin + PrimeMath.mod(k + nMod4 - xBegin, 4);
					}
				}
				// for kn = 1 mod 3 we know x = 0 mod 3 -> try to adjust xBegin
				// the mod operation is expensive it only pays out if the xRange is wide
				// -> maxFactorMultiplier mist be big e.g. 3
				while (xStepMultiplier == 3 && xBegin <= xEnd && xBegin % 3 != 0) {
					xBegin += xStep;
				}
				//                assertEquals(k * n4,kn4);
				int right = (int) ((long)xBegin * xBegin - kn4);
				final int xStepMultiplied = xStep * xStepMultiplier;
				//				for (long x = xBegin; x <= xEnd; x += xStepMultiplied) {
				for (int x = xBegin; x <= xEnd; x += xStepMultiplied) {
					// here we start using long values
					// here the main work is done: one multiplication, one subtraction, one array access
					// in java the trick to replace the multiplication with an addition does not help
					//					final long right2 = x * x - k * n4;
					//	                assertEquals(right,right2);
					// instead of taking the square root (which is a very expensive operation)
					// and squaring it, we do some mod arguments to filter out non squares
					if (isProbableSquare(right)) {
						final long y = (long) Math.sqrt(right);
						if (y * y == right) {
							final long factor = PrimeMath.gcd(n, x - y);
							if (factor != 1) {
								if (findsPrimesOnly())
								{
									primeFactors.add(factor);
									if (n != factor)
										primeFactors.add(n / factor);
									return 1;
								}
								else
									return factor;
							}
						}
					}
					right += xStepMultiplied * (2*x + xStepMultiplied);
				}
				kn4 += n12;
			}
		}
		// if we have not found a factor we still have to do the trial division phase
		if (!doTrialFirst)
			smallFactoriser.setMaxFactor(maxTrialFactor);
		n = smallFactoriser.findFactors(n, primeFactors);

		return n;
	}

	protected static boolean isProbableSquare(final long number) {
		if (isSquareMod_1024[(int) (number & 1023)]) {
			// using mod9_5_7_11 instead of hard 3465 coded number causes double run time
			// the % is expensive, but saves ~ 10% of the time, since the sqrt takes even more time
			// another mod filter gives not gain, in the YAFU impl it is used
			if (isSquareMod_9_5_7_11[(int) (number % 3465)])
			{
				final long y = (long) Math.sqrt(number);
				return y * y == number;
			}
		}
		return false;
	}
	@Override
	public boolean findsPrimesOnly(){
		return doTrialFirst;
	}

	@Override
	public String toString() {
		return "LehmanFactorFinderMod12{" +
				"maxFactorMultiplier=" + maxTrialFactorMultiplier +
				'}';
	}
}
