package factoring.math;

import static org.junit.Assert.assertEquals;

import java.math.BigInteger;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by Thilo Harich on 07.03.2018.
 */
public class ModTest {

	private static final int BITS = 32;

	long n = BigInteger.valueOf(1l << 40).nextProbablePrime().longValue();

	static long oneThird = (2^BITS + 2) / 3;
	static long inv255 = (2^BITS + 254) / 255;


	@Test
	public void testBarretPerformance ()
	{
		final Random rand = new Random();

		final int range = 10000;

		long start = System.currentTimeMillis();
		doDoubleMod(rand, range);
		long end  = System.currentTimeMillis();
		System.out.println(" time multi long    : " + (end- start));

		start = System.currentTimeMillis();
		doBrrett(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time barrett       : " + (end- start));

		start = System.currentTimeMillis();
		doDoubleRound(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time multi double  : " + (end- start));

		start = System.currentTimeMillis();
		doDoubleMod(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time multi long    : " + (end- start));

		start = System.currentTimeMillis();
		doBrrett(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time barrett       : " + (end- start));

		start = System.currentTimeMillis();
		doDoubleRound(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time multi double  : " + (end- start));

		start = System.currentTimeMillis();
		doDoubleMod(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time multi long    : " + (end- start));

		start = System.currentTimeMillis();
		doBrrett(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time barrett       : " + (end- start));

		start = System.currentTimeMillis();
		doDoubleRound(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time multi double  : " + (end- start));


	}

	@Test
	public void testMod3 ()
	{
		final Random rand = new Random();

		final int range = 30000;

		long start = System.currentTimeMillis();
		mod3Fast(rand, range);
		long end  = System.currentTimeMillis();
		System.out.println(" time mod 3    : " + (end- start));

		start = System.currentTimeMillis();
		mod3(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time mod fast : " + (end- start));

		start = System.currentTimeMillis();
		sqrt(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" sqrt          : " + (end- start));
	}
	@Test
	public void testMod255 ()
	{
		final Random rand = new Random();

		final int range = 30000;

		long start = System.currentTimeMillis();
		final long prod = mod255(rand, range);
		long end  = System.currentTimeMillis();
		System.out.println(" time mod 255      : " + (end- start));

		start = System.currentTimeMillis();
		final long prod2 = mod255MultInt(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time mod 255 inv  : " + (end- start));
		assertEquals(prod, prod2);

		start = System.currentTimeMillis();
		final long prod3 = mod255Jones (rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time mod 255 jones: " + (end- start));
		assertEquals(prod, prod2);

		start = System.currentTimeMillis();
		mod256 (rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time mod 256      : " + (end- start));
	}
	private long mod255 (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				prod += j % 255;
			}
		}
		return prod;
	}
	private long mod256 (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				prod += j >> 8;
			}
		}
		return prod;
	}
	private long mod255MultInt (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				prod += mod255MultInt(j);
			}
		}
		return prod;
	}
	private long mod255Jones (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				prod += mod255Jones(j);
			}
		}
		return prod;
	}
	private long mod255Jones(long a) {
		a = (a >> 16) + (a & 0xFFFF); /* sum base 2**16 digits */
		a = (a >>  8) + (a & 0xFF);   /* sum base 2**8 digits */
		if (a < 255) return a;
		if (a < (2 * 255)) return a - 255;
		return a - (2 * 255);
	}

	private long mod3 (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				g = j*j - n;
				prod += g % 3;
			}
		}
		return prod;
	}

	private long mod3Fast (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				g = j*j - n;
				prod += mod3(g);
			}
		}
		return prod;
	}

	private long sqrt (Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			long g = i;
			final long sqrtN = (long) Math.sqrt(n);
			for (long j = sqrtN; j < range + sqrtN; j++) {
				g = j*j - n;
				prod += Math.sqrt(g);
			}
		}
		return prod;
	}

	private long mod3(long x) {
		final long q = (oneThird * x) >> BITS;
			return x - q * 3;
	}

	private long mod255MultInt(long x) {
		final long q = (inv255 * x);
		final long xMod255 = x - q * 255;
		assertEquals(x % 255, xMod255);
		return xMod255;
	}

	private long doBrrett(Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final int [] barrettParams = PrimeMath.barrettParams(i);
			long g = i;
			for (long j = 1; j < range; j++) {
				g = PrimeMath.mod(j*j + 2, i, barrettParams[0], barrettParams[1]);
				prod = PrimeMath.mod(prod*g, i, barrettParams[0], barrettParams[1]);
			}
		}
		return prod;
	}
	private long doDoubleMod(Random rand, int range) {
		long prod = 1;
		for (long i = 1; i < range; i++) {
			final double iInv = 1.0d / i;
			long g = i;
			for (long j = 1; j < range; j++) {
				g = PrimeMath.mod(j*j + 2, i, iInv);
				prod = PrimeMath.mod(prod*g, i, iInv);
			}
		}
		return prod;
	}
	private double doDoubleRound(Random rand, int range) {
		double prod = 1;
		for (double i = 1; i < range; i++) {
			final double iInv = 1.0d / i;
			double g = i;
			for (int j = 1; j < range; j++) {
				g = PrimeMath.mod(g*g + 2, i, iInv);
				prod = PrimeMath.mod(prod*g, i, iInv);
			}
		}
		return prod;
	}
	@Test
	public void testPerformance ()
	{
		final Random rand = new Random();

		final int range = 10000;

		long start = System.currentTimeMillis();
		doThilo(rand, range);
		long end  = System.currentTimeMillis();
		System.out.println(" time thilo  : " + (end- start));

		start = System.currentTimeMillis();
		doCompare(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time nativ  : " + (end- start));

		start = System.currentTimeMillis();
		doThilo(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time thilo  : " + (end- start));

		start = System.currentTimeMillis();
		doCompare(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time till   : " + (end- start));

		start = System.currentTimeMillis();
		doThilo(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time thilo  : " + (end- start));

		start = System.currentTimeMillis();
		doCompare(rand, range);
		end  = System.currentTimeMillis();
		System.out.println(" time till   : " + (end- start));

	}

	@Test
	public void pModFixedValue ()
	{
		final Random rand = new Random();

		final int range = 80000000;

		// double has 52 bit due to squaring we can only use 26 bits
		//
		final long scaleMod = 1l << (64 - 24);
		final long scaleNumber = scaleMod;
		long fixedMod = Math.abs(rand.nextLong() / scaleMod );
		fixedMod = 1813382;
		long offset = Math.abs(rand.nextLong() / scaleNumber);
		offset = 2285584;


		for (int i = 0; i < 4; i++) {
			long start = System.currentTimeMillis();
			//			doLongArray(fixedMod, offset, range);
			long end  = System.currentTimeMillis();
			//			System.out.println(" time long array   : " + (end- start));

			start = System.currentTimeMillis();
			final long prod = (long) doDoubleDiv(fixedMod, offset, range);
			end  = System.currentTimeMillis();
			System.out.println(" time double div   : " + (end- start));

			start = System.currentTimeMillis();
			final long prod4 = doDouble4(fixedMod, offset, range);
			end  = System.currentTimeMillis();
			System.out.println(" time double unrol : " + (end- start));

			assertEquals(prod, prod4);

			start = System.currentTimeMillis();
			doDoubleArray(fixedMod, offset, range);
			end  = System.currentTimeMillis();
			System.out.println(" time double array : " + (end- start));

			start = System.currentTimeMillis();
			doDoubleCores(fixedMod, offset, range);
			end  = System.currentTimeMillis();
			System.out.println(" time double cores   : " + (end- start));

			System.out.println();
		}
	}

	private long doLongMod(long fixedMod, long offset, int range) {
		long number = offset;
		long loop = 0;
		long prod = 1;
		final long xFixed = offset;
		while ( loop < range) {
			number = (number * number + 2) % fixedMod;
			final long xDiff = Math.abs(number - xFixed);
			prod = (prod * xDiff) % fixedMod;
			loop++;
		}
		return prod;
	}

	private long doLongArray(long fixedMod, long offset, int range) {
		long number = offset;
		long loop = 0;
		final long prod = 1;
		final long xFixed = offset;
		while ( loop < range) {
			number = (number * number + 2) % fixedMod;
			//			final long xDiff = Math.abs(number - xFixed);
			//			prod = (prod * xDiff);

			//						do {
			//			final long bitsDiff = Long.numberOfLeadingZeros(fixedMod) - Long.numberOfLeadingZeros(prod)-1;
			//			if (bitsDiff >= 0)
			//				prod -= (1 << bitsDiff) * fixedMod;
			//			else if (prod > fixedMod)
			//				prod -= fixedMod;
			//						}
			//						while (Long.numberOfLeadingZeros(prod) < 32);
			loop++;
		}
		return prod;
	}

	private long doDoubleArray(long fixedMod, long offset, int range) {
		final double fixedDouble = fixedMod;
		final double fixedDoubleInv = 1.0d / fixedDouble;
		final long xFixed = 678326l;
		long loop = 0;
		long prod = 1;
		final int cores = 2;
		final long [] numbers = new long [cores ];
		final long [] xProds = new long [cores];
		final long [] xDiffs = new long [cores];
		while ( loop < range) {
			numbers [0]= PrimeMath.mod(numbers[cores-1] * numbers[cores-1] + 2l, fixedMod, fixedDoubleInv);
			for (int i = 1; i < numbers.length; i++) {
				numbers [i]= PrimeMath.mod(numbers[i-1] * numbers[i-1] + 2l, fixedMod, fixedDoubleInv);
			}
			//			numbers [0]= PrimeMath.mod(numbers[3] * numbers[3] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [1]= PrimeMath.mod(numbers[0] * numbers[0] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [2]= PrimeMath.mod(numbers[1] * numbers[1] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [3]= PrimeMath.mod(numbers[2] * numbers[2] + 2l, fixedMod, fixedDoubleInv);
			for (int i = 0; i < xDiffs.length; i++) {
				xDiffs[i] = Math.abs(numbers[i] - xFixed);
			}
			for (int i = 0; i < xProds.length; i++) {
				xProds[i] = PrimeMath.mod(xProds[i] * xDiffs[i], fixedMod, fixedDoubleInv);
			}
			prod = xProds[0];
			for (int i = 1; i < xProds.length; i++) {
				prod = PrimeMath.mod(prod * xProds[i], fixedMod, fixedDoubleInv);
			}

			loop += cores;
		}
		return prod;
	}

	private long doDouble4(long fixedMod, long offset, int range) {
		final double fixedDouble = fixedMod;
		final double fixedDoubleInv = 1.0d / fixedDouble;
		final long xFixed = fixedMod;
		long loop = 0;
		final int cores = 2;
		long numbers0 = offset;
		long numbers1 = offset;
		final long numbers2 = 0;
		final long numbers3 = offset;
		long xProd = 1;
		while ( loop < range) {
			//			numbers [0]= PrimeMath.mod(numbers[cores-1] * numbers[cores-1] + 2l, fixedMod, fixedDoubleInv);
			//						for (int i = 1; i < numbers.length; i++) {
			//							numbers [i]= PrimeMath.mod(numbers[i-1] * numbers[i-1] + 2l, fixedMod, fixedDoubleInv);
			//						}
			numbers0= PrimeMath.mod(numbers1 * numbers1 + 2l, fixedMod, fixedDoubleInv);
			numbers1= PrimeMath.mod(numbers0 * numbers0 + 2l, fixedMod, fixedDoubleInv);
			//			numbers2= PrimeMath.mod(numbers1 * numbers1 + 2l, fixedMod, fixedDoubleInv);
			//			numbers3= PrimeMath.mod(numbers2 * numbers2 + 2l, fixedMod, fixedDoubleInv);

			final long xDiffs0 = Math.abs(numbers0 - xFixed);
			final long xDiffs1 = Math.abs(numbers1 - xFixed);
			//			final long xDiffs2 = Math.abs(numbers2 - xFixed);
			//			final long xDiffs3 = Math.abs(numbers3 - xFixed);

			xProd = PrimeMath.mod(xProd * xDiffs0, fixedMod, fixedDoubleInv);
			xProd = PrimeMath.mod(xProd * xDiffs1, fixedMod, fixedDoubleInv);
			//			xProd = PrimeMath.mod(xProd * xDiffs2, fixedMod, fixedDoubleInv);
			//			xProd = PrimeMath.mod(xProd * xDiffs3, fixedMod, fixedDoubleInv);

			loop += cores;
		}
		return xProd;
	}

	private double doDoubleDiv(long fixedMod, long offset, long range) {
		final double fixedDoubleInv = 1.0d / fixedMod;
		//		double number = offset;
		//		double numberDouble = offset;
		long prod = 1;
		long loop = 0;
		long number = offset;
		final long xFixed = fixedMod;
		final long c = 2;
		while ( loop < range) {
			number = PrimeMath.mod(number * number + c, fixedMod, fixedDoubleInv);
			final long xDiff = Math.abs(number - xFixed);
			prod = PrimeMath.mod(prod * xDiff, fixedMod, fixedDoubleInv);

			//			assertEquals((long)mod, ((long)number) % fixedMod);
			//			final long mod = ((long)number) - (long)(number * fixedDoubleInv);
			loop++;
		}
		return prod;
	}
	private long doDoubleCores(long fixedMod, long offset, int range) {
		final double fixedDouble = fixedMod;
		final double fixedDoubleInv = 1.0d / fixedDouble;
		final long xFixed = 678326l;
		long loop = 0;
		final long mod = 0;
		final int cores = 8;
		final long [] numbers = new long [cores ];
		final long [] xProds = new long [cores];
		final long [] xDiffs = new long [cores];
		final boolean found = false;
		while ( loop < range) {
			for (int i = 0; i < numbers.length; i++) {
				numbers [i]= PrimeMath.mod(numbers[i] * numbers[i] + 2l, fixedMod, fixedDoubleInv);
				//				if (numbers[i] == numbers[(i+1)%cores])
				//					found = true;

			}
			//			numbers [0]= PrimeMath.mod(numbers[2] * numbers[2] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [1]= PrimeMath.mod(numbers[3] * numbers[3] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [2]= PrimeMath.mod(numbers[0] * numbers[0] + 2l, fixedMod, fixedDoubleInv);
			//			numbers [3]= PrimeMath.mod(numbers[1] * numbers[1] + 2l, fixedMod, fixedDoubleInv);
			for (int i = 0; i < xDiffs.length; i++) {
				xDiffs[i] = Math.abs(numbers[i] - xFixed);
			}
			for (int i = 0; i < xProds.length; i++) {
				xProds[i] = PrimeMath.mod(xProds[i] * xDiffs[i], fixedMod, fixedDoubleInv);
			}
			loop += cores;
		}
		return mod;
	}

	public void doThilo(final Random rand, final int range) {
		for ( long countA = 0; countA < range; countA++) {
			//			final long a = Math.abs(rand.nextLong());
			final long a = Math.abs(rand.nextInt());
			//			final long bStart = Math.abs(rand.nextLong());
			final long bStart = Math.abs(rand.nextInt());
			for (long bCount = 1; bCount < range; bCount++) {
				final long mod = PrimeMath.modByShiftAndSub(a, bStart+bCount);
				//				System.out.println(a + ",");
			}
		}
	}

	public void doCompare(final Random rand, final int range) {
		for ( long countA = 0; countA < range; countA++) {
			//			final long a = Math.abs(rand.nextLong());
			final long a = Math.abs(rand.nextInt());
			//			final long bStart = Math.abs(rand.nextLong());
			final long bStart = Math.abs(rand.nextInt());
			for (long bCount = 1; bCount < range; bCount++) {
				final long mod = a % (bStart+bCount);
			}
		}
	}

	@Test
	public void testCorrect ()
	{
		final Random rand = new Random();
		final long range = 500;
		long a,b;


		for ( long countA = 0; countA < range; countA++) {
			a = Math.abs(rand.nextLong());
			a = Math.abs(rand.nextInt(100))+1;
			//			a = 52;
			for (  long countB = 0; countB < range; countB++) {
				b = Math.abs(rand.nextLong());
				b = 71*a + Math.abs(rand.nextInt(100));
				//				b = 55;
				//            long gcd1 = PrimeMath.gcd(begin + offset, other);
				final long modAct = PrimeMath.modByShiftAndSub(b, a);
				final long modExpect = b % a;
				//            Assert.assertEquals (gcd1, gcd2);
				Assert.assertEquals (b + ", " + a, modExpect, modAct);
			}
		}
	}

}
