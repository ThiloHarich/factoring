package factoring.math;

import org.junit.Test;

/**
 * Created by Thilo Harich on 08.01.2018.
 */
public class SqrtTest {

	double [] sqrts;
	int loops = 40000;

	@Test
	public void testSqrtGen () {
		final int range = 1 << (41/3);
		sqrts = new double [range+1];

		long start, end;

		start = System.currentTimeMillis();
		multiply(range);
		end  = System.currentTimeMillis();
		System.out.println(" time mult : " + (end- start));

		start = System.currentTimeMillis();
		sqrt(range);
		end  = System.currentTimeMillis();
		System.out.println(" time sqrt  : " + (end- start));

		start = System.currentTimeMillis();
		clever(range);
		end  = System.currentTimeMillis();
		System.out.println(" time clever : " + (end- start));

	}

	private double clever(int range) {
		double prod = 1;
		final double sqrt2 = Math.sqrt(2);
		final double sqrt3 = Math.sqrt(3);
		for (int j = 0; j <loops; j++) {
			for (int i = 0; i < range; i++) {
				sqrts[i] = 0;
			}
			for (int i = 0; i < range; i++) {
				if (sqrts[i] == 0)
				{
					final double sqrt = Math.sqrt(i);
					sqrts[i] = sqrt;
					prod *= sqrt;
					final int i2 = 2*i;
					if (i2 < range) {
						if(sqrts[i2] == 0)
						{
							sqrts[i2] = sqrt * sqrt2;
							prod *= sqrts[i2];
						}
						final int i3 = 3*i;
						if (i3 < range) {
							if(sqrts[i3] == 0) {
								sqrts[i3] = sqrt * sqrt3;
								prod *= sqrts[i3];
							}
							final int i4 = 4*i;
							if (i4 < range && sqrts[i4] == 0) {
								sqrts[i4] = sqrt * 2;
								prod *= sqrts[i4];
							}

						}
					}
				}
			}
		}
		return prod;

	}
	private double cleverMask(int range) {
		double prod = 1;
		final double sqrt2 = Math.sqrt(2);
		final double sqrt3 = Math.sqrt(3);
		for (int j = 0; j <loops; j++) {
			for (int i = 0; i < range; i++) {
				sqrts[i] = 0;
			}
			for (int i = 0; i < range; i++) {
				if (sqrts[i] == 0)
				{
					final double sqrt = Math.sqrt(i);
					sqrts[i] = sqrt;
					prod *= sqrt;
					final int i2 = 2*i;
					if (i2 < range) {
						if(sqrts[i2] == 0)
						{
							sqrts[i2] = sqrt * sqrt2;
							prod *= sqrts[i2];
						}
						final int i3 = 3*i;
						if (i3 < range) {
							if(sqrts[i3] == 0) {
								sqrts[i3] = sqrt * sqrt3;
								prod *= sqrts[i3];
							}
							final int i4 = 4*i;
							if (i4 < range && sqrts[i4] == 0) {
								sqrts[i4] = sqrt * 2;
								prod *= sqrts[i4];
							}

						}
					}
				}
			}
		}
		return prod;

	}

	private double sqrt(int range) {
		double prod = 1;
		for (int j = 0; j <loops; j++) {
			for (int i = 0; i < range; i++) {
				prod *= Math.sqrt(i);
			}

		}
		return prod;
	}
	private double multiply(int range) {
		double prod = 1;
		for (int j = 0; j <loops; j++) {
			for (int i = 0; i < range; i++) {
				prod *= i*i;
			}
		}
		return prod;
	}

	@Test
	public void testPerformance ()
	{
		final int bits = 16;
		final long begin = (1L << bits) + 3534;
		final long range = 10000000L;
		final int n3 = 1 << (bits / 3);

		long start = System.currentTimeMillis();
		hart(begin, range, n3);
		long end  = System.currentTimeMillis();
		System.out.println(" time hart  : " + (end- start));

		start = System.currentTimeMillis();
		floor(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time floor : " + (end- start));

		start = System.currentTimeMillis();
		hart(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time hart  : " + (end- start));

		start = System.currentTimeMillis();
		floor(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time floor : " + (end- start));

	}

	public void hart(long begin, long range, int n3) {
		boolean found = false;
		//        long  n = begin;
		for(long n = begin; n < begin + range; n++)
		{
			final long sqrt = PrimeMath.sqrt(n);
			long x = sqrt;
			for (; x< sqrt + n3; x++) {
				final long x2 = x*x;
				final long right = x2 - n;
				if (PrimeMath.isSquare(right)) {
					found = ! found;
				}
			}
		}
	}
	public void floor(long begin, long range, int n3) {
		boolean found = false;
		final int onLevel = 3;
		//        long  n = begin;
		for(long n = begin; n < begin + range; n++)
		{
			final int sqrt = (int) PrimeMath.sqrt(n);
			final int x = sqrt;
			final long x2 = x*x;
			int right = (int) (x2 - n);
			int x21 = 2 * x + 1;
			for (; x21< (sqrt + n3)*2; x21+=2) {
				if (PrimeMath.isSquare(right)) {
					found = ! found;
				}
				right += x21;
			}
		}
	}

}
