package factoring.math;

import org.junit.Test;

/**
 * Created by Thilo Harich on 08.01.2018.
 */
public class SqrtTest {

	double [] sqrts;
	double [] sqrtInvs;
	int loops = 40000;

	@Test
	public void testSqrtGen () {
		final long range = 1l << (41/3);
		//		sqrts = new double [range+1];

		sqrts = new double[(int) range];
		sqrtInvs = new double[(int) range];
		for (int i = 0; i < sqrts.length; i++) {
			sqrts[i] = Math.sqrt(i);
			sqrtInvs[i] = 1.0/Math.sqrt(i);
		}

		long start, end;

		start = System.currentTimeMillis();
		cast(range);
		end  = System.currentTimeMillis();
		System.out.println(" time cast  : " + (end- start));

		start = System.currentTimeMillis();
		floor(range);
		end  = System.currentTimeMillis();
		System.out.println(" time floor : " + (end- start));

		start = System.currentTimeMillis();
		div(range);
		end  = System.currentTimeMillis();
		System.out.println(" time div  : " + (end- start));

		start = System.currentTimeMillis();
		primemod(range);
		end  = System.currentTimeMillis();
		System.out.println(" time prime: " + (end- start));

		//		start = System.currentTimeMillis();
		//		array2(range);
		//		end  = System.currentTimeMillis();
		//		System.out.println(" time array2 : " + (end- start));
		//
		//		start = System.currentTimeMillis();
		//		array(range);
		//		end  = System.currentTimeMillis();
		//		System.out.println(" time array  : " + (end- start));
		//
	}

	private double array(long range) {
		double prod = 1;
		final long n = 1l << 21;
		final double sqrtN =  Math.sqrt(n);
		for (int j = 1; j <loops; j++) {
			for (int i = 1; i < range; i++) {
				final double sqrt = sqrts[i] * sqrtN;
				prod += sqrt  + 1d/sqrt ;
			}

		}
		return prod;
	}
	private double array2(long range) {
		double prod = 1;
		final long n = 1l << 21;
		final double sqrtN =  Math.sqrt(n);
		final double sqrtNInv =  1d/Math.sqrt(n);
		for (int j = 1; j <loops; j++) {
			for (int i = 1; i < range; i++) {
				final double sqrt = sqrts[i] * sqrtN;
				final double sqrtInv = sqrtInvs[i] * sqrtNInv;
				prod += sqrt  + sqrtInv ;
			}
		}
		return prod;
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

	private double sqrt(long range) {
		double prod = 1;
		final long n = 1l << 21;
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final long sqrt = (long) Math.sqrt(i*n);
				prod += sqrt;
			}
		}
		return prod;
	}

	private double primemod(long range) {
		double prod = 1;
		final double inv63 = 1.0 / 63;
		final long n = 1l << 21;
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final long mod = PrimeMath.mod(i*n, 63, inv63 );
				prod += mod;
			}
		}
		return prod;
	}

	private double mod(long range) {
		long prod = 1;
		final long n = 1l << 21;
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final long mod = i*n % 8190;
				prod += mod;
			}
		}
		return prod;
	}

	private double cast(long range) {
		final double prod = .001324;
		final long n = 1l << 21;
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final long mod = (long) (i*j*prod);
			}
		}
		return prod;
	}
	private double floor(long range) {
		final double prod = .001324;
		final long n = 1l << 21;
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final double mod = Math.floor(i*j*prod);
			}
		}
		return prod;
	}

	private double div(long range) {
		double prod = 1;
		final long bigVal = 1l << 40 + 78931;
		final long n = 1l << 21;
		for (int j = 1; j <loops; j++) {
			for (long i = 1; i < range; i++) {
				final long sqrt = bigVal/(i*n);
				prod += sqrt;
			}

		}
		return prod;
	}

	/**
	 * 1 -> 1,2
	 * 3 -> 3,6
	 * 4
	 * 5 -> 5, 10
	 * 7 -> 7, 14
	 * 8
	 * 9 -> 9, 18
	 * 11->11, 22
	 * 12
	 *
	 * @param range
	 * @return
	 */
	private double sqrt2(long range) {
		double prod = 1;
		final long n = 1l << 21;
		final double sqrt2 = Math.sqrt(2);
		final double sqrt2Inv = 1d / sqrt2;
		for (int j = 1; j <loops; j++) {
			for (long i = 1; i < range; i +=2) {
				final double sqrt = Math.sqrt(i*n);
				final double sqrtInv = 1.0/sqrt;
				prod += sqrt  + sqrtInv;
				if (2*j < loops) {
					final double sqrt2IN = sqrt2 * sqrt;
					//					prod += sqrt2IN  + sqrt2Inv *  sqrtInv;
					prod += sqrt2IN  + 1d / sqrt2IN;
				}
			}
			for (long i = 4; i < range; i+=4) {
				final double sqrt = Math.sqrt(i*n);
				prod += sqrt  + 1/sqrt;
			}
		}
		return prod;
	}
	private double sqrtBig(long range) {
		double prod = 1;
		final long n = 1l << 21;
		final double sqrtN =  Math.sqrt(n);
		for (int j = 0; j <loops; j++) {
			for (long i = 0; i < range; i++) {
				final double sqrt = Math.sqrt(i) * sqrtN;
				prod += sqrt  + 1/sqrt ;
			}

		}
		return prod;
	}
	private double multiply(long range) {
		double prod = 1;
		for (int j =1; j <loops; j++) {
			for (int i = 1; i < range; i++) {
				prod += i*i;
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
