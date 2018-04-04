package factoring.trial;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class Performance {


	@Test
	public void testDivisioin()
	{
		final int nMax = 100000;
		int countWarmup = 0;
		int warmup = 5;

		for (int j=0; j< warmup; j++)
		{
			for (long n = 10; n<nMax; n++)
			{
				double inversIApprox = 1./2;
				int nDivI;
				//since we have n/i we can compare with i; we do not have to call i*i <= n or i <= sqrt(n), which is more expensive
				for (int i = 3; i*i < n; i+= 2) {
					// try to approximate 1/i from 1/(i-1) with the newton method to get around the expensive n/i calculation
					inversIApprox = inversIApprox * (2 - i * inversIApprox);
					// This provides better precision
//					final double error = 1 - i * inversIApprox;
					//					inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
					nDivI = (int) Math.round(n*inversIApprox);
					int nDivIMultiplyI = nDivI*i;

					// if the error is to high we can fix it by doing an other iteration or may just in-/decrease nDivIMuliplyIAppox by i
					while (Math.abs(nDivIMultiplyI-n) >= i )
					{
//						nDivIMultiplyI += i;  // n += i;
						inversIApprox = 1/(double)i;
						// do an other iteration
//						inversIApprox = inversIApprox * (2 - i * inversIApprox);
//						inversIApprox = 1/i;
						//						inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
						nDivI = (int) n/i;
						nDivIMultiplyI = nDivI*i;
//						nDivIMultiplyI = ((int)(n*inversIApprox))*i;
					}
					if (nDivIMultiplyI == n)
						countWarmup++;
				}
			}
		}
		assertTrue (countWarmup > 0);
		int count2 = 0;
		long start = System.nanoTime();

		for (int j=0; j< warmup; j++) {
			for (long n = 1000; n < nMax; n++) {
				double inversIApprox = 1. / 2;
				int nDivI;
				//since we have n/i we can compare with i; we do not have to call i*i <= n or i <= sqrt(n), which is more expensive
				for (int i = 3; i * i < n; i += 2) {
					// try to approximate 1/i from 1/(i-1) with the newton method to get around the expensive n/i calculation
					inversIApprox = inversIApprox * (2 - i * inversIApprox);
					// This provides better precision
//					final double error = 1 - i * inversIApprox;
					//					inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
					nDivI = (int) Math.round(n*inversIApprox);
//					nDivI = (int) n*inversIApprox;
					int nDivIMultiplyI = nDivI*i;

					// if the error is to high we can fix it by doing an other iteration or may just in-/decrease nDivIMuliplyIAppox by i
					while (Math.abs(nDivIMultiplyI-n) >= i )
//					while (nDivIMultiplyI < n )
					{
//						nDivIMultiplyI += i;  // n += i;
						inversIApprox = 1/(double)i;
						// do an other iteration
//						inversIApprox = inversIApprox * (2 - i * inversIApprox);
//						inversIApprox = 1/i;
						//						inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
						nDivI = (int) n/i;
						nDivIMultiplyI = nDivI*i;
//						nDivIMultiplyI = ((int)(n*inversIApprox))*i;
					}
					if (nDivIMultiplyI == n)
						count2++;
				}
			}
		}
		final long time2 = System.nanoTime() - start;

		for (int j=0; j<warmup; j++)
		{
			for (long n = 10; n<nMax; n++)
			{
				//				long nDivI = n/3;
				//				for (int i = 3; i <= nDivI; i+=2) {
				for (int i = 3; i*i <= n; i+=2) {
					//					nDivI = n/i;
					if (n%i == 0)
						countWarmup++;
				}
			}
		}
		assertTrue (countWarmup > 0);
		int count = 0;
		start = System.nanoTime();
		for (int j=0; j< warmup; j++) {
			for (long n = 10; n < nMax; n++) {
				//			long nDivI = n/3;
				//			for (int i = 3; i <= nDivI; i+=2) {
				for (int i = 3; i * i <= n; i += 2) {
					//				nDivI = n/i;
					if (n % i == 0)
						count++;
				}
			}
		}

		final long time = System.nanoTime() - start;
		System.out.println("division : \t" +time);
		System.out.println("newton : \t" +time2);
		System.out.println("Speedup : \t" + (time + 0.0)/time2);

		assertEquals(count, count2);
	}



	/**
	 *
	 */
	@Test
	public void testSqrt()
	{
		final long length = 1000;
		long n = 78497238l;
		double r = 0;
		int warmup = 10000;
		int loop = 100000;
		// Warmup
		for (int j=0; j<warmup; j++)
		{
			for (long i = 1; i < length; i++) {
				r += n*(1/(double) i);
			}
		}
		long start = System.nanoTime();
		for (int j=0; j<loop; j++)
		{
			for (long i = 1; i < length; i++) {
				r += n*(1/(double) i);
			}
		}
		long time = System.nanoTime() - start;
		System.out.println("1/ : \t" +time);
		assertTrue(r != 0);

		// Warmup
		for (int j=0; j<warmup; j++)
		{
			for (int i = 10; i < length; i++) {
				r += n/i;
			}
		}
		start = System.nanoTime();
		for (int j=0; j<loop; j++)
		{
			for (int i = 10; i < length; i++) {
				r += n/i;
			}
		}
		time = System.nanoTime() - start;
		System.out.println("n/i : \t" +time);
		assertTrue(r != 0);


		// Warmup
		for (int j=0; j<warmup; j++)
		{
			double inversIApprox = 1./2.;
			for (int i = 1; i < length; i++) {
				inversIApprox = inversIApprox * (2 - i * inversIApprox);
				// This provides better precision
				//					inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
				double nDivI = Math.round(n*inversIApprox);
				r+= nDivI*i;
			}
		}
		start = System.nanoTime();
		for (int j=0; j<loop; j++)
		{
			double inversIApprox = 1./9.;
			for (int i = 10; i < length; i++) {
				inversIApprox = inversIApprox * (2 - i * inversIApprox);
				// This provides better precision
//									inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
				int nDivI = (int) (n*inversIApprox);
//				r+= inversIApprox*i;
				int nDivIMulI = nDivI*i;
				while(nDivIMulI - n <= -i)
					nDivIMulI += i;
				while(nDivIMulI - n >= i)
					nDivIMulI -= i;
				r+= nDivI*i ;

//				r+= nDivI*i + i;
			}
		}
		time = System.nanoTime() - start;
		System.out.println("*d : \t" +time);
		assertTrue(r != 0);
	}


}
