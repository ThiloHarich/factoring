package factoring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class Performance {


	@Test
	public void testDivisioin()
	{
		int countWarmup = 0;
		final int nMax = 100000;

		for (int j=0; j<3; j++)
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
		int count = 0;
		long start = System.nanoTime();
		for (long n = 10; n<nMax; n++)
		{
			//			long nDivI = n/3;
			//			for (int i = 3; i <= nDivI; i+=2) {
			for (int i = 3; i*i <= n; i+=2) {
				//				nDivI = n/i;
				if (n% i == 0)
					count++;
			}
		}

		final long time = System.nanoTime() - start;
		System.out.println("division : \t" +time);
		for (int j=0; j<3; j++)
		{
			for (long n = 10; n<nMax; n++)
			{
				double inversIApprox = 1./2;
				long nDivI = n/3 + 1;
				//since we have n/i we can compare with i; we do not have to call i*i <= n or i <= sqrt(n), which is more expensive
				for (int i = 3; i < nDivI; i+= 2) {
					// try to approximate 1/i from 1/(i-1) with the newton method to get around the expensive n/i calculation
					inversIApprox = inversIApprox * (2 - i * inversIApprox);
					// This provides better precision
					final double error = 1 - i * inversIApprox;
					//					inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
					nDivI = Math.round(n*inversIApprox);
					double nDivIMuliplyI = nDivI*i;

					// if the error is to high we can fix it by doing an other iteration or may just in-/decrease nDivIMuliplyIAppox by i
					while (Math.abs(nDivIMuliplyI - n) >= i)
					{
						// do an other iteration
						inversIApprox = inversIApprox * (2 - i * inversIApprox);
						//						inversIApprox = inversIApprox + inversIApprox * (1 - i * inversIApprox);
						nDivI = Math.round(n*inversIApprox);
						nDivIMuliplyI = nDivI*i;
					}
					if (nDivIMuliplyI == n)
						countWarmup++;
				}
			}
		}
		int count2 = 0;
		start = System.nanoTime();

		for (long n = 10; n<nMax; n++)
		{
			double inversIApprox = 1./2;
			long nDivI = n/3 + 1;
			//since we have n/i we can compare with i; we do not have to call i*i <= n or i <= sqrt(n), which is more expensive
			for (int i = 3; i < nDivI; i+=2) {
				// try to approximate 1/i from 1/(i-1) with the newton method to get around the expensive n/i calculation
				inversIApprox = inversIApprox * (2 - i * inversIApprox);
				nDivI = Math.round(n*inversIApprox);
				double nDivIMuliplyI = nDivI*i;

				// if the error is to high we can do an other iteration or may just in-/decrease nDivIMuliplyIAppox by i
				while (Math.abs(nDivIMuliplyI - n) >= i)
				{
					// do an other iteration
					inversIApprox = inversIApprox * (2 - i * inversIApprox);
					nDivI = Math.round(n*inversIApprox);
					nDivIMuliplyI = nDivI*i;
				}
				if (nDivIMuliplyI == n)
					count2++;
			}
		}
		final long time2 = System.nanoTime() - start;
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
		int r = 0;
		// Warmup
		for (int j=0; j<3; j++)
		{
			for (long i = 1; i < length; i++) {
				r += (int) Math.sqrt(i);
			}
		}
		long start = System.nanoTime();
		for (int j=0; j<100000; j++)
		{
			for (long i = 1; i < length; i++) {
				r += (int) Math.sqrt(i);
			}
		}
		long time = System.nanoTime() - start;
		System.out.println("sqrt : \t" +time);
		assertTrue(r != 0);

		// Warmup
		for (int j=0; j<3; j++)
		{
			for (long i = 1; i < length; i++) {
				r += r%i;
			}
		}
		start = System.nanoTime();
		for (int j=0; j<100000; j++)
		{
			for (long i = 1; i < length; i++) {
				r += r%i;
			}
		}
		time = System.nanoTime() - start;
		System.out.println("mult : \t" +time);
		assertTrue(r != 0);


		// Warmup
		for (int j=0; j<3; j++)
		{
			for (double i = 1; i < length; i++) {
				r += i*i;
			}
		}
		start = System.nanoTime();
		for (int j=0; j<100000; j++)
		{
			for (double i = 1; i < length; i++) {
				r += i*i;
			}
		}
		time = System.nanoTime() - start;
		System.out.println("flot : \t" +time);
		assertTrue(r != 0);
	}


}
