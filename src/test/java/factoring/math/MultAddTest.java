package factoring.math;

import org.junit.Test;

/**
 * Created by Thilo Harich on 08.01.2018.
 */
public class MultAddTest {

	@Test
	public void testPerformance ()
	{
		final int bits = 21;
		final long begin = (1L << bits) + 3534;
		final long range = 10000000L;
		final int n3 = 1 << (bits / 3);

		long start = System.currentTimeMillis();
		mult(begin, range, n3);
		long end  = System.currentTimeMillis();
		System.out.println(" time mult  : " + (end- start));

		start = System.currentTimeMillis();
		add(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time add   : " + (end- start));

		start = System.currentTimeMillis();
		mult(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time mult  : " + (end- start));

		start = System.currentTimeMillis();
		add(begin, range, n3);
		end  = System.currentTimeMillis();
		System.out.println(" time add   : " + (end- start));

	}

	public void mult(long begin, long range, int n3) {
		boolean found = false;
		//        long  n = begin;
		for (int k = 1; k< n3; k++)
		{
			for(long x = begin; x < begin + range; x++)
			{
				final long x2 = x*x;
				final long right = x2 - begin;
				if (PrimeMath.isSquare(right)) {
					found = ! found;
				}
			}
		}
	}
	public void add(long begin, long range, int n3) {
		boolean found = false;
		//        long  n = begin;
		for (int k = 1; k< n3; k++){
			int x = (int) begin;
			int right = (int) (x*x - begin);
			for(; x < begin + range; x++)
			{
				right = right + 2*x + 1;
				if (PrimeMath.isSquare(right)) {
					found = ! found;
				}
			}
		}
	}
	//	public void addMod(long begin, long range, int n3) {
	//		boolean found = false;
	//		final int mod = 64;
	//		//        long  n = begin;
	//		for (int k = 1; k< n3; k++){
	//			int x = (int) begin;
	//			byte right = (byte) (x*x - begin);
	//			for(; x < begin + range; x++)
	//			{
	//				right = (byte) ((right + 2*xByte + 1) % mod);
	//				if (PrimeMath.isSquare(right)) {
	//					found = ! found;
	//				}
	//			}
	//		}
	//	}

}
