package factoring.math;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * Created by Thilo Harich on 07.03.2018.
 */
public class GdcTest {

	@Test
	public void testPerformance ()
	{
		final Random rand = new Random();

		final int range = 3000;

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

	public void doThilo(final Random rand, final int range) {
		for ( long countA = 0; countA < range; countA++) {
			final long a = Math.abs(rand.nextInt());
			final long bStart = Math.abs(rand.nextInt());
			for (long bCount = 1; bCount < range; bCount++) {
				final long gdc = PrimeMath.gcdCached(a, bStart+bCount);
			}
		}
	}

	public void doCompare(final Random rand, final int range) {
		final Gcd63 gcd2 = new Gcd63();
		for ( long countA = 0; countA < range; countA++) {
			final long a = Math.abs(rand.nextInt());
			final long bStart = Math.abs(rand.nextInt());
			for (long bCount = 1; bCount < range; bCount++) {
				//				gcd2.gcd(a, bStart+bCount);
				PrimeMath.gcdDouble((double)a, bStart+bCount);
			}
		}
	}

	@Test
	public void testCorrect ()
	{
		final Random rand = new Random();
		final long range = 50;
		final Gcd63 gcd63 = new Gcd63();
		long a,b;


		for ( long countA = 0; countA < range; countA++) {
			a = Math.abs(rand.nextInt());
			//			a = 2;
			for (  long countB = 0; countB < range; countB++) {
				b = Math.abs(rand.nextInt());
				//				b = 3;
				//            long gcd1 = PrimeMath.gcd(begin + offset, other);
				final long gcd2 = gcd63.gcd(a, b);
				//				final long gcd3 = PrimeMath.gcdByMod((double)a, (double)b);
				final long gcd3 = PrimeMath.gcdCached(a, b);
				//            Assert.assertEquals (gcd1, gcd2);
				Assert.assertEquals (a + ", " + b, gcd2, gcd3);
			}
		}
	}

}
