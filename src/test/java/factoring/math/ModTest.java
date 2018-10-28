package factoring.math;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import de.tilman_neumann.jml.gcd.Gcd63;

/**
 * Created by Thilo Harich on 07.03.2018.
 */
public class ModTest {

	@Test
	public void testPerformance ()
	{
		final Random rand = new Random();

		final int range = 30000;

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
			final long a = Math.abs(rand.nextLong());
			final long bStart = Math.abs(rand.nextLong());
			for (long bCount = 1; bCount < range; bCount++) {
				final long mod = PrimeMath.modByShift(a, bStart+bCount);
				//				System.out.println(a + ",");
			}
		}
	}

	public void doCompare(final Random rand, final int range) {
		for ( long countA = 0; countA < range; countA++) {
			final long a = Math.abs(rand.nextLong());
			final long bStart = Math.abs(rand.nextLong());
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
		final Gcd63 gcd63 = new Gcd63();
		long a,b;


		for ( long countA = 0; countA < range; countA++) {
			a = Math.abs(rand.nextLong());
			//			a = 52;
			for (  long countB = 0; countB < range; countB++) {
				b = Math.abs(rand.nextLong());
				//				b = 55;
				//            long gcd1 = PrimeMath.gcd(begin + offset, other);
				final long gcd2 = PrimeMath.binaryGcd(a, b);
				final long gcd3 = PrimeMath.gcd(a, b);
				//            Assert.assertEquals (gcd1, gcd2);
				Assert.assertEquals (a + ", " + b, gcd2, gcd3);
			}
		}
	}

}
