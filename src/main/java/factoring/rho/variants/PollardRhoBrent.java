package factoring.rho.variants;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

/**
 * Uses the brent variant of the loop detection algorithm.
 *
 * @author thiloharich
 *
 */
public class PollardRhoBrent implements FactorizationOfLongs, FactorFinder {

	Random rnd = new Random();

	@Override
	public boolean findsPrimesOnly() {
		return false;
	}

	@Override
	public boolean returnsCompositeOnly() {
		return true;
	}


	@Override
	/**
	 * Here we do the brent style cycle detection.
	 * So we look for loops of length 2^k.
	 * We compare the first
	 */
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (PrimeMath.isPrime32((int) n))
			//			if (PrimeMath.isPrime41Bit(n))
			return n;
		long xFixed = 2;
		int fails = 0;
		long factor;

		do {
			factor = 1;
			long size = 1;
			final long c = Math.abs(rnd.nextInt() % n);
			long x = rnd.nextInt() % n;

			do {
				for (int count = 1; (count <= size) && (factor <= 1); count++) {
					x = g(n, x, c);
					factor = PrimeMath.gcd(Math.abs(x - xFixed), n);
				}
				size = 2 * size;
				xFixed = x;
			}while (factor == 1);
			fails++;
		}while (factor == n && fails < 10);
		return factor;
	}

	long g (long n, long x, long c) {
		return (x * x  + c) % n;
	}

}
