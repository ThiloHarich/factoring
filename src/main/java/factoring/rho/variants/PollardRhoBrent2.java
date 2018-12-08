package factoring.rho.variants;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

public class PollardRhoBrent2 implements FactorizationOfLongs, FactorFinder {

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
		int fails = 0;
		long factor;
		int count = 0;
		long size = 2;


		do {
			factor = 1;
			//			final long c = Math.abs(rnd.nextInt() % n);
			final long c = 1;
			long x = rnd.nextInt() % n;
			long xFixed = x;

			do {
				x = g(n, x, c);
				factor = PrimeMath.gcd(Math.abs(x - xFixed), n);
				count++;
				if (count == size) {
					size = 2 * size;
					xFixed = x;
				}
			}while (factor == 1);
			fails++;
		}while (factor == n && fails < 10);
		//		System.out.println("brent " + loops);
		return factor;
	}

	long g (long n, long x, long c) {
		return (x * x  + c) % n;
	}

}
