package factoring.rho.variants;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

/**
 * Classical Pollard Rho implementation.
 * Can only factor numbers up to 31 bits, since it uses long values
 * @author thiloharich
 *
 */
public class PollardRho implements FactorizationOfLongs, FactorFinder {

	Random rnd = new Random();

	private int maxFactorSqrt = Integer.MAX_VALUE;



	public PollardRho() {
	}

	public PollardRho(int maxFactor) {
		setMaxFactor(maxFactor);
	}

	@Override
	public boolean findsPrimesOnly() {
		return false;
	}

	@Override
	public boolean returnsCompositeOnly() {
		return true;
	}
	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (PrimeMath.isPrime32((int) n))
			//		if (PrimeMath.isPrime41Bit(n))
			return n;
		long factor = 1;
		int fails = 0;
		long c = 1;
		do {
			long x = rnd.nextInt() % n;
			x = 0;
			//			final long c = rnd.nextInt() % n;
			c = (c + 2) % n;
			long y = x;

			do {
				x = g(n, x, c);
				y = g(n, y, c);
				y = g(n, y, c);

				factor = PrimeMath.gcd(Math.abs(y-x), n);
			}while (factor == 1);
			fails++;
		}
		while (factor == n && fails < 10);
		return factor;
	}

	long g (long n, long x, long c) {
		return (x * x  + c) % n;
	}

	@Override
	public void setMaxFactor(int maxFactor) {
		maxFactorSqrt = (int) Math.sqrt(maxFactor);
	}

}
