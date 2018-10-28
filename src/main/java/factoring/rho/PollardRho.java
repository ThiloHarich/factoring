package factoring.rho;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

public class PollardRho implements FactorizationOfLongs, FactorFinder {

	Random rnd = new Random();

	@Override
	public boolean findsPrimes() {
		return false;
	}

	@Override
	public boolean findsOneFactor() {
		return true;
	}


	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		long gcd = 1;
		do {
			long x = rnd.nextInt() % n;
			final long c = rnd.nextInt() % n;
			long y = x;

			do {
				x = g(n, x, c);
				y = g(n, y, c);
				y = g(n, y, c);

				gcd = PrimeMath.gcd(Math.abs(y-x), n);
			}while (gcd == 1);
		}
		while (gcd == n);
		return gcd;
	}

	long g (long n, long x, long c) {
		return (x * x  + c) % n;
	}

}
