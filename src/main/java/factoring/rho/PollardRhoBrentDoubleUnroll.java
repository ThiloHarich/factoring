package factoring.rho;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

/**
 * This is a variant of the Brent like Pollard Rho integer factoring algorithm.
 * Instead of computing the modulus of long values by "x % n" which is done
 * internally by a division, we multiply the long value x by the double value of 1/n.
 * This speeds up the hole algorithm by a factor of 3.
 * It simultaneously checks two values g(x) and g(g(x)) in one loop,
 * since this gives something like 10% performance gain.
 * It only works for numbers up to 67108864, which is 26 Bits (half of the precision of double).
 *
 * From a performance view the brent variant has the advantage, that it replaces the
 * very slow GCD operation by a fast product operation. So the bottleneck
 * is the modulus operation. This is replaced here by a multiplication with double.
 * There are other variations to calculate the mod operation fast :
 * https://en.wikipedia.org/wiki/Barrett_reduction
 * and https://en.wikipedia.org/wiki/Montgomery_modular_multiplication.
 * Even if the Barret reduction is faster when applied to single multiplications,
 * they do not seem to be faster when applied to the recursion x = x^2 + c mod n.
 *
 * @author thiloharich
 *
 */
public class PollardRhoBrentDoubleUnroll implements FactorizationOfLongs, FactorFinder {

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
		final double nInv = 1.0d/ n;
		//		long c = 1;
		int fails = 0;
		long factor;
		long xFixed;
		final int productSize = 16;
		int intervalLimit = 2;
		final int productLimit = productSize;
		int limit = Math.min(productLimit, intervalLimit);


		do {
			factor = 1;
			int numbersChecked = 0;
			final long c = Math.abs(rnd.nextInt() % n);
			//			c = (c + 2) % n;
			long x2 = Math.abs(rnd.nextInt()) % n;
			//			x = 0;

			long x1;
			do {
				xFixed = x2;
				do {
					long xProd = 1;
					// iterate as long as we do not have reached the maximal number of
					// factors needed for the gcd
					while (numbersChecked < limit) {
						x1 = PrimeMath.mod(x2 * x2  + c,  n, nInv);
						x2 = PrimeMath.mod(x1 * x1  + c,  n, nInv);
						// we might exit early if xDiff == 0
						final long xDiff1 = Math.abs(x1 - xFixed);
						final long xDiff2 = Math.abs(x2 - xFixed);
						xProd = PrimeMath.mod(xProd * xDiff1, n, nInv);
						xProd = PrimeMath.mod(xProd * xDiff2, n, nInv);
						numbersChecked+=2;
					}
					factor = PrimeMath.gcd(xProd, n);
					limit += productSize;
				}
				// check if we do not have already exceeded the search interval
				// xDiff = 0 -> factor == n
				while (factor == 1 && numbersChecked < intervalLimit);

				intervalLimit = 2 * intervalLimit;
				limit = Math.min(productLimit, intervalLimit);
				numbersChecked = 0;
			}while (factor == 1);
			fails++;
			if (factor == n)
			{
				// look for the first gcd different from 1
				x1 = xFixed;
				do {
					x1 = g(x1, c, n, nInv);
					final long xDiff = Math.abs(x1 - xFixed);
					factor = PrimeMath.gcd(xDiff, n);
				}while(factor == 1);
			}
		}while ((factor == 1 || factor == n) && fails < 10);
		//		System.out.println("brent " + loops);
		return factor;
	}


	/**
	 * This function should be
	 * a) fast
	 * b) pseudorandom
	 *
	 * To be fast we should avoid the % and replace by iterated subtractions.
	 * To be pseudorandom
	 *
	 *
	 *
	 * @param n
	 * @param x
	 * @param c
	 * @return
	 */
	long g (long x, long c, long n, double nInv) {
		return PrimeMath.mod(x * x  + c,  n, nInv);
	}

}
