package factoring.rho.variants;

import java.util.Collection;
import java.util.Random;

import factoring.FactorFinder;
import factoring.FactorizationOfLongs;
import factoring.math.BigDouble;
import factoring.math.PrimeMath;

/**
 * This is a variant of the Brent like Pollard Rho integer factoring algorithm.
 * Instead of computing the modulus of long values by "x % n" which is done
 * internally by a division, we multiply the long value x by the double value of 1/n.
 * This speeds up the hole algorithm by a factor of 3.
 * It simultaneously checks two values g(x) and g(g(x)) in one loop,
 * since this gives something like 10% performance gain.
 * It only works for numbers up to 4503599627370496, which is 52 Bits the precision of a double.
 *
 * @author thiloharich
 *
 */
public class PollardRhoBrentDoubleUnroll52 implements FactorizationOfLongs, FactorFinder {

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
		if (PrimeMath.isPrime41Bit(n))
			//			if (PrimeMath.isPrime41Bit(n))
			return n;
		final double nInv = 1.0d/ n;
		//		final double [] nInvD = BigDouble.invert2Double(n);
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
					long[] xProds = new long [2];
					xProds[1] = 1;
					// iterate as long as we do not have reached the maximal number of
					// factors needed for the gcd
					while (numbersChecked < limit) {
						x1 = g (x2, c,  n, nInv);
						x2 = g (x1, c,  n, nInv);
						// we might exit early if xDiff == 0
						final long xDiff1 = Math.abs(x1 - xFixed);
						final long xDiff2 = Math.abs(x2 - xFixed);
						// using a represe double/long
						xProds = BigDouble.multiply(xProd, xDiff1);
						xProd  = BigDouble.mod      (xProds, n, nInv);
						xProds = BigDouble.multiply(xProd, xDiff2);
						xProd  = BigDouble.mod      (xProds, n, nInv);
						//						xProd = PrimeMath.mod(xProd * xDiff1, n, nInv);
						//						xProd = PrimeMath.mod(xProd * xDiff2, n, nInv);
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
	 * we multiply lowerBits (x) * higherBits (x) such that the result is at most 52 Bits,
	 * -> lowerBits, higherBits have a length of 26
	 *
	 * @param n
	 * @param x
	 * @param c
	 * @return
	 */
	long g (long x, long c, long n, double nInv) {
		final long[] gs = BigDouble.multiply(x, x);
		long g = BigDouble.mod(gs, n, nInv);
		g += c;
		if (g > n)
			g -= n;
		return g;
	}

}
