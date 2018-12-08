package factoring.rho;

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
 * It only works for numbers up to 4503599627370496, which is 52 Bits the precision of a double.
 * It is slower (by a factor of ~1.8) then the {@link PollardRhoBrentDoubleUnroll}
 * implementation, which is only working for up to 26 Bits.
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
public class PollardRhoBrentDouble52 implements FactorizationOfLongs, FactorFinder {

	Random rnd = new Random();
	private int maximalFactor = Integer.MAX_VALUE;
	private int maxLoops = Integer.MAX_VALUE;


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
	public long findFactors(final long n, Collection<Long> primeFactors) {
		// TODO we need a check for at lest 52 bit.
		if (PrimeMath.isPrime41Bit(n))
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
		int loops = 0;


		do {
			factor = 1;
			int numbersChecked = 0;
			final long c = Math.abs(rnd.nextInt() % n);
			//			c = (c + 2) % n;
			long x = Math.abs(rnd.nextInt()) % n;
			//			x = 0;

			do {
				xFixed = x;
				do {
					long xProd = 1;
					long[] xProds = new long [2];
					xProds[1] = 1;
					// iterate as long as we do not have reached the maximal number of
					// factors needed for the gcd
					while (numbersChecked < limit) {
						x = g(x, c, n, nInv);
						// we might exit early if xDiff == 0
						final long xDiff = Math.abs(x - xFixed);
						xProds = BigDouble.multiply(xProd, xDiff);
						xProd  = BigDouble.mod(xProds, n, nInv);
						numbersChecked++;
						loops++;
					}
					factor = PrimeMath.gcd(xProd, n);
					xProd = 1;
					limit += productSize;
				}
				// check if we do not have already exceeded the search interval
				// xDiff = 0 -> factor == n
				while (factor == 1 && numbersChecked < intervalLimit);

				intervalLimit = 2 * intervalLimit;
				limit = Math.min(productLimit, intervalLimit);
				numbersChecked = 0;
				// exit if we have reached the maximal factor to search for
				if (loops > maxLoops)
					return n;
			}while (factor == 1);
			fails++;
			if (factor == n)
			{
				// look for the first gcd different from 1
				x = xFixed;
				do {
					x = g(x, c, n, nInv);
					final long xDiff = Math.abs(x - xFixed);
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
	 *
	 * @param x
	 * @param c
	 * @param n
	 * @param nInv
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

	@Override
	public void setMaxFactor(int maximalFactor) {
		this.maximalFactor = maximalFactor;
		maxLoops = (int) (6 * Math.sqrt(maximalFactor));
	}



}
