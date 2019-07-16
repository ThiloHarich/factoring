/*
 * java-math-library is a Java library focused on number theory, but not necessarily limited to it. It is based on the PSIQS 4.0 factoring project.
 * Copyright (C) 2018 Tilman Neumann (www.tilman-neumann.de)
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program;
 * if not, see <http://www.gnu.org/licenses/>.
 */
package factoring.hart;

import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;
import factoring.math.PrimeMath;

/**
 * Pretty simple yet fast variant of Hart's one line factorizer.
 *
 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
 * But the smaller possible factors get, it will become slower and slower.
 *
 * For any kind of test numbers except very hard semiprimes, Hart_TDiv_Race will be faster.
 *
 * @authors Thilo Harich & Tilman Neumann
 */
public class HartTraining extends FactorAlgorithm {


	//	private static final String FACTOR_SEQUENCE_FILE = "/tmp/factorSequence.txt";
	private static final String FACTOR_SEQUENCE_FILE = "./factorSequence.txt";

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	//	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 3*5*7; // 105
	//	private static final int K_MULT = 3*3*5; // 45
	private static final int K_MULT = 1; // 315

	/** Size of arrays */
	private static final int I_MAX = 1<<17;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();

	public static int [] ks = new int[I_MAX];
	public static int [] hits;

	boolean train;

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public HartTraining(boolean doTDivFirst, boolean train) {
		this.doTDivFirst = doTDivFirst;
		this.train = train;
		// Precompute sqrts for all k < I_MAX
		readKSequence(train);
	}

	@Override
	public String getName() {
		return "Hart_FastT(" + doTDivFirst + ")";
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		return BigInteger.valueOf(findSingleFactor(N.longValue()));
	}

	/**
	 * Find a factor of long N.
	 * @param N
	 * @return factor of N
	 */
	public long findSingleFactor(long N) {
		// do trial division before the Hart loop ?
		long factor;
		if (doTDivFirst) {
			tdiv.setTestLimit((int) Math.cbrt(N));
			if ((factor = tdiv.findSingleFactor(N))>1) return factor;
		} else {
			// at least do trial division by multiplier factors
			if ((N%3)==0) return 3;
			if ((N%5)==0) return 5;
			if ((N%7)==0) return 7;
		}

		final long fourN = N<<2;
		final double sqrt4N = Math.sqrt(fourN);
		long a, b, test, gcd;
		try {
			for (int i=1; ;i++) {
				final int k = ks[i];
				// odd k -> adjust a mod 8
				a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				a = adjustA(N, a, k);
				test = a*a - k * fourN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
						if (train)
							hits[i]++;
						//						else
						return gcd;
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// this may happen if this implementation is tested with doTDivFirst==false and N having
			// very small factors, or if N is too big
		}
		return 0;
	}

	/**
	 * Increases x to return the next possible solution for x for x^2 - 4kn = b^2.
	 * Due to performance reasons we give back solutions for this equations modulo a
	 * power of 2, since we can determine the solutions just by additions and binary
	 * operations.
	 *
	 * if k is even x must be odd.
	 * if k*n == 3 mod 4 -> x = k*n+1 mod 8
	 * if k*n == 1 mod 8 -> x = k*n+1 mod 16 or -k*n+1 mod 16
	 * if k*n == 5 mod 8 -> x = k*n+1 mod 32 or -k*n+1 mod 32
	 *
	 * @param N
	 * @param x
	 * @param k
	 * @return
	 */
	private long adjustA(long N, long x, int k) {
		if ((k&1)==0)
			return x | 1;
		final long kNp1 = k*N+1;
		if ((kNp1 & 3) == 0)
		{
			return x + ((kNp1 - x) & 7);
		}
		else if ((kNp1 & 7) == 2) {
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			final long diff = adjust1<adjust2 ? adjust1 : adjust2;
			return x + diff;
		}
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1<adjust2 ? adjust1 : adjust2);
	}


	public void readKSequence(boolean train) {
		final Path path = Paths.get(FACTOR_SEQUENCE_FILE);
		sqrt = new double[I_MAX];
		ks = new int[I_MAX];
		if (!train && Files.exists(path)) {
			List<String> lines = null;
			final boolean [] added = new boolean[I_MAX];
			int i = 0;
			try {
				lines = Files.readAllLines(path);
				final String line = lines.get(0);
				final String[] numbers = line.split(",");
				for (; i < numbers.length; i++) {
					final String number = numbers[i];
					if (!number.isEmpty()) {
						final int k = Integer.parseInt(number);
						ks[i] = k * K_MULT;
						sqrt[i] = Math.sqrt(k*K_MULT);
						added[k] = true;
					}
				}
			} catch (final IOException e) {
				e.printStackTrace();
			}
			// fill up the k Array with k's that were not hit by the training
			for (int k = 1; k < added.length; k++) {
				if (!added[k]) {
					ks[i] = k * K_MULT;
					sqrt[i++] = Math.sqrt(k*K_MULT);
				}
			}
		}else {
			hits = new int[I_MAX];
			for (int k = 0; k < ks.length; k++) {
				ks[k] = k * K_MULT;
				sqrt[k] = Math.sqrt(k*K_MULT);
			}
		}
	}

	/**
	 * Test.
	 * @param args ignored
	 * @throws IOException
	 */
	public static void main(String[] args) {

		final HartTraining hart = new HartTraining(false, true);
		final int numPrimes = 200_001;
		final long[] semiprimes = makeSemiPrimesListReal(50, numPrimes);
		for (final long l : semiprimes) {
			hart.findSingleFactor(BigInteger.valueOf(l));
		}
		final List<Pair> pairs = new ArrayList<>();
		for (int i = 0; i < HartTraining.hits.length; i++) {
			final int k = HartTraining.hits[i];
			//			if (k > HartTraining.hits[1] / 10)
			{
				//			if (k > 0) {
				final Pair pair = hart.new Pair(i, k);
				pairs.add (pair);
				if (i < 100)
					System.out.println(i + "\t" + pair + "," + hart.factor(BigInteger.valueOf(pair.k)));
			}
		}
		System.out.println("found " + pairs.size() + " factors");
		final Pair [] pairArr = new Pair[pairs.size()];
		pairs.toArray(pairArr);
		Arrays.sort(pairArr);
		int pos = 1;
		for (final Pair pair : pairArr) {
			//			if (pos < 100) {
			if (pos < 1000) {
				//				if (pair.k > 20000 && pair.hits > 0) {
				System.out.println(pos++ + "\t" + pair + "," + hart.factor(BigInteger.valueOf(pair.k)));
			}
		}
		String factorSequence = "";
		for (final Pair pair : pairArr) {
			factorSequence += pair.k + ",";
		}
		System.out.println("wrote " + pairArr.length + " factors ");
		final Path path = Paths.get(FACTOR_SEQUENCE_FILE);
		final byte[] strToBytes = factorSequence.getBytes();

		try {
			Files.write(path, strToBytes);
		} catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	public static long[] makeSemiPrimesListReal(int bits, int numPrimes) {
		final long[] semiPrimes = new long[numPrimes];
		final int limit = (int) Math.pow(2,bits / 3.0);
		final TDiv63Inverse factorizer = new TDiv63Inverse(limit);

		long candidate = (1l << bits);
		//		candidate -= (candidate % 315) - 4;
		//		assertEquals(candidate % 315, 4);
		int j = 0;
		for (int i=0; i< numPrimes; candidate++)
		{
			final BigInteger findSingleFactor = factorizer.findSingleFactor(BigInteger.valueOf(candidate));
			if (findSingleFactor == BigInteger.ONE && !PrimeMath.isPrime(candidate)) {
				semiPrimes[j++] = candidate;
				i++;
			}
		}

		return semiPrimes;
	}
	public class Pair implements Comparable<Pair>{
		public Pair(int k, int hits) {
			this.k = k;
			this.hits = hits;
		}
		int k;
		int hits;
		@Override
		public int compareTo(Pair o) {
			return o.hits - hits;
		}
		@Override
		public String toString() {
			return "Pair [k=" + k + ", hits=" + hits+ "]";
		}
	}

}
