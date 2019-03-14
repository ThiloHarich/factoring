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
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.log4j.Logger;

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


	private static final Logger LOG = Logger.getLogger(HartTraining.class);

	/**
	 * We only test k-values that are multiples of this constant.
	 * Best values for performance are 315, 45, 105, 15 and 3, in that order.
	 */
	private static final int K_MULT = 3*3*5*7; // 315
	//	private static final int K_MULT = 1; // 315

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

	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 */
	public HartTraining(boolean doTDivFirst, boolean train) {
		this.doTDivFirst = doTDivFirst;
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
			for (int i=0; ;i++) {
				final int k = ks[i];
				// odd k -> adjust a mod 8
				a = (long) (sqrt4N * sqrt[i] + ROUND_UP_DOUBLE);
				a = adjustA(N, a, k, k);
				test = a*a - k * fourN;
				b = (long) Math.sqrt(test);
				if (b*b == test) {
					if ((gcd = gcdEngine.gcd(a+b, N))>1 && gcd<N) {
						if (hits != null)
							hits[i]++;
						return gcd;
					}
				}
			}
		} catch (final ArrayIndexOutOfBoundsException e) {
			// this may happen if this implementation is tested with doTDivFirst==false and N having
			// very small factors, or if N is too big
			return 0;
		}
	}

	private long adjustA(long N, long a, int k, int i) {
		if ((i & 1) == 0)
			a |= 1;
		else {
			final long kPlusN = k + N;
			if ((kPlusN & 3) == 0) {
				a += ((kPlusN - a) & 7);
			} else {
				final long adjust1 = (kPlusN - a) & 15;
				final long adjust2 = (-kPlusN - a) & 15;
				a += adjust1<adjust2 ? adjust1 : adjust2;
			}
		}
		return a;
	}

	public void readKSequence(boolean train) {
		final Path path = Paths.get("/tmp/factorSequence.txt");
		sqrt = new double[I_MAX];
		ks = new int[I_MAX];
		if (!train && Files.exists(path)) {
			List<String> lines = null;
			try {
				lines = Files.readAllLines(path);
				final String line = lines.get(0);
				final String[] numbers = line.split(",");
				for (int i = 0; i < numbers.length; i++) {
					final String number = numbers[i];
					if (!number.isEmpty()) {
						final int k = Integer.parseInt(number);
						ks[i] = k * K_MULT;
						sqrt[i] = Math.sqrt(k*K_MULT);
					}
				}
			} catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
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
		final int numPrimes = 1_000_000;
		final long[] semiprimes = makeSemiPrimesListReal(40, numPrimes);
		for (final long l : semiprimes) {
			hart.findSingleFactor(BigInteger.valueOf(l));
		}
		final Pair[] pairs = new Pair[HartTraining.hits.length];
		for (int i = 0; i < HartTraining.hits.length; i++) {
			final int k = HartTraining.hits[i];
			final Pair pair = hart.new Pair(i, k);
			pairs[i] = pair;
		}
		Arrays.sort(pairs);
		String factorSequence = "";
		for (final Pair pair : pairs) {
			if (pair.hits >0) {
				factorSequence += pair.k + ",";
				//				System.out.print(pair.k + ",");
			}
		}
		final Path path = Paths.get("/tmp/factorSequence.txt");
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

		//		final long offset = 0;
		final Random rnd = new Random();
		final long offset = Math.abs(rnd.nextInt());
		long candidate = (1l << bits) + offset;
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
			return "Pair [k=" + k + ", hits=" + hits + "]";
		}
	}

}
