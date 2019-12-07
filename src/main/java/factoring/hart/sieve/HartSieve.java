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
package factoring.hart.sieve;

import static org.junit.Assert.assertEquals;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;

import com.google.common.collect.ImmutableSet;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.jml.gcd.Gcd63;
import de.tilman_neumann.util.SortedMultiset;
import factoring.math.PrimeMath;
import factoring.smooth.SmoothNumbers;

/**
 * A factoring algorithm for numbers up to 2^63 numbers. For increasing k it generates numbers
 * x^2 - kn. For each k the primefactors below a limit m of the numbers were calculated by sieving.
 * For smooth numbers over m we try to find a product out of these numbers which is a square.
 * The main difference to the Quadratic sieve is that we do not focus on numbers x^2 - n.
 * The numbers will grow in O(sqrt(n)). 
 * We take into consideration numbers x^2 - k*n. The k will be increased during the process,
 * as long as we do not have enough smooth numbers. 
 * 
 * @author Thilo Harich
 */
public class HartSieve extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartSieve.class);

	// product of 11 primes <= 37 
//	static int [] primes = {2,3,5,7,11,13,17,19,23,29,31,37};
//	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
//	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163};
	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,   277,   281,   283,   293,   307,   311};
//	static int [] primes = {2,     3,     5,     7,    11,    13,    17,    19,    23,    29,    31,    37,    41,    43,
//	   47,    53,    59,    61,    67,    71,    73,    79,    83,    89,    97,   101,   103,   107,
//	  109,   113,   127,   131,   137,   139,   149,   151,   157,   163,   167,   173,   179,   181,
//	  191,   193,   197,   199,   211,   223,   227,   229,   233,   239,   241,   251,   257,   263,
//	  269,   271,   277,   281,   283,   293,   307,   311,   313,   317,   331,   337,   347,   349,
//	  353,   359,   367,   373,   379,   383,   389,   397,   401,   409,   419,   421,   431,   433,
//	  439,   443,   449,   457,   461,   463,   467,   479,   487,   491,   499,   503,   509,   521};
	static double [] primesReciprocal = new double[primes.length];
//	squareRoots[n][j] = i <-> i^2 - n = 0 mod primes[j]
	private static int[][] squareRoots = new int [primes[primes.length-1]][];
	
	/** 
	 * Size of arrays: this is around 4*n^1/3.
	 * 2^21 should work for all number n up to 2^52.
	 */
	private static final int I_MAX = 1<<21;

	/** This constant is used for fast rounding of double values to long. */
	private static final double ROUND_UP_DOUBLE = 0.9999999665;

	private final boolean doTDivFirst;
	private final double[] sqrt1;
//	private final double[] sqrt2;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);
	private final Gcd63 gcdEngine = new Gcd63();



	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
	 * But the smaller possible factors get, it will become slower and slower.
	 */
	public HartSieve(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt1 = new double[I_MAX];
//		sqrt2 = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt1[i] = Math.sqrt(i);
//			if (i%7 != 0) {
//				sqrt2[i] = Math.sqrt(i*K_MULT2);
//			}
		}
		// for all the primes we have an index to a long value which stores the exponent of this
		// prime in the prime factorization. As long as they do not interfere we are fine
		for (int i = primes.length-1; i >= 0; i--) {
//			for (int i = 0; i < primes.length; i++) {
			int p = primes[i];
			primesReciprocal[i] = 1.0 / p;
			squareRoots[i] = new int [p];
			for (int n = 1; n < p; n++) {
				boolean solFound = false;
				for (int j = 0; j <= (p+1) / 2; j++) {
					if ((j*j - n) % p == 0) {
						squareRoots[i][n] = j;
						solFound = true;
					}
				}
				if (!solFound)
					squareRoots[i][n] = -1;						
			}
		}
	}

	@Override
	public String getName() {
		return "Hart_Fast2Mult(" + doTDivFirst + ")";
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
		if (doTDivFirst) {
			// do trial division before the Hart loop
			tdiv.setTestLimit((int) Math.cbrt(N));
			final long factor = tdiv.findSingleFactor(N);
			if (factor > 1) return factor;
		}
		// test for exact squares
		final double sqrtN = Math.sqrt(N);
		final double sqrt4N = sqrtN;
		SmoothNumbers smooth = new SmoothNumbers(primes);
		int smoothCount = 0;
		long a;
//		int primeBound = 40;
		int primeBound = primes.length-1;
		int sieveBound = primes[primeBound];
//		int sieveBound = 36 * 271;
		byte [] factorLength = new byte [2*sieveBound];
		long [] smoothMatrix = new long [primeBound + 20];
		double candidates = 0;
		for (int k=1; smoothCount <= primeBound; k++) {
			long beginXIndex = (long) (sqrt4N * sqrt1[k] + ROUND_UP_DOUBLE);
			double candidates4Level;
//			do {
				candidates4Level = 0;
				System.out.println("k : " + k + "  smooth : " + smoothCount);
				int step = 1;
				// bring beginXIndex to some k to value b = kn+1 mod 2^l
				if ((k&1) == 1){
					beginXIndex |= 1;
					step = 2;
				}
				for (int l = 0; l < sieveBound * step; l++) {
					a = (long) (beginXIndex + l );
					long test = a*a - k * N;
					byte lenght = (byte) (63 - Long.numberOfLeadingZeros(test));
					factorLength[l] = lenght;
				}
	
				for (int j=primeBound; j > 10 ; j--) 
				{
					final int p = primes[j];
					final int nMod = (int) PrimeMath.mod(k * N, p, primesReciprocal[j]);
					final int nMod1 = (int) ((k * N) % p);
					assertEquals(nMod1, nMod);
					int sol = squareRoots[j][nMod];
					if (sol >= 0) {					
						int primeLenght = 32 - Integer.numberOfLeadingZeros(p);
						int beginXIndexMod = (int) PrimeMath.mod(beginXIndex, p, primesReciprocal[j]);
						int solIndex = subtractMod(sol, beginXIndexMod, p);
						if ((k&1) == 1 && (solIndex&1) == 1)
							solIndex += p;
						for (int l = solIndex; l < sieveBound*step; l+= p* step) {
							checkSieving(N, k, beginXIndex, p, sol, l);
							factorLength[l] -= primeLenght;
						}
						solIndex = subtractMod(p - sol, beginXIndexMod, p);
						// make sieve start such that it hits even pos. With factor 2
						if ((k&1) == 1 && (solIndex&1) == 1)
							solIndex += p;
						for (int l = solIndex; l < sieveBound*step; l+= p*step) {
	//							checkSieveing(fourN, k1, beginXIndex, p, sol, l);
							factorLength[l] -= primeLenght;
						}
					}
				}
				for (int l = 0; l < sieveBound * step; l++) {
					if (factorLength[l] < 14) {
						candidates4Level++;
						candidates++;
						a = beginXIndex + l;
						long test = a*a - k * N;
						tdiv.setTestLimit(primes[primeBound]);
						final SortedMultiset<BigInteger> factor = tdiv.factor(BigInteger.valueOf(test));
						System.out.print("a " + a+" kn + 1 " + (k* N) % 8);
						System.out.println(factor);
	//						final boolean isSmooth = smooth.isSmoothBy3Inv(test, PRIMORIAL);
						int[] factors = smooth.factorizeIfSmooth1(test, primes[primeBound]);
						long factorsMod2 = 0;
						for (int j = 0; j < factors.length; j++) {
							factorsMod2 += factors[j] %2 == 1 ? 1l << j : 0;						
						}
						smoothMatrix[smoothCount] = factorsMod2;
						smoothCount += factors[primes.length];
					}
				}
				beginXIndex += sieveBound;
//			}
//			while(candidates4Level >  candidates/k);
		}
		Set<Integer> deleteRowsSet = printMatrix(smoothMatrix);
		deleteRowsSet.stream().map(i -> i + ",").forEach(System.out::print);
		do {
	//		long[] reducedMatrix = IntStream.range(0, smoothMatrix.length).filter(i -> !deleteRowsSet.contains(i)).mapToObj(i -> smoothMatrix[i]).collect(Collectors.toList()).stream().mapToLong(l -> l).toArray();
			long[] reducedMatrix = new long[smoothMatrix.length - deleteRowsSet.size()];
			for (int row=0,insertIndex=0; row < smoothMatrix.length; row++) {
				if (!deleteRowsSet.contains(row))
					reducedMatrix[insertIndex++] = smoothMatrix[row];
			}
			smoothMatrix = reducedMatrix;
			deleteRowsSet = printMatrix(smoothMatrix);
			System.out.println();
			deleteRowsSet.stream().map(i -> i + ",").forEach(System.out::print);
		}while (deleteRowsSet.size() > 0);
		
//		long[] reducedMatrix = new long[smoothMatrix.length - deleteRowsSet.size()];
//		for (int row=0,insertIndex=0; row < smoothMatrix.length; row++) {
//			if (!deleteRowsSet.contains(row))
//				reducedMatrix[insertIndex++] = smoothMatrix[row];
//		}

		for (int col = primeBound; col > 0; col--) {
			int row = smoothCount - col;
			for (; row <= primeBound+1 && (smoothMatrix[row] & (1l << col)) == 0; row++);
			if ((smoothMatrix[row] & (1l << col)) != 0) {
				long exchange = smoothMatrix[row];
				smoothMatrix[row] = smoothMatrix[smoothCount - col];
				smoothMatrix[smoothCount - col] = exchange;
				printMatrix(smoothMatrix);		
				row = smoothCount - col +1;
				for (; row <= primeBound+1; row++) {
					smoothMatrix[row] = smoothMatrix[row] ^ smoothMatrix[smoothCount - col];
				}
				printMatrix(smoothMatrix);		
			}
		}

		return 1;
	}

	protected Set<Integer> printMatrix(long[] smoothMatrix) {
		int [] colCounts = new int [primes.length];
		int [] deleteRows = new int [primes.length];
		for (int j = 0; j < smoothMatrix.length; j++) {
			long row = smoothMatrix[j];
			String line = String.format("%02d", j);
			if (row != 0) {
				for (int l = 0; l < primes.length; l++) {
					final boolean rowContainsPrime = (row & (1l << l)) != 0;
					line += rowContainsPrime ? "o|" : " |";
					if (rowContainsPrime) {
						if (colCounts[l] == 0)
							deleteRows[l] = j+1;
						else
							deleteRows[l] = 0;
						colCounts[l]++;
					}
				}
			}
			if (row != 0)
			System.out.println(line);
		}
		System.out.print(" ");
		for (int i = 0; i < colCounts.length; i++) {
			System.out.print(String.format("%02d", i));			
		}
		System.out.println();
		System.out.print(" ");

		for (int i = 0; i < colCounts.length; i++) {
			if (colCounts[i] == 1)
				System.out.print("XX");		
			else
				System.out.print(String.format("%02d", colCounts[i]));		
		}

		Set<Integer> deleteRowsSet = Arrays.stream(deleteRows).filter(r-> r != 0).mapToObj(i -> Integer.valueOf(i-1)).collect(Collectors.toSet());
		return deleteRowsSet;
	}

	protected int subtractMod(int sol, int beginXIndexMod, final int p) {
		int solIndex = (int) (sol - beginXIndexMod); 
		solIndex = solIndex < 0 ? solIndex + p : solIndex;
		return solIndex;
	}

	protected void checkSieving(final long fourN, int k1, final long beginXIndex, final int p, int sol, int l) {
		long a;
		a = (long) (beginXIndex + l);
		assertEquals(sol, a % p);
		long test = a*a - k1 * fourN;
		assertEquals(0, test % p);
	}
	

	private int findSolution(long kn, int primeIndex) {
		int p = primes[primeIndex];
		if (squareRoots [(int) (kn % p)][primeIndex] != 0)
			return squareRoots[(int) (kn % p)][primeIndex];
		return -1;
	}

	private String factors(long a, long[] primeFacts, int[] composites) {
		String primefactors = "";
		for (int j = 0; j < 12; j++) {
			int exponent = (int) ((primeFacts[0] >> 5*j) & 31);
			if (exponent > 0)
				primefactors += primes[j] + "^" + exponent + " * ";
		}
		final String compString = Arrays.stream(composites).filter(i -> i != 0).mapToObj(Integer::toString).collect(Collectors.joining(" * "));
		return "" + a + " = " +primefactors + compString;
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
		if ((k&1)==0) return x | 1;
		
		final long kNp1 = k*N+1;
		if ((kNp1 & 3) == 0) return x + ((kNp1 - x) & 7);
		
		if ((kNp1 & 7) == 2) {
			final long adjust1 = ( kNp1 - x) & 15;
			final long adjust2 = (-kNp1 - x) & 15;
			return x + (adjust1 < adjust2 ? adjust1 : adjust2);
		}
		
		final long adjust1 = ( kNp1 - x) & 31;
		final long adjust2 = (-kNp1 - x) & 31;
		return x + (adjust1 < adjust2 ? adjust1 : adjust2);
	}

}
