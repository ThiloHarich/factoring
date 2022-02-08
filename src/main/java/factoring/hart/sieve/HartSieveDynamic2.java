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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
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
public class HartSieveDynamic2 extends FactorAlgorithm {
	private static final Logger LOG = Logger.getLogger(HartSieveDynamic2.class);

	// product of 11 primes <= 37 
//	static int [] primes = {2,3,5,7,11,13,17,19,23,29,31,37};
//	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
//	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163};
//	static int [] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,   277,   281,   283,   293,   307,   311};
	static int [] primes = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
			103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
			199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
			313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
			433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
			563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
			673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
			811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
			941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,
			1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,
			1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,
			1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,
			1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,
			1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,
			1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,
			1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
			1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,
			1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063/*,
			2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,
			2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,
			2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,
			2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,
			2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,
			2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,
			2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,
			2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,
			3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,
			3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,
			3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,
			3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,
			3527,3529,3533,3539,3541,3547,3557,3559,3571,3581,3583,3593,3607,3613,3617,
			3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733,
			3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,
			3877,3881,3889,3907,3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,
			4003,4007,4013,4019,4021,4027,4049,4051,4057,4073,4079,4091,4093,4099,4111,
			4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,4241,
			4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,
			4373,4391,4397,4409,4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,4507,
			4513,4517,4519,4523,4547,4549,4561,4567,4583,4591,4597,4603,4621,4637,4639,
			4643,4649,4651,4657,4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,4759,
			4783,4787,4789,4793,4799,4801,4813,4817,4831,4861,4871,4877,4889,4903,4909,
			4919,4931,4933,4937,4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,5009,
			5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101,5107,5113,5119,5147,
			5153,5167,5171,5179,5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,5281,
			5297,5303,5309,5323,5333,5347,5351,5381,5387,5393,5399,5407,5413,5417,5419,
			5431,5437,5441,5443,5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,5527,
			5531,5557,5563,5569,5573,5581,5591,5623,5639,5641,5647,5651,5653,5657,5659,
			5669,5683,5689,5693,5701,5711,5717,5737,5741,5743,5749,5779,5783,5791,5801,
			5807,5813,5821,5827,5839,5843,5849,5851,5857,5861,5867,5869,5879,5881,5897,
			5903,5923,5927,5939,5953,5981,5987,6007,6011,6029,6037,6043,6047,6053,6067,
			6073,6079,6089,6091,6101,6113,6121,6131,6133,6143,6151,6163,6173,6197,6199,
			6203,6211,6217,6221,6229,6247,6257,6263,6269,6271,6277,6287,6299,6301,6311,
			6317,6323,6329,6337,6343,6353,6359,6361,6367,6373,6379,6389,6397,6421,6427,
			6449,6451,6469,6473,6481,6491,6521,6529,6547,6551,6553,6563,6569,6571,6577,
			6581,6599,6607,6619,6637,6653,6659,6661,6673,6679,6689,6691,6701,6703,6709,
			6719,6733,6737,6761,6763,6779,6781,6791,6793,6803,6823,6827,6829,6833,6841,
			6857,6863,6869,6871,6883,6899,6907,6911,6917,6947,6949,6959,6961,6967,6971,
			6977,6983,6991,6997,7001,7013,7019,7027,7039,7043,7057,7069,7079,7103,7109,
			7121,7127,7129,7151,7159,7177,7187,7193,7207,7211,7213,7219,7229,7237,7243,
			7247,7253,7283,7297,7307,7309,7321,7331,7333,7349,7351,7369,7393,7411,7417,
			7433,7451,7457,7459,7477,7481,7487,7489,7499,7507,7517,7523,7529,7537,7541,
			7547,7549,7559,7561,7573,7577,7583,7589,7591,7603,7607,7621,7639,7643,7649,
			7669,7673,7681,7687,7691,7699,7703,7717,7723,7727,7741,7753,7757,7759,7789,
			7793,7817,7823,7829,7841,7853,7867,7873,7877,7879,7883,7901,7907,7919,7927,
			7933,7937,7949,7951,7963,7993,8009,8011,8017,8039,8053,8059,8069,8081,8087,
			8089,8093,8101,8111,8117,8123,8147,8161,8167,8171,8179,8191,8209,8219,8221,
			8231,8233,8237,8243,8263,8269,8273,8287,8291,8293,8297,8311,8317,8329,8353,
			8363,8369,8377,8387,8389,8419,8423,8429,8431,8443,8447,8461,8467,8501,8513,
			8521,8527,8537,8539,8543,8563,8573,8581,8597,8599,8609,8623,8627,8629,8641,
			8647,8663,8669,8677,8681,8689,8693,8699,8707,8713,8719,8731,8737,8741,8747,
			8753,8761,8779,8783,8803,8807,8819,8821,8831,8837,8839,8849,8861,8863,8867,
			8887,8893,8923,8929,8933,8941,8951,8963,8969,8971,8999,9001,9007,9011,9013,
			9029,9041,9043,9049,9059,9067,9091,9103,9109,9127,9133,9137,9151,9157,9161,
			9173,9181,9187,9199,9203,9209,9221,9227,9239,9241,9257,9277,9281,9283,9293,
			9311,9319,9323,9337,9341,9343,9349,9371,9377,9391,9397,9403,9413,9419,9421,
			9431,9433,9437,9439,9461,9463,9467,9473,9479,9491,9497,9511,9521,9533,9539,
			9547,9551,9587,9601,9613,9619,9623,9629,9631,9643,9649,9661,9677,9679,9689,
			9697,9719,9721,9733,9739,9743,9749,9767,9769,9781,9787,9791,9803,9811,9817,
			9829,9833,9839,9851,9857,9859,9871,9883,9887,9901,9907,9923,9929,9931,9941,
			9949,9967,9973,10007,10009,10037,10039,10061,10067,10069,10079,10091,10093*/};
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
	private final double[] sqrt;
	private final TDiv63Inverse tdiv = new TDiv63Inverse(I_MAX);




	/**
	 * Full constructor.
	 * @param doTDivFirst If true then trial division is done before the Lehman loop.
	 * This is recommended if arguments N are known to have factors < cbrt(N) frequently.
	 * With doTDivFirst=false, this implementation is pretty fast for hard semiprimes.
	 * But the smaller possible factors get, it will become slower and slower.
	 */
	public HartSieveDynamic2(boolean doTDivFirst) {
		this.doTDivFirst = doTDivFirst;
		// Precompute sqrts for all k < I_MAX
		sqrt = new double[I_MAX];
//		sqrt2 = new double[I_MAX];
		for (int i=1; i<I_MAX; i++) {
			sqrt[i] = Math.sqrt(i);
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
		double sieveMult = 1;
//		int primeBound = (int) (.17 * (primes.length) -1);
		int primeBound = 30;
		int sieveBound = (int) (sieveMult * primes[primeBound]);
//		int sieveBound = 36 * 271;
		System.out.println("primes : " + primeBound + " sieve bound : " + sieveBound);
		byte [] factorLength = new byte [2*sieveBound];
		List<BitSet> smoothMatrix = new ArrayList<>();
		Set<Long> smoothNumbers = new HashSet<>();
		Set<Integer> primesUsed = new HashSet<>();
		TreeSet<Level2> bestLevel = new TreeSet<>();
		int maxLevel = 1;
		int sieveCount = 0;
		Level2 empyLevel = new Level2(1);
		bestLevel.add(empyLevel);
		while (smoothCount <= primesUsed.size()) {
			final Level2 level = bestLevel.pollFirst();
			int k = level.level;
			int newColumns = 0;
			int newRows = 0;
			int smooth4Level = 0;
			int step = (k&1) == 1 ? 2 : 1;
			long beginXIndex = (long) (sqrt4N * sqrt[k] + ROUND_UP_DOUBLE);
			beginXIndex += sieveBound * step * level.sieveBeginIndex;
			// bring beginXIndex to some k to value b = kn+1 mod 2^l
			if (step == 2){
				beginXIndex |= 1;
			}
			for (int l = 0; l < sieveBound * step; l++) {
				a = (long) (beginXIndex + l );
				long test = a*a - k * N;
				byte lenght = (byte) (63 - Long.numberOfLeadingZeros(test));
				factorLength[l] = lenght;
			}

			for (int j=primeBound; j > 5 ; j--) 
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
					a = beginXIndex + l;
					long test = a*a - k * N;
					tdiv.setTestLimit(primes[primeBound]);
					int[] factors = smooth.factorizeIfSmooth1(test, primeBound);
					if(factors[primeBound+1] == 1) {
						if ((test / k) * k == test && smoothNumbers.contains(test/k)) {
							System.out.println("Duplicate : " + test/ k);
						}
						else {
							final SortedMultiset<BigInteger> factor = tdiv.factor(BigInteger.valueOf(test));
							System.out.print("a " + a+" kn + 1 " + (k* N) % 8);
							System.out.println(factor);
							smoothNumbers.add(test);
							BitSet factorsMod2 = new BitSet();
							for (int j = 0; j < factors.length; j++) {
								if (factors[j] %2 == 1) {
									factorsMod2.set(j);
									if (!level.primesUsed.contains(j)) {
										level.primesUsed.add(j);
										newColumns++;
									}
									if (!primesUsed.contains(j)) {
										primesUsed.add(j);
									}
								}
							}
							smoothMatrix.add(factorsMod2);
							smoothCount++;	
							smooth4Level++;
							newRows++;
						}
					}
				}
			}
			sieveCount++;
			System.out.println("k : " + k + "  new rows : " + newRows + " new colums " + newColumns);
			level.smoothCountEst = (((double)level.smoothCountInLevel)/(level.sieveBeginIndex+1) + smooth4Level * 2) / 3;
			level.smoothCountInLevel += smooth4Level;
//			level.count = smooth4Level;
			if ((level.level == maxLevel) && (level.sieveBeginIndex == 0)) {
				Level2.smoothCountTrunk += smooth4Level;
				// average smooth count with sieveBeginIndex = 0 to give it a good chance to be selected
				double est = ((double)Level2.smoothCountTrunk) / maxLevel;
				Level2 newLevel = new Level2(++maxLevel, 0, 0, est);
				bestLevel.add(newLevel);
			}
			level.sieveBeginIndex++;
			bestLevel.add(level);
		}
		System.out.println("#sieves : " + sieveCount);
		System.out.println("Work sieveBound * number sieves : " + sieveBound * Math.log(Math.log(sieveBound)) * sieveCount);
		List<BitSet> reducedMatrix = smoothMatrix;
		do {
			smoothMatrix = reducedMatrix;
			reducedMatrix = reduceMatrix(smoothMatrix, primeBound);
		}
		while (reducedMatrix.size() < smoothMatrix.size());
		
		return 1;
	}

	protected List<BitSet> reduceMatrix(List<BitSet> smoothMatrix, int primeBound) {
		int [] colCounts = new int [primeBound];
		List<List<Integer>> rows4Columns = new ArrayList<>();
		for (int l = 0; l < primeBound; l++) {
			rows4Columns.add(new ArrayList<>());
		}		
		for (int j = 0; j < smoothMatrix.size(); j++) {
			BitSet row = smoothMatrix.get(j);
			String line = String.format("%02d", j);
			for (int l = 0; l < primeBound; l++) {
				final boolean rowContainsPrime = row.get(l);
				line += rowContainsPrime ? "o|" : " |";
				if (rowContainsPrime) {
					rows4Columns.get(l).add(j);
					colCounts[l]++;
				}
			}
			System.out.println(line);
		}
		System.out.print(" ");
		for (int i = 0; i < colCounts.length; i++) {
			System.out.print(String.format("%02d", i));			
		}
		System.out.println();
		System.out.print(" ");

		int colEmpty = 0;
		for (int i = 0; i < colCounts.length; i++) {
			if (colCounts[i] <= 2) {
				System.out.print("_" +String.format("%01d", colCounts[i]));	
				colEmpty += colCounts[i] == 0 ? 1 : 0;
			}
			else
				System.out.print(String.format("%02d", colCounts[i]));		
		}
		System.out.println();
		System.out.println("colums non empty : " + (colCounts.length - colEmpty));
		HashSet<Integer> deleteRowsSet = new HashSet<>();
		for (List<Integer> rows4Column : rows4Columns) {
			if (rows4Column.size() <= 2) {
				if (rows4Column.size() == 2) {
					smoothMatrix.get(rows4Column.get(1)).xor(smoothMatrix.get(rows4Column.get(0)));
				}
				if (rows4Column.size() > 0)
					deleteRowsSet.add(rows4Column.get(0));
			}
		}
		System.out.print("rows to delete : ");
		deleteRowsSet.stream().map(i -> i + ",").forEach(System.out::print);
		List<BitSet> reducedMatrix = new ArrayList<>();
		for (int row=0; row < smoothMatrix.size(); row++) {
			if (!deleteRowsSet.contains(row))
				reducedMatrix.add(smoothMatrix.get(row));
		}
		return reducedMatrix;
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
	


}
