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
package factoring;

import static de.tilman_neumann.jml.base.BigIntConstants.I_0;
import static de.tilman_neumann.jml.base.BigIntConstants.I_1;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.TestMode;
import de.tilman_neumann.jml.factor.TestNumberNature;
import de.tilman_neumann.jml.factor.TestsetGenerator;
import de.tilman_neumann.jml.factor.hart.Hart_Fast;
import de.tilman_neumann.jml.factor.hart.Hart_TDiv_Race;
import de.tilman_neumann.jml.factor.hart.Hart_TDiv_Race_Unsafe;
import de.tilman_neumann.jml.factor.lehman.Lehman_Fast;
import de.tilman_neumann.jml.factor.lehman.Lehman_Fast2;
import de.tilman_neumann.jml.factor.lehman.Lehman_Fast3;
import de.tilman_neumann.jml.factor.pollardRho.PollardRhoBrentMontgomeryR64Mul63;
import de.tilman_neumann.util.ConfigUtil;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.SortedMultiset_BottomUp;
import de.tilman_neumann.util.TimeUtil;

/**
 * Main class to compare the performance of factor algorithms.
 * @author Tilman Neumann
 */
@SuppressWarnings("unused") // suppress warnings on unused imports
public class FactorizerTest {
	private static final Logger LOG = Logger.getLogger(FactorizerTest.class);

	// algorithm options
	/** number of test numbers */
	private static final int N_COUNT = 100000;
	/** the bit size of N to start with */
	private static final int START_BITS = 30;
	/** the increment in bit size from test set to test set */
	private static final int INCR_BITS = 1;
	/** maximum number of bits to test (no maximum if null) */
	private static final Integer MAX_BITS = null;
	/** each algorithm is run REPEATS times for each input in order to reduce GC influence on timings */
	private static final int REPEATS = 1;
	/** Nature of test numbers */
	private static final TestNumberNature TEST_NUMBER_NATURE = TestNumberNature.QUITE_HARD_SEMIPRIMES;
	/** Test mode */
	private static final TestMode TEST_MODE = TestMode.FIRST_FACTOR;

	/**
	 * Algorithms to compare. Non-static to permit to use Loggers in the algorithm constructors.
	 */
	private final FactorAlgorithm[] algorithms;

	public FactorizerTest() {
		algorithms = new FactorAlgorithm[] {

				// Trial division
				//new TDiv31(),
				//new TDiv31Preload(),
				//			new TDiv31Inverse(), // Fastest algorithm for N <= 24 bit
				//			new TDiv31Inverse_NoDoubleCheck(),
				//			new TDiv31Inverse_NoDoubleCheck_Unroll(),
				//			new TDiv63Inverse(1<<21),
				//			new TDiv63Inverse_NoDoubleCheck(1<<21),
				//			new TDiv63Inverse_NoDoubleCheck_Unroll(1<<21), // very good trial division algorithm

				// Hart's one line factorizer
				//new Hart_Simple(),
				new Hart_Fast(false), // best algorithm for hard semiprimes
				//			new Hart_Fast(true),
				new Hart_TDiv_Race(), // best safe algorithm for any N with 25 to 49 bits, best for moderate semiprimes <= 44 bit
				new Hart_TDiv_Race_Unsafe(), // best algorithm for moderate semiprimes >= 45 bit (but fails for some N having small factors)

				new factoring.hart.Hart_TDiv_Race2(),
				new factoring.hart.Hart_FastT(false),

				// Lehman
				//new Lehman_Simple(false),
				//new Lehman_Smith(false),
				new Lehman_Fast(false), // the variant implemented by bsquared
				//			new Lehman_Fast(true),
				new Lehman_Fast2(false), // best Lehman for moderate semiprimes
				new Lehman_Fast3(false), // best Lehman for hard semiprimes

				// PollardRho
				//new PollardRho(),
				//new PollardRho_ProductGcd(),
				//new PollardRhoBrent(),
				//new PollardRho31(),
				//new PollardRhoBrent31(),
				//new PollardRhoBrentMontgomery63(), // first long version, not optimized any further
				new PollardRhoBrentMontgomeryR64Mul63(), // best algorithm for N from 50 to 56 bit
				//			new PollardRhoBrentMontgomery64(), // best algorithm for N from 57 to 62 bit

				// SquFoF variants
				// * pretty good, but never the best algorithm
				// * SquFoF31 works until 52 bit and is faster there than SquFoF63
				// * best multiplier sequence = 1680 * {squarefree sequence}
				// * best stopping criterion = O(5.th root(N))
				//			new SquFoF63(),
				//new SquFoF31(),
				//			new SquFoF31Preload(),

				// CFrac
				// * never the best algorithm: SquFoF63 is better for N <= 65 bit, SIQS is better for N >= 55 bits
				// * stopRoot, stopMult: if big enough, then a second k is rarely needed; (5, 1.5) is good
				// * TDiv_CF01 is good for N < 80 bits; for N > 90 bit we need TDiv_CF02
				// * ksAdjust: Must be <=3 for N=20bit, <=6 for N=30 bit etc. // TODO this implies some optimization potential
				//new CFrac(true, 5, 1.5F, 0.152F, 0.253F, new TDiv_CF01(), 10, new MatrixSolver01_Gauss(), 5, false),
				//new CFrac(true, 5, 1.5F, 0.152F, 0.253F, new TDiv_CF02(), 10, new MatrixSolver01_Gauss(), 5, false),
				//			new CFrac63(true, 5, 1.5F, 0.152F, 0.25F, new TDiv_CF63_01(), 10, new MatrixSolver01_Gauss(), 3),
				//new CFrac63(true, 5, 1.5F, 0.152F, 0.25F, new TDiv_CF63_02(), 10, new MatrixSolver01_Gauss(), 12),

				// ECM
				//			new EllipticCurveMethod(),

				// SIQS:
				// * best until 220 bit: Sieve03gU + smallPowers + TDiv1L + Gauss
				// * best for 230, 240 bit: Sieve03gU + smallPowers + TDivnL + BL
				// * best for >= 250 bit: (Sieve03gU or SingleBlockHybridSieve) + (noPowers or smallPowers) + (TDiv2L or TDivnL) + BL
				//new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SimpleSieve(), new TDiv_QS_1Large(), 10, new MatrixSolver01_Gauss(), false),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),

				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_2Large_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_nLarge(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),

				// sieving with prime powers: best sieve for small N!
				//			new SIQS(0.32F, 0.37F, null, null, new PowerOfSmallPrimesFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),
				//			new SIQS(0.32F, 0.37F, null, null, new PowerOfSmallPrimesFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new AllPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),

				// segmented sieve: slower; best block size is about 32k
				//			new SIQS(0.32F, 0.385F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockSieve(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
				//			new SIQS(0.32F, 0.385F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockSieveU(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.41F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockSieve(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),

				// hybrid sieves:
				// * single block hybrid is level with Sieve03g
				// * best block size is about 32k (my processor has 16kB L1 cache, 128kB L2-cache per thread)
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockHybridSieve01(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockHybridSieveU(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockHybridSieve01(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
				//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockHybridSieveU(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),

				// Multi-threaded SIQS:
				// * 4/6 threads takes over at N around 100 bit (more exact estimates: 4 threads best for N>=88 bits, 6 threads for N>=112 bits)
				// * we need 0.14 < maxQRestExponent < 0.2; everything else is prohibitive; use null for dynamic determination
				// * BlockLanczos is better than Gauss solver for N > 200 bit
				//			new PSIQS(0.32F, 0.37F, null, null, 6, new NoPowerFinder(), new MatrixSolver02_BlockLanczos(), false),
				//			new PSIQS_U(0.32F, 0.37F, null, null, 6, new NoPowerFinder(), new MatrixSolver02_BlockLanczos(), true),
				//			new PSIQS_U(0.32F, 0.37F, null, null, 6, new PowerOfSmallPrimesFinder(), new MatrixSolver02_BlockLanczos(), true),
				//			new PSIQS_U(0.32F, 0.37F, null, null, 6, new AllPowerFinder(), new MatrixSolver02_BlockLanczos(), true),
				//			new PSIQS_SBH_U(0.32F, 0.37F, null, null, 32768, 6, new PowerOfSmallPrimesFinder(), new MatrixSolver02_BlockLanczos(), true), // best for large N

				// Combination of best algorithms for all factor argument sizes
				//			new CombinedFactorAlgorithm(6, 1<<16, true, false),
		};
	}

	private void testRange(int bits) {
		final BigInteger N_min = I_1.shiftLeft(bits-1);
		// Compute test set
		final BigInteger[] testNumbers = TestsetGenerator.generate(N_COUNT, bits, TEST_NUMBER_NATURE);
		final BigInteger[] factors = new BigInteger[N_COUNT];
		@SuppressWarnings("unchecked")
		final
		SortedMultiset<BigInteger>[] factorSetArray = new SortedMultiset_BottomUp[N_COUNT];

		LOG.info("Test N with " + bits + " bits, i.e. N >= " + N_min);

		// take REPEATS timings for each algorithm to be quite sure that one timing is not falsified by garbage collection
		final TreeMap<Long, List<FactorAlgorithm>> ms_2_algorithms = new TreeMap<>();
		for (int i=0; i<REPEATS; i++) {
			for (final FactorAlgorithm algorithm : algorithms) {
				// exclude special size implementations
				final String algName = algorithm.getName();
				if (bits<53 && algName.startsWith("SIQS")) continue; // unstable for smaller N // TODO the bound has been much smaller some time ago
				if (bits<57 && algName.startsWith("PSIQS")) continue; // unstable for smaller N
				if (bits>98 && algName.startsWith("CFrac63")) continue; // unstable for N>98 bits
				if (bits>52 && algName.startsWith("SquFoF31")) continue; // int implementation
				if (bits>60 && algName.startsWith("Lehman")) continue;
				if (bits>31 && algName.startsWith("TDiv31")) continue; // int implementation
				if (bits>31 && algName.startsWith("PollardRho31")) continue; // long implementation
				if (bits>42 && algName.startsWith("TDiv63Inverse")) continue; // not enough primes stored

				System.gc(); // create equal conditions for all algorithms

				int failCount = 0;
				BigInteger failExample = null;
				long duration;
				switch (TEST_MODE) {
				case FIRST_FACTOR: {
					// test performance
					final long startTimeMillis = System.currentTimeMillis();
					for (int j=0; j<N_COUNT; j++) {
						try {
							factors[j] = algorithm.findSingleFactor(testNumbers[j]);
						} catch (final ArithmeticException e) {
							LOG.error("Algorithm " + algorithm + " threw Exception while searching for a factor of N=" + testNumbers[j] + ": " + e);
						}
					}
					final long endTimeMillis = System.currentTimeMillis();
					duration = endTimeMillis - startTimeMillis; // duration in ms
					//LOG.debug("algorithm " + algName + " finished test set with " + bits + " bits");

					// verify
					for (int j=0; j<N_COUNT; j++) {
						final BigInteger N = testNumbers[j];
						final BigInteger factor = factors[j];
						// test correctness
						if (factor==null || factor.equals(I_0) || factor.equals(I_1) || factor.mod(N).equals(I_0)) {
							//LOG.error("FactorAlgorithm " + algorithm.getName() + " did not find a factor of N=" + N + ", it returned " + factor);
							failExample = N;
							failCount++;
						} else {
							// not null, not trivial -> test division
							final BigInteger[] test = N.divideAndRemainder(factor);
							if (!test[1].equals(I_0)) {
								//LOG.error("FactorAlgorithm " + algorithm.getName() + " returned " + factor + ", but this is not a factor of N=" + N);
								failExample = N;
								failCount++;
							}
						}
					}
					break;
				}
				case PRIME_FACTORIZATION: {
					// test performance
					final long startTimeMillis = System.currentTimeMillis();
					for (int j=0; j<N_COUNT; j++) {
						try {
							factorSetArray[j] = algorithm.factor(testNumbers[j]);
						} catch (final ArithmeticException e) {
							LOG.error("Algorithm " + algorithm.getName() + " threw Exception while factoring N=" + testNumbers[j] + ": " + e);
							factorSetArray[j] = null; // to have correct fail count
						}
					}
					final long endTimeMillis = System.currentTimeMillis();
					duration = endTimeMillis - startTimeMillis; // duration in ms
					//LOG.debug("algorithm " + algName + " finished test set with " + bits + " bits");

					// verify
					for (int j=0; j<N_COUNT; j++) {
						final BigInteger N = testNumbers[j];
						final SortedMap<BigInteger, Integer> factorSet = factorSetArray[j];
						// test correctness
						if (factorSet!=null) {
							for (final BigInteger factor : factorSet.keySet()) {
								if (factor==null || factor.equals(I_0) || factor.equals(I_1) || factor.mod(N).equals(I_0)) {
									//LOG.error("FactorAlgorithm " + algorithm.getName() + " did not find a factor of N=" + N + ", it returned " + factor);
									failExample = N;
									failCount++;
								} else {
									// not null, not trivial -> test division
									final BigInteger[] test = N.divideAndRemainder(factor);
									if (!test[1].equals(I_0)) {
										//LOG.error("FactorAlgorithm " + algorithm.getName() + " returned " + factor + ", but this is not a factor of N=" + N);
										failExample = N;
										failCount++;
									}
								}
							}
						} else {
							failCount++;
						}
					}
					break;
				}
				default: throw new IllegalArgumentException("TestMode = " + TEST_MODE);
				}

				List<FactorAlgorithm> algList = ms_2_algorithms.get(duration);
				if (algList==null) algList = new ArrayList<>();
				algList.add(algorithm);
				ms_2_algorithms.put(duration, algList);
				if (failCount>0) {
					LOG.error("FactorAlgorithm " + algorithm.getName() + " failed at " + failCount + "/" + N_COUNT + " test numbers, e.g. for N = " + failExample);
				}
			}
		}

		// log best algorithms first
		int rank=1;
		for (final long ms : ms_2_algorithms.keySet()) {
			final List<FactorAlgorithm> algList = ms_2_algorithms.get(ms);
			int j=0;
			for (final FactorAlgorithm algorithm : algList) {
				final String durationStr = TimeUtil.timeStr(ms);
				LOG.info("#" + rank + ": Algorithm " + algorithm.getName() + " took " + durationStr);
				j++;
			}
			rank += j;
		}
	}

	/**
	 * Test factor algorithms for sets of factor arguments of growing size and report timings after each set.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		ConfigUtil.initProject();
		final FactorizerTest testEngine = new FactorizerTest();
		int bits = START_BITS;
		while (true) {
			// test N with the given number of bits, i.e. 2^(bits-1) <= N <= (2^bits)-1
			testEngine.testRange(bits);
			bits += INCR_BITS;
			if (MAX_BITS!=null && bits > MAX_BITS) break;
		}
	}
}