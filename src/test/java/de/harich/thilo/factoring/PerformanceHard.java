package de.harich.thilo.factoring;

import de.harich.thilo.factoring.hart.HartFactorizationAlgorithm;
import de.harich.thilo.factoring.hart.calculator.*;
import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.tdiv.TDiv31Barrett;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import factoring.hart.HartSmoothMultiplier;
import factoring.hart.Hart_Fast2Mult2;
import factoring.math.PrimeMath;
import factoring.primes.Primes;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard {


//	2^(36/3) = 2^12 ~ 10^4
	final static int bits = 41	;
	final static int numPrimes = 20000;
//	final static int loop = 9000;
		final static int loop = 1;
	static long[] semiprimes;

    public static void main(String[] args) {
		final PerformanceHard test = new PerformanceHard();
		test.singleFactor();
		//		factorize();
		//		singleFactor();
	}

	private void analyzeMultipliers(){


	}

	private void singleFactor() {

//        final HartFactorizer factorizer1 = new HartFactorizer(new NativeDatatypeCalculator(1));
//        final HartFactorizer factorizer1 = new HartFactorizer(new NativeDatatypeCalculator(3*3*5*7));
//        final HartFactorizationAlgorithm factorizer1 = new HartFactorizationAlgorithm(new Vector3MultiplierCalculator());
//        final HartFactorizer factorizer2 = new HartFactorizer(new Simd2ArraysCalculator());
//        final HartFactorizationAlgorithm factorizer1 = new HartFactorizationAlgorithm(new Mod8GroupsMultiplierCalculator(false, 4096));
        final HartFactorizationAlgorithm factorizer1 = new HartFactorizationAlgorithm(new AdjustXModPow2Calculator());
//        final HartFactorizationAlgorithm factorizer2 = new HartFactorizationAlgorithm(new Vector3MultiplierCalculator());
//        final HartFactorizationAlgorithm factorizer1 = new HartFactorizationAlgorithm(new VectorNo9Calculator());
//        final HartFactorizationAlgorithm factorizer1 = new HartFactorizationAlgorithm(new RemainderGroupArraysCalculator());
//        final HartFactorizer factorizer2 = new HartFactorizer(new Simd2ArraysCalculator());
//        final FactorAlgorithm factorizer2 = new Hart_FastFmaSimple(false);
//        final HartFactorizer factorizer2 = new HartFactorizer(new FMACalculator());
//		final FactorAlgorithm factorizer2 = new Hart_Fast2Mult(false, bits);
		final FactorAlgorithm factorizer2 = new Hart_Fast2Mult2(false);

        boolean readFromFile = true;
        semiprimes = Primes.makeSemiPrimesList(bits, numPrimes, readFromFile);
//        semiprimes = Primes.makeNumberList(bits, 20000, numPrimes, readFromFile);
//				semiprimes = makeSemiPrimesListReal(bits, numPrimes);
//		semiprimes = new long[]{574000473499l};
		final long min1 = factorize(factorizer1);
//		final long min2 = test2(factorizer1);

		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, loop);

		final long min2 = factorize(factorizer2);

		System.out.println("Speedup : " + min1 / (double) min2);
		System.out.println("Speedup : " + min2 / (double) min1);
		//		final double completeWork = Arrays.asList(HartMod8.foundInStep).stream().flatMapToInt(IntStream::of).sum();

//        if (factorizer1.getCalculator() instanceof VectorFillupGroupsCalculator) {
//            VectorFillupGroupsCalculator calculator = ((VectorFillupGroupsCalculator) factorizer1.getCalculator());
//            for (int i = 0; i < 4 ; i++) {
////                for (int remainder = 0; remainder < MODULUS_FOR_ADJUST_X; remainder++) {
////                    int multiplier = MODULUS_FOR_ADJUST_X * i + remainder;
////                    final int count = calculator.multiplierFactorizationCount[i].factorisations.size();
////                    System.out.println(calculator.multiplierFactorizationCount[i].multiplier + " : " + count);
//                final int count = calculator.multiplierFactorizationCount[0][i];
//                System.out.println(i + " : " + count);
////                }
//            }
//        }
//        if (factorizer2.getCalculator() instanceof SimdCalculator) {
//            SimdCalculator calculator = ((SimdCalculator) factorizer2.getCalculator());
//            for (int i = 0; i < calculator.mod; i++) {
//                System.out.println("mod : " + i + " : " + calculator.solutionsMod[i]);
//            }
//        }

		int completeSolutions = 0;
		int solutionsMod8 = 0;
		final int mod = 4;
		int ignoreCount = 0;
		int hartRange = (1 << (bits/3)) * 3;
		double entryLevel = .1;
		List<HartSmoothMultiplier.Score> scores = new ArrayList<>();
		int goodMultipliers = 0;
//		for (int i = 0; i < Hart_FastMult2Analyze.I_MAX; i++) {
//			if (Hart_FastMult2Analyze.multiplierScore[i] > entryLevel){
//				HartSmoothMultiplier.Score score = new HartSmoothMultiplier.Score(Hart_FastMult2Analyze.multiplierScore[i], i);
//				scores.add(score);
//				goodMultipliers++;
//			}
//		}
		System.out.println("Found : " + goodMultipliers + " good multipliers");
		HartSmoothMultiplier.Score[] scoreArr = scores.toArray(new HartSmoothMultiplier.Score[goodMultipliers]);
		Arrays.sort(scoreArr, (a, b) -> (int) Math.signum(b.score - a.score));

		FactorAlgorithm factorizer = new TDiv31Barrett();
		for (int i = 0; i < goodMultipliers; i++) {
			System.out.println(scoreArr[i].index + " : " + scoreArr[i].score + " : " + factorizer.factor(BigInteger.valueOf(scoreArr[i].index)));
		}

		System.out.println("ignore     : \t" + ignoreCount);
		System.out.println("speedup    : \t" + 1000/ (1000.0 - ignoreCount) );

		System.out.println("sol over all  : \t" + completeSolutions);


//		for (int i = 0; i < Hart_Fast2Mult2Analyze.modCase.length; i++) {
//			final double work = (Hart_Fast2Mult2Analyze.modCase[i] +0.0) / completeWork;
//			System.out.println("mod case " + i + " : " + work + " avg k : " + (Hart_Fast2Mult2Analyze.modCaseK[i] / Hart_Fast2Mult2Analyze.modCase[i]));
//		}
		//		System.out.println("fma    : " + factorizer1.fma);
		//		System.out.println("adjust : " + factorizer1.adjust);
		//		System.out.println("test   : " + factorizer1.test);
		//		System.out.println("sqrt   : " + factorizer1.sqrtTime);
	}


	public static long factorize(HartFactorizationAlgorithm factorizer) {
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		long minTime = factorize(factorizer, semiprimes, loop);
//        if (factorizer.getCalculator() instanceof VectorGoodMultiplierCalculator) {
//            VectorGoodMultiplierCalculator mod8GroupMultiplierCalculator = (VectorGoodMultiplierCalculator) factorizer.getCalculator();
//            mod8GroupMultiplierCalculator.calculateMultipliersFromCounts(semiprimes);
//        }


		minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));

		return minTime;
	}

    public static long factorize(FactorAlgorithm factorizer) {
        final long start = System.currentTimeMillis();
        final long end = System.currentTimeMillis();
        System.out.println("time for setup : " + (end - start));

        long minTime = factorize(factorizer, semiprimes, loop);
        minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
        minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
        minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));
        minTime = Math.min(minTime, factorize(factorizer, semiprimes, loop));

        return minTime;
    }

    protected static long factorize(final HartFactorizationAlgorithm factorizer1, final long[] semiprimes, int loop) {
        final long start = System.nanoTime();
        int j = 0;
        int mod = 4;
        for (int i = 0; i < loop; i++) {
            for (final Long semiprime : semiprimes) {
//                if (semiprime % mod == 3)
                    factorizer1.findSingleFactor(semiprime);
//                }
            }
        }
        final long time = System.nanoTime() - start;
        final String name = String.format("%-45s", factorizer1.getClass().getSimpleName() + "." + factorizer1.getCalculator().getClass().getSimpleName());
        System.out.println(name + " :    \t" +  time);
        return time;
    }

    protected static long factorize(final FactorAlgorithm factorizer1, final long[] semiprimes, int loop) {
        final long start = System.nanoTime();
        int j = 0;
        int mod = 4;
        for (int i = 0; i < loop; i++) {
            for (final Long semiprime : semiprimes) {
//                if (semiprime % mod == 1) {
//                if (semiprime != null)
                    factorizer1.findSingleFactor(BigInteger.valueOf(semiprime));
                //                factorizer1.factor(BigInteger.valueOf(semiprime));
//				if (++j % 10000 == 0){
//					System.out.println((100.0 * j / numPrimes) + "% of factorization done. factorized " + j + " semi primes");
//				}
            }
        }
        final long time = System.nanoTime() - start;
        final String name = String.format("%-45s", factorizer1.getClass().getSimpleName());
        System.out.println(name + " :    \t" +  time);
        return time;
    }

	public static long[] makeSemiPrimesListReal(int bits, int numPrimes) {
		final long[] semiPrimes = new long[numPrimes];
		final int limit = (int) Math.pow(2,bits / 3.0);
		final TDiv63Inverse factorizer = new TDiv63Inverse(limit);

		final long offset = 0;
		//		final Random rnd = new Random();
		//		final long offset = Math.abs(rnd.nextInt());
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




}
