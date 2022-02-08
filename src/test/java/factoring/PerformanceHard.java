package factoring;

import java.math.BigInteger;
import java.util.*;

import de.tilman_neumann.jml.factor.FactorAlgorithm;

import de.tilman_neumann.jml.factor.tdiv.TDiv31Barrett;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;

import factoring.hart.*;
import factoring.math.PrimeMath;
import factoring.primes.Primes;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard {


//	2^(36/3) = 2^12 ~ 10^4
	final static int bits = 42	;
	final static int numPrimes = 20000;
//	final static int loop = 9000;
		final static int loop = 20;
	static long[] semiprimes;
	private static boolean readFromFile = false;

	public static void main(String[] args) {
		final PerformanceHard test = new PerformanceHard();
		test.singleFactor();
		//		factorize();
		//		singleFactor();
	}

	private void analyzeMultipliers(){


	}

	private void singleFactor() {
		//		final FactorAlgorithmBase factorizer2 = new factoring.fermat.lehman.Lehman_Fast(false);
		//		final FactorAlgorithmBase factorizer2 = new SquFoF31();
		//		final FactorizationOfLongs factorizer2 = new TrialInvFact2(1 << (bits/2));
		//		final FactorAlgorithmBase factorizer1 = new LehmanSimple();
//		final FactorAlgorithm factorizer1 = new TDiv52InverseFMA();
//		final FactorAlgorithm factorizer2 = new TDiv31Barrett();
//				final FactorAlgorithm factorizer2 = new TDiv23InverseFMA();
		//		final FactorAlgorithmBase factorizer2 = new de.tilman_neumann.jml.factor.lehman.Lehman_Fast(false);
		//		final FactorizationOfLongs factorizer1 = new PollardRhoBrentDouble53();
		//		final FactorAlgorithmBase factorizer1 = new Lehman_FastJones(true);
		//				final FactorAlgorithmBase factorizer1 = new Lehman_Fast6(true);
		//		final FactorizationOfLongs factorizer1 = new TrialMultiplyCorrection(1 << (bits/2));
		//		final FactorAlgorithm factorizer2 = new TDiv63Inverse(1 << (bits/2));
		//		final FactorAlgorithm  factorizer1 = new TDiv63Inverse_NoDoubleCheck_Unroll(1 << (bits/2));
		//		final FactorAlgorithm  factorizer2 = new TDiv63Inverse_NoDoubleCheck_Unroll(1 << (bits/2));
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorFinder(bits, 2.f, false);
		//		final FactorizationOfLongs factorizer2 = new PollardRhoBrentDouble52();
		//		final FactorAlgorithmBase factorizer1 = new LehmanMultiplier6_5_7_11(true);
		//		final FactorAlgorithm factorizer2 = new Lehman_Fast(false);
		//		final FactorAlgorithm factorizer2 = new Lehman_FastOrig(false);
		//		final FactorAlgorithm factorizer2 = new Lehman315(false);
		//		final FactorAlgorithm factorizer2 = new Lehman_Fast(false);
		//		final FactorizationOfLongs factorizer1 = new LehmanMod30(false);
		//		final FactorAlgorithm factorizer1 = new LehmanMultiplier(false);
		//				final FactorAlgorithmBase factorizer2 = new LehmanMidRange(false, 1.);
		//		final FactorAlgorithmBase factorizer2= new LehmanMidRange5(1);
		//		final FactorAlgorithm factorizer2 = new LehmanMidRange7(2, 3);
		//				final FactorAlgorithm factorizer1 = new Hart_TDiv_Race();
		//		final FactorAlgorithm factorizer2 = new Hart_TDiv_Race();
		//		final FactorAlgorithm factorizer1 = new Hart_TDiv_Race2();
		//		final FactorAlgorithm factorizer2 = new de.tilman_neumann.jml.factor.hart.Hart_TDiv_Race();
		//		final FactorAlgorithm factorizer2 = new Lehman_CustomKOrder(true);
		//		final FactorAlgorithm factorizer2 = new Lehman_CustomKOrder(false);
		//				final FactorAlgorithm factorizer1 = new Lehman_CustomKOrderTh(true);
		//		final FactorAlgorithm factorizer2 = new Lehman_CustomKOrderTh3(false);
		//		final FactorAlgorithm factorizer2 = new Lehman_CustomKOrderTh(false);
		//		final FactorAlgorithm factorizer1 = new HartTraining(false, false);
//		final FactorAlgorithm factorizer2 = new Hart_Fast2MultOrig(false);
//		final FactorAlgorithm factorizer2 = new Hart_FastOddFirst(false);
//		final FactorAlgorithm factorizer2 = new Hart_FastNo9(false);
//		final FactorAlgorithm factorizer1 = new Hart_FastNarrow(false);
//		final FactorAlgorithm factorizer1 = new Hart_FastAdjustMap(false);
//		final FactorAlgorithm factorizer1 = new Hart_Fast2Mult4(false);
//		final FactorAlgorithm factorizer2 = new Hart_Fast2MultAdjustALookup(false);
		final FactorAlgorithm factorizer2 = new Hart_Fast2Mult(false, bits);
		final FactorAlgorithm factorizer1 = new HartSmoothMultiplier(false, bits);
//		final FactorAlgorithm factorizer1 = new Hart_FastMult2Analyze(false, bits);
//				final FactorAlgorithm factorizer2 = new SmoothNumbersSieve();
		//		final FactorAlgorithmBase factorizer1 = new LehmanHart(0);
		//		final FactorAlgorithmBase factorizer1 = new LehmanHart2();
		//		final FactorAlgorithm factorizer1 = new Hart_FastNo9(false);
		//		final FactorAlgorithm factorizer1 = new Hart_FastAdjustMap(false);
		//		final FactorAlgorithm factorizer1 = new Hart_FastOdd(false);


		//		final FactorAlgorithm factorizer1 = new Hart_FastTRound(false);
		//		final FactorAlgorithm factorizer1 = new HartSimple4();
		//		final FactorAlgorithm factorizer1 = new HartSimpleMin();
		//		final FactorAlgorithm factorizer2 = new HartSimpleMin();
		//		final FactorAlgorithm factorizer1 = new HartMod8(true);
		//		final FactorAlgorithm factorizer2 = new HartMod8(true);
		//		final FactorAlgorithm factorizer1 = new LehmanMidRange5(1);
		//		final FactorAlgorithmBase factorizer1 = new LehmanMultiplier6_5_7(true);
		semiprimes = Primes.makeSemiPrimesList(bits, numPrimes, readFromFile);
		//		semiprimes = makeSemiPrimesListReal(bits, numPrimes);
//		semiprimes = new long[]{574000473499l};
		final long min1 = test2(factorizer1);
//		final long min2 = test2(factorizer1);

		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, loop);

		final long min2 = test2(factorizer2);

		System.out.println("Speedup : " + min1 / (double) min2);
		System.out.println("Speedup : " + min2 / (double) min1);
		//		final double completeWork = Arrays.asList(HartMod8.foundInStep).stream().flatMapToInt(IntStream::of).sum();
		int completeSolutions = 0;
		int solutionsMod8 = 0;
		final int mod = 4;
		int ignoreCount = 0;
		int hartRange = (1 << (bits/3)) * 3;
		double entryLevel = .1;
		List<HartSmoothMultiplier.Score> scores = new ArrayList<>();
		int goodMultipliers = 0;
		for (int i = 0; i < Hart_FastMult2Analyze.I_MAX; i++) {
			if (Hart_FastMult2Analyze.multiplierScore[i] > entryLevel){
				HartSmoothMultiplier.Score score = new HartSmoothMultiplier.Score(Hart_FastMult2Analyze.multiplierScore[i], i);
				scores.add(score);
				goodMultipliers++;
			}
		}
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

	private String bar(int count, String barChar) {
		double logCount = Math.log(count) * 10;
		String bar = "";
		for (int i=0; i < logCount; i++)
			bar += barChar;
		return bar;
	}
	//	private static void factorize() {
	//		final LehmanSimple factorizer1 = new LehmanSimple();
	//		final FactorAlgorithm factorizer2 = new Lehman_FastOrig(true);
	//		//		final FactorAlgorithmBase factorizer2 = new Lehman_Fast33(true);
	//		semiprimes = makeSemiPrimesList(bits, numPrimes);
	//		factorize(factorizer1);
	//		factorize(factorizer2);
	//	}

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

	public static long test2(FactorAlgorithm factorizer) {
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		findFactors(factorizer, semiprimes, loop);
		//		for (final Entry<BigInteger, Integer> entry : Hart_FastT2.factors.entrySet()) {
//		final Iterator<Entry<BigInteger, Integer>> iterator = Hart_FastT2.factors.entrySet().iterator();
//		final Iterator<Entry<BigInteger, Integer>> iterator2 = Hart_FastT2.factorsExpo.entrySet().iterator();
//		for (int i = 0; i < Hart_FastT2.factors.entrySet().size(); i++) {
//			final Entry<BigInteger, Integer> entry = iterator.next();
//			final Entry<BigInteger, Integer> entry2 = iterator2.next();
//			System.out.println("factor : \t" + entry.getKey() + " : \t" + (entry.getValue() * entry.getKey().longValue()) / (double) (Hart_FastT2.factorNumber)+
//					"\tfactorExpo : \t" + entry2.getKey() + " : \t" + (entry2.getValue() * entry2.getKey().longValue()) / (double) (Hart_FastT2.factorNumber - 1));
//		}
		//		for (int i = 0; i < 1000; i++) {
		//			System.out.println(i   + " \t:"  + Hart_FastT.ks[i] + " \t:"  + factorizer.factor(BigInteger.valueOf(i)));
		//		}

		long minTime = findFactors(factorizer, semiprimes, loop);
		minTime = Math.min(minTime, findFactors(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, findFactors(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, findFactors(factorizer, semiprimes, loop));
		minTime = Math.min(minTime, findFactors(factorizer, semiprimes, loop));
//
		return minTime;
	}


	public static void factorize(FactorizationOfLongs factorizer){
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
	}


	public static void test2(FactorizationOfLongs factorizer){
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
		findFactors(factorizer, semiprimes, loop);
	}

	protected static long factorize(final FactorAlgorithm factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.factor(BigInteger.valueOf(semiprime));
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-45s", factorizer1.getClass().getSimpleName());
		System.out.println(name + " :    \t" +  time);
		return time;
	}

	protected static long findFactors(final FactorAlgorithm factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		int j = 0;
		for (int i = 0; i < loop; i++) {
			for (final Long semiprime : semiprimes) {
				if (semiprime != null)
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

	protected static long factorize(final FactorizationOfLongs factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.factorization(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-50s", factorizer1);
		System.out.println(name + " :    \t" + time);

		return time;
	}

	protected static long findFactors(final FactorizationOfLongs factorizer1, final long[] semiprimes, int loop) {
		final long start = System.nanoTime();
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findFactor(semiprime);
				//				factorizer1.findFactors(semiprime, null);
				//				factorizer1.factorization(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-50s", factorizer1);
		System.out.println(name + " :    \t" + time);

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
