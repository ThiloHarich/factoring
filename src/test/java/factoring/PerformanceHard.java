package factoring;

import java.math.BigInteger;
import java.util.Random;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.hart.Hart_Fast;
import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import factoring.fermat.FermatXSquareMultiplied;
import factoring.fermat.FermatXSquareShift;

import factoring.fermat.FermatXNonSquare;
import factoring.hart.HartTwoSieves;
import factoring.hart.HartTwoSievesSkip;
import factoring.hart.Hart_Fast2MultOrig;
import factoring.hart.Hart_FastFmaSimple;
import factoring.hart.sieve.HartSieve;
import factoring.hart.sieve.HartSieveBitSet;
import factoring.hart.sieve.HartSieveDynamic;
import factoring.hart.sieve.HartSieveDynamic2;
import factoring.hart.sieve.HartSieveFindK;
import factoring.hart.smooth.HartSmooth2;
import factoring.math.PrimeMath;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard {


	final static int bits = 40;
	final static int numPrimes = 4400;
	final static int loop = 90;
	//	final static int loop = 1;
	static long[] semiprimes;

	public static void main(String[] args) {
		final PerformanceHard test = new PerformanceHard();
		test.singleFactor();
		//		factorize();
		//		singleFactor();
	}

	private void singleFactor() {
		//		final FactorAlgorithmBase factorizer2 = new factoring.fermat.lehman.Lehman_Fast(false);
		//		final FactorAlgorithmBase factorizer2 = new SquFoF31();
		//		final FactorizationOfLongs factorizer2 = new TrialInvFact2(1 << (bits/2));
		//		final FactorAlgorithmBase factorizer1 = new LehmanSimple();
		//		final FactorAlgorithm factorizer1 = new TDiv63Inverse(1 << (bits/2));
		//		final FactorAlgorithm factorizer2 = new TDiv63Inverse(1 << (bits/2));
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
		//		final FactorAlgorithm factorizer1 = new HartTraining(false, false);
		final FactorAlgorithm factorizer2 = new Hart_Fast(false);
				final FactorAlgorithm factorizer1 = new Hart_FastFmaSimple(false);
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
		semiprimes = makeSemiPrimesList(bits, numPrimes);
		//		semiprimes = makeSemiPrimesListReal(bits, numPrimes);
		final long min1 = test2(factorizer1);

		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, loop);

		final long min2 = test2(factorizer2);

		System.out.println("Speedup : " + min1 / (double) min2);
		System.out.println("Speedup : " + min2 / (double) min1);
		//		final double completeWork = Arrays.asList(HartMod8.foundInStep).stream().flatMapToInt(IntStream::of).sum();
		//		for (int i = 0; i < HartMod8.foundInStep.length; i++) {
		//			final double work = (HartMod8.foundInStep[i] +0.0) / completeWork;
		//			System.out.println("step " + i + " : " + work);
		//		}
		//		System.out.println("fma    : " + factorizer1.fma);
		//		System.out.println("adjust : " + factorizer1.adjust);
		//		System.out.println("test   : " + factorizer1.test);
		//		System.out.println("sqrt   : " + factorizer1.sqrtTime);
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
		for (int i = 0; i < loop; i++) {
			for (final long semiprime : semiprimes) {
				factorizer1.findSingleFactor(BigInteger.valueOf(semiprime));
				//                factorizer1.factor(BigInteger.valueOf(semiprime));
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

	public static long[] makeSemiPrimesList(int bits, int numPrimes) {
		final long[] semiPrimes = new long[numPrimes];
		Random rnd = new Random(872);
		for (int i=0; i< numPrimes; i++)
		{
			final int smallFactorBits = (bits / 2 );
			//			final int smallFactorBits = (bits / 3) - 2;

			rnd = new Random();
			final BigInteger fact1 = BigInteger.probablePrime(smallFactorBits, rnd);
			final int bigFactorBits = bits - smallFactorBits;
			rnd = new Random();
			final BigInteger fact2 = BigInteger.probablePrime(bigFactorBits, rnd);
			semiPrimes[i] = fact1.longValue() * fact2.longValue();
		}

		return semiPrimes;
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
