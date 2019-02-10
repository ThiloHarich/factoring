package factoring;

import java.math.BigInteger;
import java.util.Random;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.hart.Hart_TDiv_Race;
import factoring.hart.HartMod8;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard {

	final static int bits = 45;
	final static int numPrimes = 5550;
	final static int loop = 18;
	static long[] semiprimes;

	public static void main(String[] args) {
		singleFactor();
		//		factorize();
	}

	private static void singleFactor() {
		//		final FactorAlgorithmBase factorizer2 = new factoring.fermat.lehman.Lehman_Fast(false);
		//		final FactorAlgorithmBase factorizer2 = new SquFoF31();
		//		final FactorizationOfLongs factorizer2 = new TrialInvFact2(1 << (bits/2));
		//		final FactorAlgorithmBase factorizer1 = new LehmanSimple();
		//		final FactorAlgorithmBase factorizer1 = new PollardRhoBrentMontgomery63();
		//		final FactorAlgorithmBase factorizer2 = new de.tilman_neumann.jml.factor.lehman.Lehman_Fast(false);
		//		final FactorizationOfLongs factorizer1 = new PollardRhoBrentDouble53();
		//		final FactorAlgorithmBase factorizer1 = new Lehman_FastJones(true);
		//				final FactorAlgorithmBase factorizer1 = new Lehman_Fast6(true);
		//		final FactorizationOfLongs factorizer1 = new TrialInvFact2(1 << (bits/3));
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorFinder(bits, 2.f, false);
		//		final FactorizationOfLongs factorizer2 = new PollardRhoBrentDouble52();
		//		final FactorAlgorithmBase factorizer1 = new LehmanMultiplier6_5_7_11(true);
		//		final FactorAlgorithmBase factorizer2 = new de.tilman_neumann.jml.factor.lehman.Lehman_Fast(false);
		//				final FactorAlgorithmBase factorizer2 = new LehmanMidRange(false, 1.);
		//		final FactorAlgorithmBase factorizer2= new LehmanMidRange5(1);
		//		final FactorAlgorithm factorizer2 = new LehmanMidRange7(0, 1);
		final FactorAlgorithm factorizer2 = new Hart_TDiv_Race();
		//		final FactorAlgorithm factorizer1 = new factoring.hart.Hart_TDiv_Race();
		//		final FactorAlgorithm factorizer2 = new Hart_Fast(true);
		//		final FactorAlgorithmBase factorizer1 = new LehmanHart(0);
		//		final FactorAlgorithmBase factorizer1 = new LehmanHart2();
		final FactorAlgorithm factorizer1 = new HartMod8(true);
		//		final FactorAlgorithm factorizer1 = new HartMod8(true);
		//		final FactorAlgorithm factorizer2 = new HartMod8(true);
		//		final FactorAlgorithmBase factorizer1 = new LehmanMidRange2(false, 1.4);
		//		final FactorAlgorithmBase factorizer1 = new LehmanMultiplier6_5_7(true);
		semiprimes = makeSemiPrimesList(bits, numPrimes);
		test2(factorizer1);

		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, loop);

		test2(factorizer2);
		//		final double completeWork = Arrays.asList(HartMod8.foundInStep).stream().flatMapToInt(IntStream::of).sum();
		//		for (int i = 0; i < HartMod8.foundInStep.length; i++) {
		//			final double work = (HartMod8.foundInStep[i] +0.0) / completeWork;
		//			System.out.println("step " + i + " : " + work);
		//		}
		//		System.out.println("loop 6k      first : " + factorizer2.loop_6_1);
		//		System.out.println("loop 6      ground : " + factorizer2.loop_ground);
		//		System.out.println("loop 6k     second : " + factorizer2.loop_6_2);
		//		System.out.println("loop 6k + 3 ground : " + factorizer2.loop_3);
	}
	//	private static void factorize() {
	//		final LehmanSimple factorizer1 = new LehmanSimple();
	//		final FactorAlgorithm factorizer2 = new Lehman_FastOrig(true);
	//		//		final FactorAlgorithmBase factorizer2 = new Lehman_Fast33(true);
	//		semiprimes = makeSemiPrimesList(bits, numPrimes);
	//		factorize(factorizer1);
	//		factorize(factorizer2);
	//	}

	public static void factorize(FactorAlgorithm factorizer) {
		final long start = System.currentTimeMillis();
		final long end = System.currentTimeMillis();
		System.out.println("time for setup : " + (end - start));

		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
		factorize(factorizer, semiprimes, loop);
	}

	public static void test2(FactorAlgorithm factorizer) {
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
				factorizer1.findFactors(semiprime, null);
				factorizer1.factorization(semiprime);
			}
		}
		final long time = System.nanoTime() - start;
		final String name = String.format("%-50s", factorizer1);
		System.out.println(name + " :    \t" + time);

		return time;
	}

	public static long[] makeSemiPrimesList(int bits, int numPrimes) {
		final long[] semiPrimes = new long[numPrimes];
		for (int i=0; i< numPrimes; i++)
		{
			final Random rnd = new Random();
			final int smallFactorBits = (bits / 4 )  +  rnd.nextInt(bits / 4 ) ;
			//			final int smallFactorBits = (bits / 2 );
			//			final int smallFactorBits = (bits / 3) - 2;

			final BigInteger fact1 = BigInteger.probablePrime(smallFactorBits, rnd);
			final BigInteger fact2 = BigInteger.probablePrime(bits - smallFactorBits, rnd);
			semiPrimes[i] = fact1.longValue() * fact2.longValue();
		}

		return semiPrimes;
	}




}
