package factoring.smooth;

import org.junit.Test;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import com.google.common.collect.TreeMultiset;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

import factoring.trial.TrialInvFact;
import junit.framework.TestFailure;

public class SmoothNumbersTest {
	
    private static final int _100000 = 10000;

	public static void main(String[] args) {
    	for (int i = 0; i< 10; i++)
    	{
    		long start = System.currentTimeMillis();
    		testSmooth3();
    		long end = System.currentTimeMillis();
    		System.out.println("Time testSmooth3 : " + (end - start));
    		start = System.currentTimeMillis();
    		testSmoothInv();
    		end = System.currentTimeMillis();
    		System.out.println("Time testSmoothI : " + (end - start));
    		start = System.currentTimeMillis();
    		testSmoothByFactor1();
    		end = System.currentTimeMillis();
    		System.out.println("Time testfactor1 : " + (end - start));
//    		start = System.currentTimeMillis();
//    		testSmoothByFactor();
//    		end = System.currentTimeMillis();
//    		System.out.println("Time testSmoothF : " + (end - start));
    	}
    	
    }
	
	@Test
	public void testSmooth(){
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		
		int smoothCount = 0;
		int maybeSmoothCount = 0;
		
		SmoothNumbers smooth = new SmoothNumbers();
		TrialInvFact fact = new TrialInvFact(37);

		for (long i = (1l << 29) + 13; maybeSmoothCount < 1000; i++) {
			TreeMultiset<Long> factors = fact.factorization(i);
			long prod = factors.stream().reduce(1l, (a,b) -> a * b);
			smoothCount += prod == i ? 1 : 0;
			boolean isSmooth = smooth.isSmoothBy3Inv(i, PRIMORIAL);
//			isSmooth = smooth.factorizeIfSmooth(i, 37) > 0;
//			int[] primeFactors = smooth.factorizeIfSmooth1(i, 37);
			int[] primeFactors = smooth.factorizeIfSmooth2(i, 37);
//			boolean isSmooth = primeFactors[primeFactors.length-1] > 0;
			maybeSmoothCount += isSmooth ? 1 : 0;
			
			if (isSmooth)
				assertEquals("i : " + i, prod == i, isSmooth);	
//				final long multiplyPrimeFactors = LongStream.range(0, primeFactors.length-1)
//						.filter(j -> primeFactors[(int)j] > 0).reduce(1l, (a,b) -> a * (long)Math.pow(SmoothNumbers.primes[(int)b], primeFactors[(int)b]));
//				assertEquals("factoring " + i, prod, multiplyPrimeFactors);
		}
		System.out.println("smooth Fraction : " + (0.0 + smoothCount)/ maybeSmoothCount);
	}

	public static int testSmooth3(){
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		SmoothNumbers smooth = new SmoothNumbers();

		int count = 0;
		int smoothCount = 0;
		
		for (long i = 1l << 29; count  < _100000; i++) {
			boolean isSmooth = smooth.isSmoothBy3Modulus(i, PRIMORIAL);
			smoothCount += isSmooth ? 1 : 0;
			count++;
		}
		return smoothCount;
	}
	public static int testSmoothInv(){
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		SmoothNumbers smooth = new SmoothNumbers();

		int count = 0;
		int smoothCount = 0;
		
		for (long i = 1l << 29; count  < _100000; i++) {
			boolean isSmooth = smooth.isSmoothBy3Inv(i, PRIMORIAL);
			smoothCount += isSmooth ? 1 : 0;
			count++;
		}
		return smoothCount;
	}
	public static int testSqrt(){
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		SmoothNumbers smooth = new SmoothNumbers();

		int count = 0;
		int sqrts = 0;
		for (long i = 1l << 29; count  < _100000; i++) {
			long sqrt = (long)Math.sqrt(i);
			long sqrtSqare = sqrt * sqrt;
			if (sqrtSqare == i)
				sqrts++;
//			smoothCount += isSmooth ? 1 : 0;
			count++;
		}
		return sqrts;
	}

	public static void testSmooth5(){
		final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 
		SmoothNumbers smooth = new SmoothNumbers();

		int smoothCount = 0;
		for (long i = 1l << 29; smoothCount  < _100000; i++) {
			boolean isSmooth = smooth.isSmoothBy5Modulus(i, PRIMORIAL);
			smoothCount += isSmooth ? 1 : 0;
		}
	}

	public static void testSmoothByFactor(){
		TrialInvFact fact = new TrialInvFact(37);

		int smoothCount = 0;
		for (long i = 1l << 29; smoothCount < _100000; i++) {
			TreeMultiset<Long> factors = fact.factorization(i);
			Long prod = factors.stream().reduce(1l, (a,b) -> a * b);
			smoothCount += prod == i ? 1 : 0;
		}
	}
	public static void testSmoothByFactor1(){
		SmoothNumbers smooth = new SmoothNumbers();

		int smoothCount = 0;
		for (long i = 1l << 29; smoothCount < _100000; i++) {
			int[] factors = smooth.factorizeIfSmooth1(i, 37);
			smoothCount += factors[factors.length-1] > 0 ? 1 : 0;
		}
	}

}
