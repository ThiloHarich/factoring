package factoring.smooth;

import java.math.BigInteger;
import java.util.*;

import de.tilman_neumann.jml.factor.tdiv.TDiv63Inverse;
import de.tilman_neumann.util.SortedMultiset;
import factoring.math.PrimeMath;

import static java.lang.Math.ceil;
import static java.lang.Math.sqrt;
import static java.lang.Math.log;
import static org.junit.Assert.assertTrue;

/**
 * An algorithm for determining if a number is smooth over some (prime) Factors.
 * It is based on the idea of Daniel Bernstein http://cr.yp.to/factorization/smoothparts-20040510.pdf
 * Since it is based on long numbers it can only determine 37 - smooth numbers, 
 * because the product of the first primes has to fit into a long value. 
 * It calculates the primorial (the product of the primes) modulus the number 'x'
 * to be checked. Then this remainder is squared again modulus 'x' for a finite 
 * small time < 4 and checked weather this is '0'.
 * @author thiloharich
 *
 */
public class SmoothNumbers {

	final public static long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l;

	private int[] primes;
	private double[] reciprocals;
	
	// The allowed discriminator bit size is d <= 53 - bitLength(N/p), thus d<=23 would be safe
	// for any integer N and p>=2. d=10 is the value that performs best, determined by experiment.
	private static final double DISCRIMINATOR = 1.0/(1<<10);

	public SmoothNumbers(int ... primes) {
		this.primes = primes;
		reciprocals = new double [primes.length];
		for (int i = 0; i < primes.length; i++) {
			reciprocals[i] = 1.0 /primes[i];			
		}
	}
	/**
	 * If we have a table of all smooth numbers below n^1/3 we construct a semi smooth number
	 * x = s_1 * s_2 * f^2 where s_1,s_2, f^2 ~ n^1/3
	 * for all smooth s_2 and the given f^2 we calculate s_1 ~ n / (f^2 * s_2).
	 *
	 */

	/**
	 * Given a number n and a factor f, we are going to create a number s = z*f greater then n, where
	 * s <= n + n^1/4 * f^3/4 and s is smooth over a factor base. If no number can be found -1 will
	 * be returned.
	 *
	 * We look for z = (x+l) * (x-l) * f , where x = ceil(sqrt(n/f))
	 * x+l and x-l have to be smooth. so we try different l near the optimal solution l = sqrt(x^2*f - n).
	 *
	 * If we want to use such kind of s to be smooth over the given factor base,
	 * we want s - n ~ (n/f)^1/2  = size of the numbers we will sieve over
	 * n^1/4 * f^3/4 = (n/f)^1/2
	 * f^5/4 = n^1/4
	 * f = n^1/5
	 *
	 * By choosing f = n^1/5 the resulting number s and the numbers to be sieved on have the same size
	 * n^4/10 instead of n^1/2 in case of the QS.
	 *
	 * prob of a number lower x to be smooth over a factor base of size y is
	 * r^-r where r = log (x)/log(y)
	 * here x = n^4/10 instead of x = n^1/2 in the QS, but we have to combine 3 numbers of this size:
	 * -> r^-3r and r= 4/10 * log(n) /log y, we might choose y lower here but y growth superloynomial slower then x
	 *
	 * creating smooth numbers is done by sieving over numbers of size n' = sqrt(n/f).
	 * Thus this takes time exp(sqrt(log(sqrt(n'))*log(log(sqrt(n'))).
	 *
	 * We look for smooth numbers around sqrt(n/f) +/- floor(sqrt(n/f))^2 - n +/- searchInterval
	 * where searchInterval < exp(sqrt(log(sqrt(n'))*log(log(sqrt(n'))) = O(n^1/e) for each e.
	 * the size of the biggest factor of the factor base should be lower then the
	 * search interval to be able to search properly.
	 * The bigger the factor base will be chosen, the higher the chances that there will
	 * be a smooth number.
	 *
	 * We use the following fact: if x = ceil(sqrt(n/f)) then x^2*f - n = y <= 2 * sqrt(n/f)*f + 1
	 * <= 2 * sqrt(n*f)
	 * if l = sqrt(y/f) + e is an natural number then
	 * (x+l)(x-l)*f - n = x^2*f - l^2*f - n = y - y + 2*sqrt(y/f)*f*e + f*e^2 <= 2*sqrt(y*f)*e + f*e^2 <=
	 * 2*sqrt(sqrt(n*f)*f)*e  <= 2*n^(1/4)*f^3/4*e
	 * searchInterval < O((n*f)^(1/4+epsilon)) for any epsilon
	 * l^2 * f = y
	 * due to the chosen search interval and the size of the factor base, there should be such
	 * an l with a high probability fulfilling these conditions.
	 *
	 * @param n
	 * @param factorBase
	 * @param factor
	 * @return
	 */
	public long createSmooth(long n, int[] factorBase, long factor){
		double sqrt = sqrt(n / factor);
		long sqrtNDivF = (long) ceil(sqrt);
		long y = sqrtNDivF * sqrtNDivF * factor - n;
		double shiftD = Math.sqrt(y/(factor));
		long shift = Math.round(Math.sqrt(y/(factor)));
		double logSqrtN = log(sqrtNDivF); // = 1/2 * log(n)
//		int sieveInterval = (int) exp(sqrt(logSqrtN * log(logSqrtN))*.7);
		int sieveInterval = factorBase[factorBase.length -1];
		long xLower = sqrtNDivF - shift;
		long xHigher = sqrtNDivF + shift;
//		double xLowerD = (sqrtNDivF - shift);
//		double xHigherD = (sqrtNDivF + shift);
		double target = xLower * xHigher * factor - n;
//		double factorBaseSize = sieveInterval/2;

		int[] numberLength = new int [sieveInterval];
		long beginIntervalLower = xLower + (sieveInterval >> 1);
		long beginIntervalHigher = xHigher - (sieveInterval >> 1);
		// sieve the numbers around sqrt(n) - sqrt(y)
		for (int i = factorBase.length -1; i > 0; i--) {
			long beginPos = ((beginIntervalHigher / factorBase[i])+ 1) * factorBase[i];
			int primeLength = 31 - Integer.numberOfLeadingZeros(factorBase[i]);
			for (long pos = beginPos - beginIntervalHigher; pos < sieveInterval; pos += factorBase[i]){
				numberLength [(int) pos] -= primeLength;
				assertTrue((pos + beginIntervalHigher) % factorBase[i] == 0);
			}
		}
		final TDiv63Inverse tdiv = new TDiv63Inverse(factorBase[factorBase.length -1]);

		// find the smoothNumbers. Here we do real factoring.
		// TODO do resieving
		int lengthThreshold = 31 - Integer.numberOfLeadingZeros((int) sqrtNDivF);
		long minY = Long.MAX_VALUE;
		long minXLower = Long.MAX_VALUE;
		long minXHigher = Long.MAX_VALUE;
		SortedMultiset<BigInteger> minFactorXLower;
		SortedMultiset<BigInteger> minFactorXHigher;
		for (int i=0; i < sieveInterval; i++){
			int soothLength = -lengthThreshold + 5;
			if (numberLength[i] < soothLength){
//				smoothCandidates.add(i);
				long smoothXHigher = i + beginIntervalHigher;
				SortedMultiset<BigInteger> factorXHigher = tdiv.factor(BigInteger.valueOf(smoothXHigher));
				System.out.print("factors of x-i = " + factorXHigher);
				int product = 1;
				for (Map.Entry<BigInteger, Integer> entry:
						factorXHigher.entrySet()) {
					product *= Math.pow(entry.getKey().intValue(), entry.getValue());
				}
				if (product == smoothXHigher) {
					long smoothXLower = beginIntervalLower - i;
					SortedMultiset<BigInteger> factorXLower = tdiv.factor(BigInteger.valueOf(smoothXLower));
					System.out.print("factors of x+i = " + factorXLower);
					product = 1;
					for (Map.Entry<BigInteger, Integer> entry:
							factorXLower.entrySet()) {
						product *= Math.pow(entry.getKey().intValue(), entry.getValue());
					}
					if (product == smoothXLower) {
						long newY = smoothXLower * smoothXHigher * factor - n;
						if (Math.abs(newY) < minY){
							minY = Math.abs(newY);
//							minFactorXHigher = factorXHigher;
//							minFactorXLower = factorXLower;
							minXLower = smoothXLower;
							minXHigher = smoothXHigher;
						}
					}
				}
			}
		}
		if (minXLower != Long.MAX_VALUE){
			return minXLower * minXHigher * factor - n;
		}
		return -1;
	}

	/**
	 * from a square x^2 near n we create a smooth number over the factor base.
	 * We sieve around sqrt(n) for factors of size n^1/8. The search limit is at least the maximal factor
	 * from the database. If the smooth part of x is s. Then
	 * x^2 = s^2 * p^2 , n^(1/8 - e) > s > n^(1/8 + e) -> p <= n^(3/8+e)
	 * We find (efficiently!?) smooth p+i and p-i
	 * -> x_i = s^2 * (p-i) * (p+i)  is smooth
	 * since x_i - n = = s^2 * (p^2 - i^2) - n = (s*p)^2 - s^2*i^2 - n = x^2 - n - s^2*i^2
	 * (s*p)^2 - n - s^2*i^2 = 0 mod q
	 *
	 *
	 * if x = n^1/2 + l
	 *
	 * x_i - n = 2*l*n^1/2  - s^2*i^2
	 * i^2 = 2*l*n^1/2 / s^2
	 * i ~ 2*l^1/2 * n^1/4 / n^1/8 (1+e)
	 * i ~ 2*l^1/2 n^1/8 (1+e)
	 * x_{i+i} - x_i =  s^2*((i+1)^2 - i^2)
	 * = s^2*(2i+1)
	 * = n^1/4 * n^1/8 = n^3/8
	 *
	 * 1/2 - e = 2e + (1/4 - e)
	 * 1/2  = 2e + 1/4
	 * 1/4 = 2e
	 * e = 1/8
	 *
	 * This means (p-i), (p+i), x_i - n are all bounded by n^(3/8+e)
	 *
	 * Since p+i and y_i = x_i - n is a polynom in i we can technically sieve to find solutions.
	 * For any e' there is an n were the sieve interval is greater then the
	 * maximal factor from the factor base -> we can sieve efficiently.
	 * To find a x with a smooth part of size n^(1/4 - e) we can sieve with
	 * a factor base which is much smaller then the factor base of the regular QS.
	 * If x is smooth over the small factor base, there are many (smooth) p of size n^(1/4 + e).
	 * If e is bigger, there are many i and the probability to find an x_i near n is higher.
	 * For small n we can check if p +/- i is smooth by a precalculated lookup table.
	 * Finally we will check if y_i = x_i - n is smooth
	 * @param n
	 * @param factorBase
	 * @return
	 */
	public long createSmoothFromSquare(long n, int[] factorBase) {
		double sqrtN = sqrt(n);
		long sqrtNCeil = (long) Math.ceil(sqrtN);
//		long y = sqrtNCeil*sqrtNCeil - n;
//		long shift = (long) Math.sqrt(y);
		int sieveInterval = factorBase[factorBase.length -1]*4;
		double e = Math.log(sieveInterval)/Math.log(n);
		// e = 1/6 = 0.1666 is the theoretical limit if all p+i, p-i would be smooth
		e = Math.min(.05, e);
		// e = 1/10, smoothBound = n^0.375 + epsilon
		int smoothBound = (int) Math.pow(n, .5);
		long beginInterval = sqrtNCeil - (sieveInterval >> 1);
		int[] numberLengthSqrtN = new int [sieveInterval];
		int[] numberLength = new int [smoothBound];
		for (int i = factorBase.length -1; i >= 0; i--) {
			long beginPos = ((beginInterval / factorBase[i])+ 1) * factorBase[i];
			int primeLength = 32 - Integer.numberOfLeadingZeros(factorBase[i]-1);
			for (long pos = beginPos - beginInterval; pos < sieveInterval; pos += factorBase[i]){
				numberLengthSqrtN [(int) pos] -= primeLength;
				assertTrue((pos + beginInterval) % factorBase[i] == 0);
			}
			for (long pos = factorBase[i]; pos < smoothBound; pos += factorBase[i]){
				numberLength [(int) pos] -= primeLength;
				assertTrue((pos) % factorBase[i] == 0);
			}
		}
		int lengthThreshold;
		SortedMultiset<BigInteger> minFactorXLower;
		BitSet smoothNumber = new BitSet(smoothBound);
		final TDiv63Inverse tdiv = new TDiv63Inverse(factorBase[factorBase.length -1]);
		for (int j=smoothBound-1; j > 0 ; j--){
			lengthThreshold = 32 - Integer.numberOfLeadingZeros(j-1);
			int soothLength = -lengthThreshold + 2;
			if (numberLength[j] <= soothLength) {
				smoothNumber.set(j);
			}
		}
		int nPowOneFourth = (int) Math.pow(n, .125);
		long xi = 0;
		int lengthThresholdX = 32 - Integer.numberOfLeadingZeros(nPowOneFourth-1);
		for (int j=sieveInterval / 2; j < sieveInterval; j++) {
			int soothLength = -lengthThresholdX + 3;
			if (numberLengthSqrtN[j] < soothLength) {
				// here we found j dividing x
				long x = j + beginInterval;
				SortedMultiset<BigInteger> factorX = tdiv.factor(BigInteger.valueOf(x));
				List<Integer> factorList = new ArrayList<>();
//				long s = 1;
				for (Map.Entry<BigInteger, Integer> entry : factorX.entrySet()) {
					for (int i = 0; i < entry.getValue(); i++) {
						factorList.add(entry.getKey().intValue());
//						s *= Math.pow(entry.getKey().intValue(), entry.getValue());
					}
				}
				Set<Integer> smooth = new HashSet<>();
				// optimal value for lower factor s is n^1/8 = n^0.125
				int lower = (int) Math.pow(n, .1);
				int higher = (int) Math.pow(n, .15);
				findFactors(lower, higher, factorList, 1, smooth, 0);
				System.out.println(factorX);
				for (int s : smooth) {
					long xDivS = (int) (x / s);
					long y0 = (s * s * xDivS * xDivS - n) / (s*s);
					long iOpt = (long) Math.sqrt(y0);
//					y += s^2*((iOpt + i2)^2 - iOpt^2)
//					smoothBound - y0 = s^2*(2*iOpt * i2 + i2^2)
//					(smoothBound - y0)/(s^2 * 2*iOpt) ~ i2
					long iRange = (smoothBound -y0)/ (s*s * 2 * iOpt) + 1;
					iRange = iRange < 0 ? 1 : iRange;
					long i = iOpt- iRange;
					boolean found = false;
					while (i < iOpt + iRange && i > 0) {
						i = xDivS - smoothNumber.previousSetBit((int) (xDivS - i));
						if (smoothNumber.get((int) (xDivS + i))) {
							found = true;
							long x2 = s * s * (xDivS) * (xDivS);
							xi = s * s * (xDivS - i) * (xDivS + i);
							int y = Math.abs((int) (xi - n));
							SortedMultiset<BigInteger> factorsXi = tdiv.factor(BigInteger.valueOf(xi));
							SortedMultiset<BigInteger> factorsY = tdiv.factor(BigInteger.valueOf(y));
							System.out.print("Found candidate " + i + " fact y " + factorsY);
							System.out.println();
//							if (Math.abs(y) < 2 * dist * sqrtN) {
							if (smoothNumber.get(Math.abs(y)))
							{
								System.out.println("s : " + s + " y : " + y);
							}
//							}
						}
//						else {
							i = smoothNumber.nextSetBit((int) (xDivS + i+1)) - xDivS;
//						}
					}
				}
			}
		}
		return xi;
	}

		/**
         * we sieve around x = sqrt(n) to find some factors i ~ n^1/4 which divides x.
         * This is x=i*x'
         * if i consists out of many small primes we can build many i' which divide i and x
         * If (x-i') and (x+i') are smooth
         * (x-i')*(x+i') are smooth and (x+i')(x-i') - n = x^2 - i'^2 - n < x^2 -n
         * since x-i' are divided by i' x/i' -1 is of size n^1/2 / n^1/4 = n^1/4
         * for i = sqrt(x^2 - n) + e we have
         * (x+i)(x-i) - n = x^2 - i^2 - n = x^2 - n - (x^2 - n) + 2*sqrt(x^2 - n)*e + e^2
         * = 2*sqrt(x^2 - n)*e + e^2 = O(n^1/4 * e)
         *
         * (n^1/2+i-l)(n^1/2+i+l) - n =
         * (n^1/2+i)^2-l^2 - n = n +2i*n^1/2 + i^2 - l^2 - n
         * we need 2i*n^1/2 = l^2
         * l < i^1/2 * n^1/4
         * @param n
         * @param factorBase
         * @return
         */
	public long createSmooth(long n, int[] factorBase){
		double sqrtN = sqrt(n);
		long sqrtNCeil = (long) Math.ceil(sqrtN);
//		long y = sqrtNCeil*sqrtNCeil - n;
//		long shift = (long) Math.sqrt(y);
		int sieveInterval = factorBase[factorBase.length -1]*4;
		long beginInterval = sqrtNCeil - (sieveInterval >> 1);
		int[] numberLengthSqrtN = new int [sieveInterval];
		int[] numberLength = new int [sieveInterval];
		for (int i = factorBase.length -1; i > 0; i--) {
			long beginPos = ((beginInterval / factorBase[i])+ 1) * factorBase[i];
			int primeLength = 31 - Integer.numberOfLeadingZeros(factorBase[i]);
			for (long pos = beginPos - beginInterval; pos < sieveInterval; pos += factorBase[i]){
				numberLengthSqrtN [(int) pos] -= primeLength;
				assertTrue((pos + beginInterval) % factorBase[i] == 0);
			}
			for (long pos = factorBase[i]; pos < sieveInterval; pos += factorBase[i]){
				numberLength [(int) pos] -= primeLength;
				assertTrue((pos) % factorBase[i] == 0);
			}
		}
		int nPowOneFourth = (int) Math.pow(n, 1.0/3);
		int lengthThreshold;
		SortedMultiset<BigInteger> minFactorXLower;
		BitSet smoothNumber = new BitSet(sieveInterval);
		final TDiv63Inverse tdiv = new TDiv63Inverse(factorBase[factorBase.length -1]);
		for (int j=0; j < sieveInterval; j++){
			lengthThreshold = 31 - Integer.numberOfLeadingZeros(j);
			int soothLength = -lengthThreshold + 3;
			if (numberLength[j] < soothLength) {
				smoothNumber.set(j);
			}
		}
		lengthThreshold = 31 - Integer.numberOfLeadingZeros(nPowOneFourth);
		for (int j=sieveInterval / 2; j < sieveInterval; j++){
			int soothLength = -lengthThreshold + 3;
			if (numberLengthSqrtN[j] < soothLength) {
				// here we found j dividing x
				long x = j + beginInterval;
				SortedMultiset<BigInteger> factorX = tdiv.factor(BigInteger.valueOf(x));
				List<Integer> factorList = new ArrayList<>();
				for (Map.Entry<BigInteger, Integer> entry: factorX.entrySet()) {
					for (int i=0; i < entry.getValue(); i++)
						factorList.add(entry.getKey().intValue());
				}
				int y = (int) (x*x - n);
				double target = sqrt(y);
				long dist = x - sqrtNCeil;
				System.out.println("dist " + dist + " target : " + target + "factors of x = " + factorX);
				int yMin = (int) ((x - (long)target)*(x + (long)target) - n);
				Set<Integer> shifts = new HashSet<>();
				findFactors(target, factorList, 1, shifts, 1);
				for (int i: shifts) {
					int i2 = i;
					while (i < 2*target) {
						int xDivi = (int) (x / i);
						if (smoothNumber.get(xDivi - 1) && smoothNumber.get(xDivi + 1)) {
							y = (int) ((x - i) * (x + i) - n);
//							if (Math.abs(y) < 2 * dist * sqrtN) {
								if (smoothNumber.get(Math.abs(y))) {
									System.out.print("Found number shift " + i + " rel size " + y / sqrtN);
									System.out.println();
								}
//							}
						}
						i += i2;
					}
				}
			}
		}
		return -1;
	}

	private void findFactors(double lower, double higher, List<Integer> factorX, int prod, Set<Integer> shifts, int index) {
		if (lower < prod && prod < higher) {
			shifts.add(prod);
		}
		if (index == factorX.size())
			return;
		findFactors(lower, higher, factorX, prod, shifts, index + 1);
		findFactors(lower, higher, factorX, prod * factorX.get(index), shifts, index + 1);
	}
	private void findFactors(double target, List<Integer> factorX, int prod, Set<Integer> shifts, int index) {
		if (target / prod > .5 && target / prod < 200.0) {
			shifts.add(prod);
		}
		if (index == factorX.size())
			return;
		findFactors(target, factorX, prod, shifts, index + 1);
		findFactors(target, factorX, prod * factorX.get(index), shifts, index + 1);
	}

	/**
	 * Does not detect all smooth numbers. Uses 3 slow modulus operations.
	 * @param x - the number to be checked to be smooth over the primorial
	 * @param primorial - the product of the prime factors which should be in the prime factorization of x
	 * @return - true if the number is smooth. Might return false, even if the number is smooth.
	 */
	boolean isSmoothBy3Modulus(long x, long primorial){
		
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		long y = primorial % x;
		y = (y * y) % x;
		y = (y * y) % x;  // we loose 13% smooth numbers here for numbers ~ 2^29
//		y = (y * y) % x;  // we loose  6% smooth numbers here
//		y = (y * y) % x;  // perfect match
		return y == 0;
	}
	/**
	 * Does not detect all smooth numbers. Calculate 1/x only once and
	 * calculates the modulus values by subtractions and multiplications.
	 * @param x - the number to be checked to be smooth over the primorial
	 * @param primorial - the product of the prime factors which should be in the prime factorization of x
	 * @return - true if the number is smooth. Might return false, even if the number is smooth.
	 */
	public boolean isSmoothBy3Inv(long x, long primorial){

		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		double xInv = 1.0 / x;
		long y = PrimeMath.mod(primorial, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		long rem = PrimeMath.gcd(y, x);
		return y == 0;
	}
	public static long remainder(long x, long primorial){
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		double xInv = 1.0 / x;
		long y = PrimeMath.mod(primorial, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		y = PrimeMath.mod(y * y, x, xInv);
		long rem = PrimeMath.gcd(y, x);
		return rem << exp2;
	}
	/**
	 * @param x - the number to be checked to be smooth over the primorial
	 * @param primorial - the product of the prime factors which should be in the prime factorization of x
	 * @return
	 */
	boolean isSmoothBy5Modulus(long x, long primorial){
		
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		long y = primorial % x;
		y = (y * y) % x;
		y = (y * y) % x;  // we loose 13% smooth numbers here for numbers ~ 2^29
		y = (y * y) % x;  // we loose  6% smooth numbers here
		y = (y * y) % x;  // perfect match
		return y == 0;
	}
	
	/**
	 * 
	 * @param x
	 * @return a positive number if x can be factorized.
	 */
	long factorizeIfSmooth (long x) {
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		long singlePrimes = exp2;
		
		for (int i = 1; i < primes.length; i++) {
			final long xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
			if (xDivPrime * primes[i] == x) {
				singlePrimes += i;
				singlePrimes <<= 4;
				// TODO multiply the primes
				x = xDivPrime;
			}
		}
		long primeFactors = 0;
		int iLast = 13;
		singlePrimes >>= 4;
		int i = 0;
		
		while(x > 1 && (singlePrimes & 15) != 0) {
//			while((singlePrimes & 15) != 0) {
			i = (int) (singlePrimes & 15);
			singlePrimes >>= 4;
			primeFactors <<= 5*(iLast - i);
			primeFactors++;
			iLast = i;
			long xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
			while (xDivPrime * primes[i] == x) {
				x = xDivPrime;
				xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
				primeFactors++;
			}
			singlePrimes >>= 4;			
		}
		primeFactors <<= 5*i;

		primeFactors += exp2;
		primeFactors = x == 1 ? primeFactors : -primeFactors;
		return primeFactors;
	}
	
	/**
	 * 
	 * @param x
	 * @param maxFactorIndex
	 * @return a positive number if x can be factorized. 
	 */
	public int [] factorizeIfSmooth1 (long x, int maxFactorIndex) {
		int [] primeFactors = new int [maxFactorIndex+2]; 
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		primeFactors[0] = exp2;
		
		for (int i = 1; i < maxFactorIndex; i++) {
			while ((long) (x*reciprocals[i] + DISCRIMINATOR) * primes[i] == x) {
				primeFactors[i]++;
				x = (long) (x*reciprocals[i] + DISCRIMINATOR);
			}
//			while ((long) (x*reciprocals[i]) * primes[i] == x) {
//				primeFactors[i]++;
//				x = (long) (x*reciprocals[i]);
//			}
		}
		primeFactors[maxFactorIndex+1] = (int) x;
		return primeFactors;
	}
	
	/**
	 * 
	 * @param x
	 * @param maxFactor
	 * @return a positive number if x can be factorized. 
	 */
	int [] factorizeIfSmooth2 (long x, long maxFactor) {
		int [] primeFactorIndex = new int [primes.length]; 
		int [] primeFactors     = new int [primes.length+1]; 
		int exp2 = Long.numberOfTrailingZeros(x);
		x = x >> exp2;
		primeFactors[0] = exp2;
		int endPos = 0;
		
		for (int i = 1; i < primes.length; i++) {
			long xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
			if (xDivPrime * primes[i] == x) {
				primeFactorIndex[endPos++] = i;
				primeFactors[i] = 1;
				x = xDivPrime;
				xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
			}
		}
		int primeIndex = 0;
		// TODO check for < 37 and stop early.
		while(x > 1 && primeIndex < endPos) {
			final int i = primeFactorIndex[primeIndex++];
			long xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
			while (xDivPrime * primes[i] == x) {
				x = xDivPrime;
				xDivPrime = (long) (x*reciprocals[i] + DISCRIMINATOR);
				primeFactors[i]++;
			}		
		}
		primeFactors[primes.length] = x == 1 ? 1 : 0;
		return primeFactors;
	}
	
}
