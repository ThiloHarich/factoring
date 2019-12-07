package factoring.smooth;

import java.util.Iterator;

import factoring.math.PrimeMath;

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

	final long PRIMORIAL = 2l*3l*5l*7l*11l*13l*17l*19l*23l*29l*31l*37l; 

	private int[] primes;
	private double[] reciprocals;
	
	// The allowed discriminator bit size is d <= 53 - bitLength(N/p), thus d<=23 would be safe
	// for any integer N and p>=2. d=10 is the value that performs best, determined by experiment.
	private static final double DISCRIMINATOR = 1.0/(1<<10);

	public SmoothNumbers(int[] primes) {
		this.primes = primes;
		reciprocals = new double [primes.length];
		for (int i = 0; i < primes.length; i++) {
			reciprocals[i] = 1.0 /primes[i];			
		}
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
	 * @param maxFactor
	 * @return a positive number if x can be factorized. 
	 */
	long factorizeIfSmooth (long x, long maxFactor) {
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
	 * @param maxFactor
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
