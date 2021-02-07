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
package factoring.trial;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.primes.exact.AutoExpandingPrimesArray;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.SortedMultiset_BottomUp;
import org.apache.log4j.Logger;

import java.math.BigInteger;
import java.util.BitSet;

import static de.tilman_neumann.jml.base.BigIntConstants.I_2;

/**
 * Trial division using Barrett reduction,
 * see https://en.wikipedia.org/wiki/Barrett_reduction.
 *
 * Significantly faster than TDiv31Inverse.
 *
 * @authors Tilman Neumann + Thilo Harich
 */
public class TDiv31Barrett extends FactorAlgorithm {
	@SuppressWarnings("unused")
	private static final Logger LOG = Logger.getLogger(TDiv31Barrett.class);

	private AutoExpandingPrimesArray SMALL_PRIMES = AutoExpandingPrimesArray.get();	// "static" would be slightly slower

	private int[] primes;
	private long[] pinv;

	public TDiv31Barrett() {
		primes = new int[NUM_PRIMES_FOR_31_BIT_TDIV];
		pinv = new long[NUM_PRIMES_FOR_31_BIT_TDIV];
		for (int i=0; i<NUM_PRIMES_FOR_31_BIT_TDIV; i++) {
			int p = SMALL_PRIMES.getPrime(i);
			primes[i] = p;
			pinv[i] = (1L<<32)/p;
		}
	}

	@Override
	public String getName() {
		return "TDiv31Barrett";
	}

	@Override
	public SortedMultiset<BigInteger> factor(BigInteger Nbig) {
		SortedMultiset<BigInteger> primeFactors = new SortedMultiset_BottomUp<>();
		int N = Nbig.intValue();

		// Powers of 2 can be removed very fast.
		// This is required also because the Barrett division does not work with p=2.
		int lsb = Integer.numberOfTrailingZeros(N);
		if (lsb > 0) {
			primeFactors.add(I_2, lsb);
			N >>= lsb;
		}

		// Test odd primes
		int q;
		for (int i=1; ; i++) {
			final long r = pinv[i];
			final int p = primes[i];
			int exp = 0;
			while ((q = (1 + (int) ((N*r)>>32))) * p == N) {
				exp++;
				N = q;
			}
			if (exp>0) {
				primeFactors.add(BigInteger.valueOf(p), exp);
			}
			if (p*(long)p > N) {
				break;
			}
		}

		if (N>1) {
			// either N is prime, or we could not find all factors with p<=pLimit -> add the rest to the result
			primeFactors.add(BigInteger.valueOf(N));
		}
		return primeFactors;
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		if (N.bitLength() > 31) throw new IllegalArgumentException("TDiv31Barrett.findSingleFactor() does not work for N>31 bit, but N=" + N);
		return BigInteger.valueOf(findSingleFactor(N.intValue(), NUM_PRIMES_FOR_31_BIT_TDIV));
	}

	public boolean factor(int N, int maxIndex, SortedMultiset<BigInteger> primeFactors) {
//		SortedMultiset<BigInteger> primeFactors = new SortedMultiset_BottomUp<>();
		if (N<0) N = -N; // sign does not matter

		int lsb = Integer.numberOfTrailingZeros(N);
		if (lsb > 0) {
			primeFactors.add(I_2, lsb);
			N = N >> lsb;
		}
		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int i=0;
		int unrolledLimit = maxIndex;
		// TODO we might check if primes[i]^2 < N
		for ( ; i<unrolledLimit; i++) {
			while ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) {
				N = (1 + (int) ((N*pinv[i])>>32));
				primeFactors.add(BigInteger.valueOf(primes[i]));
			};
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
//			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
		}
//		for ( ; i<maxIndex; i++) {
//			if ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) return primes[i];
//		}
//		if (N <= primes[maxIndex-1] && N > 1){
//			primeFactors.add(BigInteger.valueOf(N));
//			return true;
//		}
		return N == 1;
	}

	public int findSingleFactor(int N, int maxIndex) {
		if (N<0) N = -N; // sign does not matter
		if ((N&1)==0) return 2; // N even

		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int i=1;
		int unrolledLimit = maxIndex-8;
		for ( ; i<unrolledLimit; i++) {
			if ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return primes[i];
		}
		for ( ; i<maxIndex; i++) {
			if ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) return primes[i];
		}
		// otherwise N is prime
		return 1;
	}
	public int findSingleFactorIndex(int N, int maxIndex) {
		if (N<0) N = -N; // sign does not matter
		if ((N&1)==0) return 0; // N even

		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int i=1;
		int unrolledLimit = maxIndex-8;
		for ( ; i<unrolledLimit; i++) {
			if ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
			if ((1 + (int) ((N*pinv[++i])>>32)) * primes[i] == N) return i;
		}
		for ( ; i<maxIndex; i++) {
			if ((1 + (int) ((N*pinv[i])>>32)) * primes[i] == N) return i;
		}
		// otherwise N is prime
		return -1;
	}

	/**
	 * Here we calculate a representation of smooth
	 * for the matrix step, which is using very very few memory.
	 * This consinsts out of three parts:
	 * <li>a bit set of the (odd) exponents below factor index splitIndex</li>
	 * <li> the highest factor</li><br>
	 * The first set is returned by the method and is essential for finding the solution of the matrix.
	 * The second is needed to find the factor after solving the maxrix.
	 * In oder to calculate the sqrt(x) and sqrt(y) for a relation  x^2 - n = y^2.
	 * It is passed over as parameter factorCounts.
	 * The representation allows fast multiplications of the representations.
	 * It is represented by an array of longs, where some (8) bits per factor represents
	 * the exponent of a factor of the prime representation.
	 *
	 * If the number is smooth over the factor bas up to factor with index maxIndex .
	 *
	 * @param factorCounts a bit set of the exponents mod 2
	 * @param N
	 * @param maxIndex
	 * @return
	 */
	public ExponentsMaxFactor fillFactorCounts(int N, int maxIndex) {
		BitSet factors = new BitSet();
		if (N<0) N = -N; // sign does not matter

		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int number = N;
		int maxFactorIndex = 0;

		long exponent = 0l;
		while ((number & 1) == 0){
			number >>= 1;
			exponent++;
		}
		if ((exponent & 1l) == 1l ) {
			factors.set(1);
			maxFactorIndex = 1;
		}
		int i=1;
		for ( ; i<maxIndex && number > 1; i++) {
			exponent = 0;
			while ((1 + (int) ((number*pinv[i])>>32)) * primes[i] == number) {
				number = 1 +  (int) ((number * pinv[i]) >> 32);
				exponent++;
			}
			int primeIndex = i+1;
			if ((exponent & 1l) == 1l )
				factors.set(primeIndex);
			maxFactorIndex = primeIndex;
			// 8 = 2^3 bits per exponent -> epxonent per prime < 256 -> number > p^256 >= 2^256
			// TODO dynamic bits ?
//			final int longIndex = primeIndex >> 3;
//			final int indexInWord = (primeIndex - (longIndex << 3)) << 3;
//			factorCounts[longIndex] |= (exponent << indexInWord);
		}
		// if the number can not be factored (number > 1) adjust max factor index
		maxFactorIndex = number == 1 ? maxFactorIndex : Short.MAX_VALUE;
		// otherwise N is prime
		ExponentsMaxFactor result = new ExponentsMaxFactor(factors, maxFactorIndex, number == 1);
		return result;
	}
	public boolean fillFactorCounts(int N, int maxIndex, long [] factorCounts) {
		BitSet factors = new BitSet();
		if (N<0){
			N = -N; // sign does not matter
			factorCounts[0] += 1;
		}
		int number = N;

		long exponent = 0l;
		while ((number & 1) == 0){
			number >>= 1;
			exponent++;
		}
		factorCounts[0] += 256 * exponent;
		int i=1;
		for ( ; i<maxIndex-1 && number > 1; i++) {
			exponent = 0;
			while ((1 + (int) ((number*pinv[i])>>32)) * primes[i] == number) {
				number = 1 +  (int) ((number * pinv[i]) >> 32);
				exponent++;
			}
			int primeIndex = i+1;
			// 8 = 2^3 bits per exponent -> epxonent per prime < 256 -> number > p^256 >= 2^256
			// TODO dynamic bits ?
			final int longIndex = primeIndex >> 3;
//			final int indexInWord = primeIndex << 3;
			final int indexInWord = (primeIndex - (longIndex << 3)) << 3;
			factorCounts[longIndex] += (exponent << indexInWord);
		}
		// otherwise N is prime
		return number == 1;
	}

	public class ExponentsMaxFactor {
		public BitSet exponents;
		public int maxFactorIndex;
		public boolean isSmooth;

		public ExponentsMaxFactor(BitSet factors, int maxFactorIndex, boolean smooth) {
			this.exponents = factors;
			this.maxFactorIndex = maxFactorIndex;
			this.isSmooth = smooth;
		}
	}
}
