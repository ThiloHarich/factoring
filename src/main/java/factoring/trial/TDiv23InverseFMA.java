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

import java.math.BigInteger;

/**
 * Trial division factor algorithm replacing division by multiplications.
 * 
 * Instead of dividing N by consecutive primes, we store the reciprocals of those primes, too,
 * and multiply N by those reciprocals. Only if such a result is near to an integer we need
 * to do a division.
 * 
 * This variant abstains from testing N%primes[i] when the discriminator test indicates a neat division,
 * and unrolls the loop in findSingleFactor().
 * 
 * Another bit faster than TDiv31Inverse_NoDoubleCheck.
 * 
 * @authors Thilo Harich + Tilman Neumann
 */
public class TDiv23InverseFMA extends FactorAlgorithm {

	// TODO "static" leads to NPE in constructor ???
	private AutoExpandingPrimesArray SMALL_PRIMES = AutoExpandingPrimesArray.get();

	// The allowed discriminator bit size is d <= 53 - bitLength(N/p), thus d<=23 would be safe
	// for any integer N and p>=2. d=10 is the value that performs best, determined by experiment.
	private static final float DISCRIMINATOR = 0.01f;

	static int NUM_PRIMES = NUM_PRIMES_FOR_31_BIT_TDIV;

	private int[] primes;
	private float[] reciprocals;

	private static final double ROUNDING_CORRECTION = 0.01;
	private int maxFactor = 65535;
	float[] primesInv;
//	int[] primes;


	public TDiv23InverseFMA() {
		primes = new int[NUM_PRIMES];
		reciprocals = new float[NUM_PRIMES];
		for (int i=0; i<NUM_PRIMES; i++) {
			int p = SMALL_PRIMES.getPrime(i);
			primes[i] = p;
			reciprocals[i] = 1.0f/p;
		}
	}
	
	@Override
	public String getName() {
		return "TDiv23InverseFMA";
	}

	@Override
	public SortedMultiset<BigInteger> factor(BigInteger Nbig) {
		SortedMultiset<BigInteger> primeFactors = new SortedMultiset_BottomUp<>();
		int N = Nbig.intValue();
		
		int q;
		for (int i=0; ; i++) {
			double r = reciprocals[i];
			int p = primes[i];
			int exp = 0;
			while ((q = (int) (N*r + DISCRIMINATOR)) * p == N) {
				exp++;
				N = q; // avoiding a division here by storing q benefits the int version but not the long version
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
	public int findSingleFactor(int N) {
		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int i=0;
		if ((int) (N*reciprocals[i]   + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];
		if ((int) (N*reciprocals[++i] + DISCRIMINATOR) * primes[i] == N) return primes[i];

		// otherwise N is prime!
		throw new IllegalArgumentException("N = " + N + " is prime!");
	}
	public int findSingleFactor(int N, int maxIndex) {
		if (N<0) N = -N; // sign does not matter
		if ((N&1)==0) return 2; // N even

		// if N is odd and composite then the loop runs maximally up to prime = floor(sqrt(N))
		// unroll the loop
		int i=1;
		int unrolledLimit = maxIndex-4;
//		double prodOfMods = 1;
//		boolean found = false;
		for ( ; i<unrolledLimit; i++) {
//			prodOfMods =  Math.fma((int) (Math.fma(N,reciprocals[i], DISCRIMINATOR) ), primes[i],-N) *
//					Math.fma((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)), primes[i], -N) *
//					Math.fma((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)), primes[i], -N) *
//					Math.fma((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)), primes[i], -N);
//			found = ((int) (Math.fma(N,reciprocals[i], DISCRIMINATOR) )* primes[i] == N) ||
//			((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) ||
//			((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) ||
//			((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N);
			if ((int) (Math.fma(N,reciprocals[i], DISCRIMINATOR) )* primes[i] == N) return primes[i];
			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
//			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
//			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
//			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];
//			if ((int) (Math.fma(N,reciprocals[++i], DISCRIMINATOR)) * primes[i] == N) return primes[i];

		}
//		if (found){
////			if (prodOfMods == 0){
//			maxIndex -=  4;
//		}
		for ( ; i<maxIndex; i++) {
			if ((int) (N*reciprocals[i]   + DISCRIMINATOR) * primes[i] == N) return primes[i];
		}
		// otherwise N is prime
		return 1;
	}
}
