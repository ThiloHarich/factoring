/*
 * PSIQS 4.0 is a Java library for integer factorization, including a parallel self-initializing quadratic sieve (SIQS).
 * Copyright (C) 2018  Tilman Neumann (www.tilman-neumann.de)
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
package factoring.squfof;

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.math.base.bigint.primes.exact.SieveFacade;
import de.tilman_neumann.util.ConfigUtil;

import static de.tilman_neumann.math.base.bigint.BigIntConstants.*;

/**
 * Sequence of multiplier * {squarefree numbers 1,2,3,5,6,7,10,11,13,...}, BigInteger implementation.
 * @author Tilman Neumann
 */
public class SquarefreeSequence implements IntegerSequence<BigInteger> {
	private static final Logger LOG = Logger.getLogger(SquarefreeSequence.class);

	private SieveFacade primeGen = SieveFacade.get();
	
	private BigInteger multiplier;
	private BigInteger next;
	
	public SquarefreeSequence(long multiplier) {
		this.multiplier = BigInteger.valueOf(multiplier);
	}

	public SquarefreeSequence(BigInteger multiplier) {
		this.multiplier = multiplier;
	}
	
	@Override
	public String getName() {
		return multiplier + "*squarefree";
	}

	public void reset() {
		this.next = ONE;
	}
	
	@Override
	public void reset(BigInteger N) {
		this.next = ONE;
	}

	@Override
	public BigInteger next() {
		BigInteger ret = next;
		// compute next square free number
		while (true) {
			next = next.add(ONE);
			boolean isSquareFree = true;
			// next must not be divisible by any prime square p^2 <= next.
			// divide by p and if necessary by another p is much faster than division by p^2 in one go!
			BigInteger test = next;
			int primeIndex = 0;
			for (BigInteger p = TWO; p.multiply(p).compareTo(test)<=0; p = BigInteger.valueOf(primeGen.getPrime(++primeIndex))) {
				// test p
				BigInteger[] div = test.divideAndRemainder(p);
				if (div[1].equals(ZERO)) {
					test = div[0];
					BigInteger[] div2 = test.divideAndRemainder(p);
					if (div2[1].equals(ZERO)) {
						// next is not square free !
						isSquareFree = false;
						break;
					}
				}
			}
			if (isSquareFree) break; // found next square-free number
		}
		return ret.multiply(multiplier);
	}
	
	// standalone test
	public static void main(String[] args) {
	   	ConfigUtil.initProject();
	   	SquarefreeSequence seqGen = new SquarefreeSequence(ONE);
	   	long start = System.currentTimeMillis();
		seqGen.reset(ONE); // >1 required!
		for (int i=1; i<=1000000; i++) {
			@SuppressWarnings("unused")
			BigInteger squarefree = seqGen.next();
			//LOG.info("squarefree(" + i + ") = " + squarefree);
		}
		LOG.info("computation took " + (System.currentTimeMillis()-start) + " ms");
	}
}
