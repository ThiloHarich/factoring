package factoring.fermat.lehman;

import java.math.BigInteger;

import com.google.common.collect.TreeMultiset;
import de.tilman_neumann.math.factor.SingleFactorFinder;
import de.tilman_neumann.types.SortedMultiset;
import de.tilman_neumann.types.SortedMultiset_BottomUp;

/**
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanSingleFactorFinder implements SingleFactorFinder {
//	private final LehmanNoSqrtFact impl;
	int bits;
	LehmanReverseFact impl;
//LehmanNoSqrtFact impl;
	public LehmanSingleFactorFinder(int bits) {
		this.bits = bits;
//		impl = new LehmanNoSqrtFact(bits, 1.01f);
		impl = new LehmanReverseFact(bits, 3f);
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		//        LehmanNoSqrtFact impl = new LehmanNoSqrtFact(41, 1.001f);
		final long factor = impl.findPrimeFactors(N.longValue(), null);
		return BigInteger.valueOf(factor);
	}

	@Override
	public SortedMultiset<BigInteger> factor(BigInteger n) {
		TreeMultiset<Long> allFactors = factorLong(n.longValue());
//		allFactors.addAll(primeFactors);
		SortedMultiset result = new SortedMultiset_BottomUp(allFactors);
		return result;
	}

	public TreeMultiset<Long> factorLong(long n) {
		// if we have a prime return an empty set
		TreeMultiset<Long> factorsEven = TreeMultiset.create();
		while ((n & 1) == 0)
		{
			factorsEven.add(2l);
			n = n >> 1;
		}
		if (n == 1) {
			return factorsEven;
		}
		TreeMultiset<Long> primeFactors = TreeMultiset.create();
		// find one factor and decomposite this factor and n/factor
		long factor1 = impl.findPrimeFactors(n, primeFactors);
		// if we do not find a divisor just return it
		if (factor1 == n){
			factorsEven.add(n);
			return factorsEven;
		}
		// also divide out the prime factorsEven
		long factor2 = n/factor1;
		for (long factor : primeFactors) {
			factor2 /= factor;
		}
		TreeMultiset<Long> subFactors1 = factorLong(factor1);
		TreeMultiset<Long> subFactors2 = factorLong(factor2);
		factorsEven.addAll(subFactors1);
		factorsEven.addAll(subFactors2);
		factorsEven.addAll(primeFactors);
		return factorsEven;
	}

	@Override
	public String getName() {
		return LehmanSingleFactorFinder.class.getSimpleName();
	}
}
