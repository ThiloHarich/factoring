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
	public LehmanSingleFactorFinder(int bits) {
		this.bits = bits;
//		impl = new LehmanNoSqrtFact(bits, 1.01f);
		impl = new LehmanReverseFact(bits, 1.5f);
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		//        LehmanNoSqrtFact impl = new LehmanNoSqrtFact(41, 1.001f);
		final long factor = impl.findPrimeFactors(N.longValue(), null);
		return BigInteger.valueOf(factor);
	}

	@Override
	public SortedMultiset<BigInteger> factor(BigInteger N) {
		TreeMultiset<Long> primeFactors = TreeMultiset.create();
		final long factor = impl.findPrimeFactors(N.longValue(), primeFactors);
		if (factor != 1)
			primeFactors.add(factor);
		SortedMultiset result = new SortedMultiset_BottomUp(primeFactors);
		return result;
	}

	@Override
	public String getName() {
		return LehmanSingleFactorFinder.class.getSimpleName();
	}
}
