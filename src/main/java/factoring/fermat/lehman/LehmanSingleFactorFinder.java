package factoring.fermat.lehman;

import java.math.BigInteger;

import com.google.common.collect.TreeMultiset;
import de.tilman_neumann.math.factor.SingleFactorFinder;
import de.tilman_neumann.types.SortedMultiset;
import de.tilman_neumann.types.SortedMultiset_BottomUp;

/**
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanSingleFactorFinder extends AbstractFactorFinder implements SingleFactorFinder, PrimeFactorFinder {
//	private final LehmanNoSqrtFact impl;
	int bits;
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

	@Override
	public String getName() {
		return LehmanSingleFactorFinder.class.getSimpleName();
	}
}
