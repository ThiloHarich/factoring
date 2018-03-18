package factoring.fermat.lehman;

import de.tilman_neumann.math.factor.SingleFactorFinder;
import de.tilman_neumann.types.SortedMultiset;

import java.math.BigInteger;

/**
 * Created by Thilo Harich on 18.03.2018.
 */
public class LehmanSingleFactorFinder implements SingleFactorFinder {
    @Override
    public BigInteger findSingleFactor(BigInteger N) {
        LehmanNoSqrtFact impl = new LehmanNoSqrtFact(41, 1.001f);
        final long factor = impl.findPrimeFactors(N.longValue(), null);
        return BigInteger.valueOf(factor);
    }

    @Override
    public SortedMultiset<BigInteger> factor(BigInteger N) {
        // TODO
        return null;
    }

    @Override
    public String getName() {
        return LehmanSingleFactorFinder.class.getSimpleName();
    }
}
