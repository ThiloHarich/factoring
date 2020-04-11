package factoring.fermat;

import java.math.BigInteger;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;
import de.tilman_neumann.jml.factor.FactorAlgorithm;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class FermatXNonSquare extends FactorAlgorithm implements FactorizationOfLongs{



	private long minFactor = 3;

	@Override
	public long findFactors(long n, Collection<Long> factors) {
		final long sqrtN = (long) Math.ceil(Math.sqrt(n));
		final long xEnd = (n / minFactor  + minFactor) / 2;
		for (long x = sqrtN; x <= xEnd; x++) {
			for(long noSquare = 1600 ; noSquare >= 1; noSquare--) {
				final long right = noSquare * x*x - n;
				if (PrimeMath.isSquare(right)) {
					final long y = PrimeMath.sqrt(right);
					final long factorHigh = x + y;
					factors.add(factorHigh);
					return x - y;
				}
			}
		}
		return n;
	}

	@Override
	public String getName() {
		return null;
	}

	@Override
	public BigInteger findSingleFactor(BigInteger N) {
		Set<Long> factors = new HashSet<>();
		return BigInteger.valueOf(findFactors(N.longValue(), factors));
	}
}
