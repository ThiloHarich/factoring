package factoring.trial;

import java.util.Collection;

import factoring.FactorizationOfLongs;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialFact implements FactorizationOfLongs {


	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		//        for (int i = 2; i <= Math.sqrt(n); i++) {
		final long initialN = n;
		for (int factor = 3; factor*factor <= initialN; factor = factor + 2) {
			while (n%factor == 0) {
				if (primeFactors == null)
					return factor;
				primeFactors.add((long)factor);
				n /= factor;
			}
		}
		return n;
	}
}
