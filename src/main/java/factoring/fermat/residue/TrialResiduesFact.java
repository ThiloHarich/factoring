package factoring.fermat.residue;

import java.util.Collection;

import factoring.FactorizationOfLongs;
import factoring.trial.variant.TrialFactMod;

/**
 * Created by Thilo Harich on 28.06.2017.
 */
public class TrialResiduesFact implements FactorizationOfLongs {

	TrialFactMod smallFactoriser = new TrialFactMod();
	int factor = 1;

	@Override
	public long findFactors(long n, Collection<Long> factors) {
		// TODO estimate the speedup should be something like
		// exp ( ln(n')/ ln(ln(n')))
		final int estSpeedup = 10;
		int limit =  (int) Math.ceil(Math.pow(n, .5))/estSpeedup;
		limit = Math.max(19, limit);
		smallFactoriser.setLimit(limit);
		n = smallFactoriser.findFactors(n, factors);

		if (n<= limit) {
			n = smallFactoriser.findFactors(n, factors);
		}
		else
		{
			//            if (n < 1600)
			final FermatResiduesRec fermatRes = new FermatResiduesRec(limit);
			//            long factor = fermatRes.findPrimeFactors(n, factors);
			//            FermatFact fermatRes = new FermatFact(limit);
			return fermatRes.findFactors(n, factors, n);
		}
		return n;
	}
}
