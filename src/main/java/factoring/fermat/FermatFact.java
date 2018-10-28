package factoring.fermat;

import java.util.Collection;

import factoring.FactorizationOfLongs;
import factoring.math.PrimeMath;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class FermatFact implements FactorizationOfLongs {

	protected int [] smallFactors = {2,3,5,7,11,13,17,19};
	protected int minFactor;

	public FermatFact(int minFactor)
	{
		this.minFactor = minFactor;
	}

	// TODO take smallFactors[smallFactors.lenght - 1]
	public FermatFact()
	{
		minFactor = 3;
	}
	//    /**
	//     * In contrast to the trial division we can not be sure that a found factor is a prime
	//     * so we have to factorize it again
	//     * @param n
	//     * @return
	//     */
	//    public TreeMultiset<Long> findAllPrimeFactors(long n) {
	//        TreeMultiset<Long> primeFactors = TreeMultiset.create();
	//        TreeMultiset<Long> factors = findAllFactors(n);
	//        while (factors.size() > 0) {
	//            Multiset.Entry<Long> entry = factors.pollLastEntry();
	//            long factor = entry.getElement();
	////            if (factor*factor > n || factor < 30)
	////                primeFactors.add(factor, entry.getCount());
	////            else
	//                {
	//                TreeMultiset<Long> newFactors = TreeMultiset.create();
	//                storeFactors(factor, newFactors);
	//                // if no factor found it is a prime factor
	//                if (newFactors.lastEntry().getElement() == factor) {
	//                    primeFactors.add(factor, entry.getCount());
	//                    newFactors.pollLastEntry();
	//                }
	//                else {
	//                    for (long newFactor : newFactors) {
	//                        factors.add(newFactor, entry.getCount());
	//                    }
	//                }
	//            }
	//        }
	//        return primeFactors;
	//    }
	//
	//    /**
	//     *
	//     * @param n
	//     * @param factors
	//     * @return
	//     */
	//    public Collection<Long> storeFactors(long n, Collection<Long> factors) {
	//        // first make the number odd, this is required by the fermat factorizations
	//        while ((n & 1) == 0)
	//        {
	//            factors.add(2l);
	//            n = n >> 1;
	//        }
	//        if (n>1) {
	//            long remainder = findFactors(n, factors);
	//            if (remainder != 1)
	//                factors.add(remainder);
	//        }
	//        return factors;
	//    }

	@Override
	public long findFactors(long n, Collection<Long> factors) {
		// This part is needed for the mod variants. It also improves performance
		for (final int factor : smallFactors)
			if (n%factor == 0)
			{
				if (n/factor > 1)
					factors.add(n/factor);
				return factor;
			}

		final long sqrtN = (long) Math.ceil(Math.sqrt(n));
		final long xEnd = (n / minFactor + minFactor) / 2;
		for (long x = sqrtN; x <= xEnd; x++) {
			final long right = x*x - n;
			if (PrimeMath.isSquare(right)) {
				final long y = PrimeMath.sqrt(right);
				final long factorHigh = x + y;
				factors.add(factorHigh);
				return x - y;
			}
		}
		return n;
	}
}
