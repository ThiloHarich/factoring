package factoring;

import factoring.fermat.lehman.LehmanFactorFinder;
import factoring.fermat.lehman.LehmanFactorFinderMod12;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard2 extends PerformanceHard{

	public static void main(String[] args) {
		final FactorizationOfLongs factorizer2 = new LehmanFactorFinderMod12(bits, 2.f, false);
		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2.f, false);
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorFinderRange(bits, 2f, false);
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2f, false);
		semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);
		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, 1);

		test2(factorizer2);
	}

}
