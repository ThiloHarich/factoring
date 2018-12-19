package factoring;

import de.tilman_neumann.jml.factor.FactorAlgorithmBase;
import factoring.fermat.lehman.Lehman_TillSimple;
import factoring.fermat.lehman.Lehman_TillSimple2;

//import de.tilman_neumann.math.factor.CombinedFactorAlgorithm;
//import de.tilman_neumann.math.factor.FactorAlgorithm;

public class PerformanceHard2 extends PerformanceHard{

	public static void main(String[] args) {
		final FactorAlgorithmBase factorizer2 = new Lehman_TillSimple2(1);
		//		final FactorizationOfLongs factorizer2 = new LehmanFactorFinderMod12(bits, 2.f, false);
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2.f, false);
		final FactorAlgorithmBase factorizer3 = new Lehman_TillSimple(1);
		//		final FactorizationOfLongs factorizer1 = new LehmanFactorFinder(bits, 2f, false);
		semiprimes = makeSemiPrimesList(bits, smallFactorBits, numPrimes);
		test2(factorizer3);
		//		test2(factorizer1);
		//		findFactors(factorizer1, semiprimes, 1);

		test2(factorizer2);
	}


}
