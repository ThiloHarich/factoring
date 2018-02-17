package factoring.math;

import java.util.Collection;

/**
 * Looks for solutions of x^2 + n = y^2 by looking for solutions modulo residueClassses[0].
 * For each solution x it calls {@link #findFactors(long, Collection, int, int)} with the next residueIndex.
 * Here the solutions x+residueClassses[residueIndex], x+2*residueClassses[residueIndex], ... will be merged with the solutions
 * Modulo residueClassses[residueIndex+1].
 * This is done by checking if the value x + i* residueClassses[residueIndex] has a solution modulo residueClassses[residueIndex+1].
 * This is done by checking against a {@link #initX(long)}}.
 *
 * @author thiloharich
 *
 */
public class FermatResiduesRecursive  extends FermatResiduesBooleanArray{

	int sqrtN = -1;
//	public static int[] residueClasses;
	long previousResidueProduct;

	/**
	 * The list of primes to be used as residues
	 */
//	static FermatResiduesRec[] factorizers = new FermatResiduesRec[residuesArray.length+1];

	public FermatResiduesRecursive nextFermatResidue;
	/**
	 * to determine the end of the for loop for x we need the biggest possible factor of the
	 * number nOrig to be factorized
	 */
	private static long minFactor = 3;

	public FermatResiduesRecursive(int mod) {
		super(mod);
	}
		/**
         * Compute the squares xArray^2 modulo mod2Pow. We avoid expensive % calculations by
         * using (xArray+1)^2 - xArray^2 = 2x + 1. Since (-xArray)^2 = xArray^2 we only have to iterate until mod2Pow/2 +1
         */
	public FermatResiduesRecursive(int mod, long previousResidueProduct, long n) {
		this.mod = mod;
		this.previousResidueProduct = previousResidueProduct;
		initSquares(mod);
		initX(n);
	}

	public FermatResiduesRecursive() {
	}


	/**
	 * initializes the initX of all factorizers (with different modulus).
	 * @param n
	 * @return If n is dividable by one of the mods this mod2Pow is returned, otherwise 0;
	 */
	private int initX4All(int n) {
		if (mod > 0) {
			int nMod = PrimeMath.mod(n, mod);
			// TODO if mod2Pow = p^i also check p^j , j<i
			// If we do a check for small factors first, this is not needed
			if (nMod == 0)
				return mod;
			initX(nMod);
			nextFermatResidue.initX(n);
		}
		return 0;
	}

	public static FermatResiduesRecursive create (long n, int ... residueClasses)
	{
		FermatResiduesRecursive fermat = new FermatResiduesRecursive(residueClasses[0], 0, n);
		FermatResiduesRecursive.residueClasses = residueClasses;
		fermat.initSquares(fermat.mod);
		FermatResiduesRecursive root = fermat;
		long modProd = 1;
		int i=0;
		for (; i < residueClasses.length-1; i++) {
			int mod = residueClasses[i];
			fermat.nextFermatResidue = new FermatResiduesRecursive(mod, modProd, n);
			if (mod > 0)
				modProd *= mod;
			fermat = fermat.nextFermatResidue;
		}
		fermat.nextFermatResidue = new FermatResiduesRecursive(0, modProd, n);
		return root;
	}


	public static FermatResiduesRecursive create(long n) {
		int [] residueClasses = calculateResidueClasses(n/6 + 2);
		return FermatResiduesRecursive.create(n, residueClasses);
	}



	/**
	 * generates all possible x values for the given mod2Pow and calls {@link #findFactors(long, Collection, int, int)} for
	 * x. if a real factor is returned, it will be returned, otherwise the next x is going to be generated.
	 * Precondtition: squares are beeing calculated
	 */
	public long findFactors(long n, Collection<Long> factors, int x, int mod) {
		if (this.mod == 0) {
			int modOld = (int) previousResidueProduct;
			long xEnd = (n / minFactor + minFactor) / 2;
			// adjust sqrtN to be a multiple of modOld
			int xBegin = (sqrtN / modOld) * modOld;
			long xL = xBegin + x;
//			if (modulus * mod2Pow  >  xRange/SEARCH_COUNT) {
			while (xL <= xEnd) {
				long right = xL * xL - n;
				if (SquaresMod.isSquare(right)) {
					long y = PrimeMath.sqrt(right);
					long factorHigh = xL + y;
					long factorLow = xL - y;
					factors.add(factorHigh);
					return factorLow;
				}
				xL+= modOld;
			}
			return n;
		}
		// TODO put the x with x^2 - n = 0 mod2Pow mod2Pow to the front
		long factor = findFactorsByMerge(n, factors, x, mod);
			if (factor != n)
				return factor;
		return n;
	}


	/**
	 * Here the solutions x modulo modOld , ... will be merged with the solutions
	 * Modulo modNew.
	 * This is done by checking if the value x + i* modOld has a solution modulo modNew.
	 * we calculate (x + i* modOld) mod2Pow modNew and check if {@link #initX(long)}}.
	 *
	 * complexity: #oldSolutions * modNew
	 *
	 * @param n
	 * @param factors
	 * @param x
	 * @param modOld
	 * @return
	 */
	public long findFactorsByMerge(long n, Collection<Long> factors, int x, int modOld) {
		// step is the lowest number greater modulus which is a multiple of modNew
		// TODO avoid division here
		final int step = (modOld / mod + 1) * mod;

		// avoid % here, maybe we can keep old xMod
		int xMod = x;
		while (xMod >= mod)
			xMod -= mod;
		while (x <  modOld * mod)
		{
			if (xMask[xMod]) {
				long factor = nextFermatResidue.findFactors(n, factors, x, mod);
				if (factor != n)
					return factor;
			}
			x += modOld;
			xMod += modOld;
			// We avoid calculating expensive % operation
			xMod -= xMod >= step ? step : (step - mod);
		}
		return n;
	}
}
