package factoring;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Looks for solutions x for  x^2 - n = y^2 mod m
 *
 * for f(x) = x^2 - n mod m we
 * define a Map_n x^2 - n  -> x
 *
 * @author thiloharich
 *
 */
public class SquaresMod {

	public static final int UNDEFINED = -1;

	public static class Pair {
		int first = UNDEFINED;
		int second = UNDEFINED;
		public Pair(int x) {
			first = x;
		}
	}

	/**
	 * TODO lamda for f?
	 * @param n
	 * @param mod
	 * @return
	 */
	public static Map<Integer, Pair> solutionQuadraticFunction (int n, int mod){
		final HashMap<Integer, Pair> squares = new HashMap<Integer, Pair>();
		for (int x=0; x<mod; x++)
		{
			final int f = mod(x*x - n, mod);
			Pair pair;
			if (!squares.containsKey(f))
			{
				pair = new Pair(x);
				squares.put(f, pair);
			}
			else
			{
				pair = squares.get(f);
				pair.second = x;
			}
		}
		return squares;
	}

	public static List<Integer> solutionsFermatMod (int n, int mod)
	{
		final Map<Integer, Pair> solutions0 = solutionQuadraticFunction(0, mod);
		final Map<Integer, Pair> solutionsN = solutionQuadraticFunction(n, mod);

		final Set<Integer> roots0 = solutions0.keySet();
		final Set<Integer> rootsN = solutionsN.keySet();

		// look for values x^2 - n that are also y^2
		rootsN.retainAll(roots0);

		// TODO solutions with y^2 = 0 are good solutions

		final List<Integer> solutions = new ArrayList<Integer>();
		for (final Pair pair : solutionsN.values()) {
			solutions.add(pair.first);
			if (pair.second != UNDEFINED)
				solutions.add(pair.second);
		}
		return solutions;
	}

	public static Collection<Integer> combineSolutions (Collection<Integer> sol1, int mod1,Collection<Integer> sol2, int mod2)
	{
		final int mod = mod1 * mod2;
		final List<Integer> solutions = new ArrayList<Integer>();
		for (final Integer i1 : sol1) {
			for (final Integer i2 : sol2) {
				final int sol = (i1 * mod2 + i2 * mod1) % mod;
				solutions.add(sol);
			}
		}

		return solutions;
	}

	public static int mod (int x, int mod)
	{
		int xMod = x % mod;
		if (xMod < 0)
			xMod += mod;
		return xMod;
	}

}
