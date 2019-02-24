package factoring.math;

import org.junit.Test;

public class QuadraticDiophantineModTest {

	@Test
	public void printMultipliers() {
		for (int m = 2; m <= 27; m++) {
			System.out.println("mod : " + m);
			final QuadraticDiophantineModBit mod = new QuadraticDiophantineModBit(m);
			for (int i=0; i< m; i++) {
				final int[] x = mod.xArray(i);
				int count = 0;
				for (final int element : x) {
					if (element < 0)
						break;
					count++;
				}
				System.out.println( i + " : " + count);
			}
		}

	}
	@Test
	public void printHart() {

		final int mod = 16;
		final QuadraticDiophantineModBit solutionsMod = new QuadraticDiophantineModBit(mod);

		solutionsMod.lehmanSolutions();
	}
	@Test
	public void print63() {

		final int mod = 5*7*64;
		final QuadraticDiophantineModBit solutionsMod = new QuadraticDiophantineModBit(mod);

		solutionsMod.lehmanSolutions();
		final int[] x = solutionsMod.xArray(0);

		for (final int element : x) {
			if (element < 0)
				break;
			System.out.print(", " + element);
		}
	}
}
