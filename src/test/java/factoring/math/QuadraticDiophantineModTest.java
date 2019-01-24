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
}
