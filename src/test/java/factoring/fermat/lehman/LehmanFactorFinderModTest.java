package factoring.fermat.lehman;

import factoring.fermat.lehman.playground.LehmanFactorFinderMod;
import factoring.math.PrimeMath;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertTrue;

/**
 * Created by Thilo Harich on 14.04.2018.
 */
public class LehmanFactorFinderModTest {

    @Test
    public void testNextX(){
        Random rand = new Random();
        LehmanFactorFinderMod test = new LehmanFactorFinderMod(1, 1);
//        test.getNextX();
        for (int j = 1; j < 16; j+=2) {
            int x = PrimeMath.mod(rand.nextInt(),  16);
            x += test.nextX[j][x%16];
            for (int i = 0; i < 16; i++) {
                int y = PrimeMath.mod(x * x - j, 16);
                assertTrue("j = " + j + " i = " + i,Math.sqrt(y * y) == y);
                x += test.nextX[j][x%16];
            }
        }
    }
}
