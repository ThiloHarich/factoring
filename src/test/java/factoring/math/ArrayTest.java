package factoring.math;

import org.junit.Test;

import static junit.framework.TestCase.assertEquals;

/**
 * Created by Thilo Harich on 15.04.2018.
 */
public class ArrayTest {


    private static int [] mod16 = new int[16];

    static {
        for (int i = 0; i < 16; i++) {
            mod16[i] = 4;
        }
    }

    @Test
    public void testModArray(){
        int loop = 1000000000;
        int sum1 = 0;
        int n = 6378463;
        long start = System.currentTimeMillis();
        for (int i = 0; i < loop; i++) {
            int x = 4;
            sum1 += x*x - n;
        }
        long end = System.currentTimeMillis();
        System.out.println("Loop  : " + (end - start));
        int sum2 = 0;
        start = System.currentTimeMillis();
        for (int i = 0; i < loop; i++) {
            int x = i%2 == 0 ? 4 : 2;
            sum2 += x*x - n;
        }
        end = System.currentTimeMillis();
        System.out.println("if    : " + (end - start));
        int sum3 = 0;
        start = System.currentTimeMillis();
        for (int i = 0; i < loop; i+=2) {
            int x = 4;
            sum3 += x*x - n;
        }
        for (int i = 1; i < loop; i+=2) {
            int x = 2;
            sum3 += x*x - n;
        }
        end = System.currentTimeMillis();
        System.out.println("2Loops : " + (end - start));
        int sum4 = 0;
        start = System.currentTimeMillis();
        for (int i = 0; i < loop; i+=2) {
            int x = 4;
            sum4 += x*x - n;
        }
        for (int i = 1; i < loop; i+=2) {
            int x = 2;
            sum4 += x*x - n;
        }
        end = System.currentTimeMillis();
        System.out.println("2Loops : " + (end - start));
        assertEquals(sum3, sum2);
    }
}
