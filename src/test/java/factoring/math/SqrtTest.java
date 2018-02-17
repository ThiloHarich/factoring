package factoring.math;

import org.junit.Test;

/**
 * Created by Thilo Harich on 08.01.2018.
 */
public class SqrtTest {

    @Test
    public void testPerformance ()
    {
        int bits = 16;
        long begin = (1L << bits) + 3534;
        long range = 10000000L;
        int n3 = 1 << (bits / 3);

        long start = System.currentTimeMillis();
        hart(begin, range, n3);
        long end  = System.currentTimeMillis();
        System.out.println(" time hart  : " + (end- start));

        start = System.currentTimeMillis();
        floor(begin, range, n3);
        end  = System.currentTimeMillis();
        System.out.println(" time floor : " + (end- start));

        start = System.currentTimeMillis();
        hart(begin, range, n3);
        end  = System.currentTimeMillis();
        System.out.println(" time hart  : " + (end- start));

        start = System.currentTimeMillis();
        floor(begin, range, n3);
        end  = System.currentTimeMillis();
        System.out.println(" time floor : " + (end- start));

    }

    public void hart(long begin, long range, int n3) {
        boolean found = false;
//        long  n = begin;
        for(long n = begin; n < begin + range; n++)
        {
            long sqrt = PrimeMath.sqrt(n);
            long x = sqrt;
            for (; x< sqrt + n3; x++) {
                long x2 = x*x;
                long right = x2 - n;
                if (PrimeMath.isSquare(right)) {
                    found = ! found;
                }
            }
        }
    }
    public void floor(long begin, long range, int n3) {
        boolean found = false;
        int onLevel = 3;
//        long  n = begin;
        for(long n = begin; n < begin + range; n++)
        {
            int sqrt = (int) PrimeMath.sqrt(n);
            int x = sqrt;
            long x2 = x*x;
            int right = (int) (x2 - n);
            int x21 = 2 * x + 1;
            for (; x21< (sqrt + n3)*2; x21+=2) {
                if (PrimeMath.isSquare(right)) {
                    found = ! found;
                }
                right += x21;
            }
        }
    }

}
