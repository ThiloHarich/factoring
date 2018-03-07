package factoring.math;

import de.tilman_neumann.math.base.bigint.Gcd63;
import org.junit.Test;

import java.util.Random;

/**
 * Created by Thilo Harich on 07.03.2018.
 */
public class GdcTest {

    @Test
    public void testPerformance ()
    {
        Random rand = new Random();
        long begin = rand.nextLong();
        long range = 5000000;
        Gcd63 gcd63 = new Gcd63();
        long offset = rand.nextLong();

        long start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
            PrimeMath.gcd10(begin, other);
        long end  = System.currentTimeMillis();
        System.out.println(" time thilo  : " + (end- start));

        start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
           gcd63.gcd(begin + offset, other);
        end  = System.currentTimeMillis();
        System.out.println(" time till : " + (end- start));

        start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
            PrimeMath.gcd10(begin, other);
        end  = System.currentTimeMillis();
        System.out.println(" time thilo  : " + (end- start));

        start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
            gcd63.gcd(begin, other);
        end  = System.currentTimeMillis();
        System.out.println(" time till : " + (end- start));

        start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
            PrimeMath.gcd10(begin, other);
        end  = System.currentTimeMillis();
        System.out.println(" time thilo  : " + (end- start));

        start = System.currentTimeMillis();
        for (long other = begin + offset; other < begin + offset + range; other++)
            gcd63.gcd(begin, other);
        end  = System.currentTimeMillis();
        System.out.println(" time till : " + (end- start));

    }

}
