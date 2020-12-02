package factoring.sieve.triple;

import org.junit.Test;

import java.math.BigInteger;

public class SmoothNumbersTest {

    SmoothNumbersSieve sieve = new SmoothNumbersSieve();

    @Test
    public void test87463(){
        sieve.factor(BigInteger.valueOf(87463));
    }
    @Test
    public void test13_17(){
        sieve.factor(BigInteger.valueOf(13*17));
    }
    @Test
    public void test13_23(){
        sieve.factor(BigInteger.valueOf(29*23));
    }
    @Test
    public void test8462977(){
        sieve.factor(BigInteger.valueOf(8462977));
    }

    @Test
    public void test769311827(){
        sieve.factor(BigInteger.valueOf(33637 * 22871));
    }
    @Test
    public void test8000009_4000037(){
        final long n = 8_000_0009l * 4_000_037l;
        sieve.factor(BigInteger.valueOf(n));
    }

    @Test
    public void test28Bit(){
 //        final long n = 150000001l * 300000007l;
        final long n = 150000001l * 280000027l;
        sieve.factor(BigInteger.valueOf(n));
    }

}
