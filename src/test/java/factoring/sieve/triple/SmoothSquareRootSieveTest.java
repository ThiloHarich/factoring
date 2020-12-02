package factoring.sieve.triple;

import org.junit.Test;

import java.math.BigInteger;

public class SmoothSquareRootSieveTest {

    @Test
    public void test87463(){
        SmoothSquareRootSieve sieve = new SmoothSquareRootSieve();
        sieve.factor(BigInteger.valueOf(87463));
    }
    @Test
    public void test8462977(){
        SmoothSquareRootSieve sieve = new SmoothSquareRootSieve();
        sieve.factor(BigInteger.valueOf(8462977));
    }

    @Test
    public void test769311827(){
        SmoothSquareRootSieve sieve = new SmoothSquareRootSieve();
        sieve.factor(BigInteger.valueOf(33637 * 22871));
    }
}
