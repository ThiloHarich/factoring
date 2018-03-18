package factoring.trial;

import factoring.FindPrimeFact;

import java.util.Collection;
import java.util.stream.IntStream;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialStream  extends FindPrimeFact {

    @Override
    public long findPrimeFactors(long n, Collection<Long> primeFactors) {
        primeFactors.add((long)IntStream.range(3, (int) Math.sqrt(n) + 1).filter(p -> n%p == 0).findFirst().orElse(-1));
        return n;
//        return IntStream.iterate(3, i -> i + 2).limit((int) Math.sqrt(n)/2 + 1).filter(p -> n%p == 0).findFirst().orElse(-1);
    }
}
