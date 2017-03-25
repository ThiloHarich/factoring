package factoring;

import java.util.OptionalInt;
import java.util.stream.IntStream;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class TrialStream implements Factorizer {

    @Override
    public int findFactor(long n) {
        return IntStream.range(3, (int) Math.sqrt(n) + 1).filter(p -> n%p == 0).findFirst().orElse(-1);
//        return IntStream.iterate(3, i -> i + 2).limit((int) Math.sqrt(n)/2 + 1).filter(p -> n%p == 0).findFirst().orElse(-1);
    }
}
