package factoring.fermat;

import factoring.math.PrimeMath;
import factoring.math.SquaresMod;

import java.util.Collection;

/**
 * Created by Thilo Harich on 09.01.2018.
 */
public class FermatAddFact extends FermatFact {


    @Override
    public long findFactors(long n, Collection<Long> factors) {
        for (int factor : smallFactors)
            if (n%factor == 0)
            {
                if (n/factor > 1)
                    factors.add(n/factor);
                return factor;
            }

        int sqrtN = (int) PrimeMath.sqrt(n);
        long xEnd = (n / minFactor + minFactor) / 2;
        int x = sqrtN;
//        long x2 = sqrtN;
        int right = (int) (x*x - n);
        for (; x <= xEnd; x++) {
//            long right2 = (x2*x2 - n);
//            if (right != right2)
//                System.out.println();
            if (SquaresMod.isSquare(right)) {
                long y = PrimeMath.sqrt(right);
                long factorHigh = x + y;
                factors.add(factorHigh);
                return x - y;
            }
            right += 2 * x + 1;
        }
        return n;
    }

}
