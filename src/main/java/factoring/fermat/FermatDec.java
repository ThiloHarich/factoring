package factoring.fermat;

import factoring.math.PrimeMath;

import java.util.Collection;

/**
 * n = (xArray+y)^2 - (xArray-y)^2 = xArray^2 - y^2 , 0 <= y <= sqrt(n)
 *  n + y^2 = xArray^2
 * Created by Thilo Harich on 02.03.2017.
 */
public class FermatDec extends FermatFact {
    @Override
    public long findFactors(long n, Collection<Long> factors) {
        long sqrtN = (long) Math.floor(Math.sqrt(n));
        for (long y = sqrtN; y >= 0; y--) {
            long right = n + y*y;
//            long sqrtR = PrimeMath.sqrt(right);
//            long error = n + sqrtR*sqrtR;
//            right = n + (y-error)*(y-error);

            if(PrimeMath.isSquare(right))
            {
                long x = (long) Math.sqrt(right);
                factors.add((x+y));
                return n/(x+y);
            }
        }
        return n;
    }
}
