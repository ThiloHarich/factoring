package factoring;

/**
 * Created by Thilo Harich on 02.03.2017.
 */
public class FermatDec extends FermatInc {
    @Override
    public int findFactor(long n) {
        long sqrtN = (long) Math.floor(Math.sqrt(n));
        for (long y = sqrtN; y >= 0; y--) {
            long right = n + y*y;
//            long sqrtR = MyMath.sqrt(right);
//            long error = n + sqrtR*sqrtR;
//            right = n + (y-error)*(y-error);

            if(MyMath.isSquare(right))
            {
                long x = (long) Math.sqrt(right);
                return (int) (x+y);
            }

        }
        return -1;
    }
}
