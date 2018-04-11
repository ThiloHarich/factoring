package factoring.math;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import factoring.trial.TrialFact;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by Thilo Harich on 16.12.2017.
 */
public class FermatResiudesTest {

    /**
     * x^2 - 23 = y^2  mod 24
     * x^2 + 1 = y^2  mod 24
     *
     * 0 -> 0
     * 1 -> 1
     * 2 -> 4
     * 3 -> 9
     * 4 -> 16
     * 5 -> 1
     * 6 -> 12
     * 7 -> 1
     * 8 -> 16
     * 9 -> 9
     * 10-> 4
     * 11-> 1
     * 12-> 0
     * x = 0, 12
     *
     * n^2 mod 12
     * 1 -> 1
     * 5 -> 1
     * 7 -> 49 -> 1
     * 11 -> 121 -> 1
     *
     * Multiply n with n mod residueclass
     */
    @Test
    public void testNumber () {
        // 16 * 27
//        int number = 2*2*2*2*2 * 3*3 * 5 * 7 * 11;
//        int number = 2*2*2 * 3*3 * 5 * 7 * 11;
//        int number = 8*9*5*7*11;
        int number = 2;
        for (int i = number; i <= number + 32; i++)
        {
            FermatResidueMergeByInversion fermat = new FermatResidueMergeByInversion(i);
            Stream<Boolean> stream = IntStream.range(0, fermat.squares.length)
                    .mapToObj(idx -> fermat.squares[idx]);
            long countSquares = stream.filter(b -> b == true).count();

            TrialFact factorizer = new TrialFact();

            long max = 0;
            long min = Long.MAX_VALUE;
            int minIndex = 0;
            double avg = 0;
            int num = 0;
            String details = "";
            for (int n = 1; n< i; n++) {
                if (PrimeMath.gcd(i, n) == 1) {
                    fermat.initX(n);
                    long countFermat = Arrays.stream(fermat.xArray).boxed().collect(Collectors.toList()).indexOf(-1);
                    max = Math.max(max, countFermat);
                    if (countFermat < min)
                        minIndex = n;
                    min = Math.min(min, countFermat);
                    avg += countFermat;
                    num++;
                    details += countFermat + " , ";
                }
            }
//            double quality = Math.pow(Math.log(i),1.2) / Math.log(min);
//            double quality = Math.log(i) / Math.log(max);
            double quality = (i + 0.0) / (min);
//            if (quality > 5)
            {
                float ifl = (float)i;
                String factors = String.format("%-20s", factorizer.printFactorization(i));
                System.out.print(i + "(" + factors + ") \t: " + countSquares + " ; ");
                System.out.print(" exp : " + (ifl+1)/(2*ifl) * ifl);
                System.out.print(" max : " + max + " min : " + min + "(" + minIndex + ") avg : " + (avg/ num)  + " quality : " + quality + "              : " + details);
                System.out.println();
//                TreeMultiset<Long> factorSet = factorizer.findAllPrimeFactors(i);
//                for (int j=1; j<i; j++)
//                {
//                    boolean hasDivisor = hasDivisor(factorSet, j);
//                    if (!hasDivisor) {
//                        //
//                        Integer multiplier = null;
//                        for (int k=0; k<i && multiplier == null; k++)
//                            if (PrimeMath.mod(k*j, i) == minIndex) {
//                                multiplier = k;
//                            }
//                        if (multiplier == null)
//                            System.out.println("no Multiplier found for " + j);
//
//                    }
//                }
//                int n = minIndex;
//                int[] xArray = fermat.xArray(n);
//                System.out.println(Arrays.stream(xArray).boxed().filter(x -> x> 0).map(Object::toString).collect(Collectors.joining(",")));
            }
        }
    }

    public boolean hasDivisor(TreeMultiset<Long> factorSet, int j) {
        Iterator<Long> iter = factorSet.elementSet().iterator();
        while(iter.hasNext()) {
            long factor = iter.next();
            if (j % factor == 0) {
                return true;
            }
        }
        return false;
    }
}
