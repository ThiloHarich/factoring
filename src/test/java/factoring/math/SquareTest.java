package factoring.math;

import com.google.common.primitives.Ints;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by Thilo Harich on 12.08.2017.
 */
public class SquareTest {

    public static void main(String[] args)
        {
            SquareTest s = new SquareTest();
            s.perf();
        }

    @Test
    public void xArrayBigTest()
    {
        int mod = 256;
        QuadraticDiophantineModBig squaresBig = new QuadraticDiophantineModBig(mod);
        SquaresModBitSet squaresBit = new SquaresModBitSet(mod);
//        QuadraticDiophantineModBit squaresBit2 = new QuadraticDiophantineModBit(mod2);

        for (int k = 0; k < mod; k++) {
            System.out.print(k + " : ");
            int nMod = PrimeMath.mod(k, mod);
//            int nMod64 = PrimeMath.mod2Pow(-k, 64);
            int[] xCandidatesModBit = squaresBit.xArray(nMod);
            List<Integer> list = Ints.asList(xCandidatesModBit);
            int[] actual = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();

            int[] xCandidatesModBig = squaresBig.xArray(nMod);
            List<Integer> list3 = Ints.asList(xCandidatesModBig);
            int[] listBig = list3.stream().limit(list3.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
            Arrays.stream(listBig).forEach(e -> System.out.print(e + ","));
            System.out.print("\n" + k + " : ");
            Arrays.stream(actual).forEach(e -> System.out.print(e + ","));
            System.out.println();
            assertTrue(Arrays.equals(listBig, actual));
//            128 -> 10
//            256 -> 32
//            4*3*5= 60 -> 8
//            4*3*5*7= 420 -> 16
        }
    }
//    @Test
//    public void mergeMultipleTest()
//    {
//        int mod1 = 81;
//        int mod2 = 49;
////        int mod3 = 7;
//        int mod2Pow = mod1*mod2;
//        FermatResiduesRec squaresBit = new FermatResiduesRec(mod2Pow);
//        FermatResiduesBooleanArray squaresBit1 = new FermatResiduesBooleanArray(mod1, mod2);
//
//        for (int k = 0; k < mod2Pow; k++) {
//            if (k%mod1 != 0 && k%mod2 !=0) {
//                int nMod = PrimeMath.mod2Pow(-k, mod2Pow);
//                squaresBit.initX(nMod);
//                int[] xCandidatesModBit = squaresBit.xArray;
//                List<Integer> list = Ints.asList(xCandidatesModBit);
//                int[] actual = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
//
//                int[] xCandidatesModBit1 = squaresBit1.x4ResidueClasses(-k);
//                List<Integer> list3 = Ints.asList(xCandidatesModBit1);
//                int[] mergeList = list3.stream().limit(list3.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
//                Arrays.sort(mergeList);
//                Arrays.sort(actual);
//                Arrays.stream(mergeList).forEach(e -> System.out.print(e + ","));
//                System.out.println();
//                Arrays.stream(actual).forEach(e -> System.out.print(e + ","));
//                System.out.println();
//                assertTrue(Arrays.equals(mergeList, actual));
//            }
//        }
//    }
    @Test
    public void mergeTest()
    {
        int mod1 = 3;
        int mod2 = 5;
        int mod3 = 4;
        int mod = 3*5*4;
        QuadraticDiophantineModBit squaresBit = new QuadraticDiophantineModBit(mod);
        QuadraticDiophantineModBit squaresBit1 = new QuadraticDiophantineModBit(mod1);
        QuadraticDiophantineModBit squaresBit2 = new QuadraticDiophantineModBit(mod2);
        QuadraticDiophantineModBit squaresBit3 = new QuadraticDiophantineModBit(mod3);

        for (int k = 0; k < mod; k++) {
            int nMod = SquaresModArray.mod(-k, mod);
            int[] xCandidatesModBit = squaresBit.xArray(nMod);
            List<Integer> list = Ints.asList(xCandidatesModBit);
            int[] actual = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();

            int nMod1 = SquaresModArray.mod(-k, mod1);
            int[] xCandidatesModBit1 = squaresBit1.xArray(nMod1);
            int nMod2 =  SquaresModArray.mod(-k, mod2);
            int[] xCandidatesModBit5 = squaresBit2.xArray(nMod2);
            int[] xCandidatesModBit2 = squaresBit2.merge(xCandidatesModBit1, squaresBit1.xLength(), mod1);
            int nMod3 = SquaresModArray.mod(-k, mod3);
            int[] xCandidatesModBit4 =squaresBit3.xArray(nMod3);
            int[] xCandidatesModBit3 = squaresBit3.merge(xCandidatesModBit2, squaresBit2.xLength(), mod1*mod2);
            List<Integer> list3 = Ints.asList(xCandidatesModBit3);
            int[] mergList = list3.stream().limit(list3.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
            Arrays.sort(mergList);
            Arrays.stream(mergList).forEach(e -> System.out.print(e + ","));
            System.out.println();
            Arrays.stream(actual).forEach(e -> System.out.print(e + ","));
            System.out.println();
            assertTrue(Arrays.equals(mergList, actual));
        }
    }
    @Test
    public void mergeIteratorTest()
    {
        int mod1 = 3;
        int mod2 = 5;
        int mod = mod1*mod2;
        QuadraticDiophantineModBit squaresBit = new QuadraticDiophantineModBit(mod);
        FermatResiduesRecursive squaresBit1 = new FermatResiduesRecursive(mod1);
        FermatResiduesBooleanArray squaresBit2 = new FermatResiduesBooleanArray(mod2);

        for (int k = 0; k < mod; k++) {
            int nMod = PrimeMath.mod(-k, mod);
            int[] xCandidatesModBit = squaresBit.xArray(nMod);
            List<Integer> list = Ints.asList(xCandidatesModBit);
            List<Integer> expected = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).collect(Collectors.toList());

            int nMod1 = PrimeMath.mod(-k, mod1);
            squaresBit1.initX(nMod1);
            int[] xCandidatesModBit1 = squaresBit1.xArray;
            int nMod2 =  PrimeMath.mod(-k, mod2);
            squaresBit2.initX(nMod2);
            FermatResiduesBooleanArray.MergeIterator xCandidatesModBit2 = squaresBit2.mergeIterator(xCandidatesModBit1, squaresBit1.xLength(), mod1);
            List<Integer> mergeList = new ArrayList<>();
            while (xCandidatesModBit2.hasNext()) {
                int element = xCandidatesModBit2.next();
                mergeList.add(element);
            }
            Collections.sort(mergeList);
            mergeList.forEach(e -> System.out.print(e + ","));
            System.out.println();
            Collections.sort(expected);
            expected.forEach(e -> System.out.print(e + ","));
            System.out.println();
            assertEquals(expected, mergeList);
        }
    }
    @Test
    public void mergePerf()
    {
        int mod1 = 3;
        int mod2 = 5;
        int mod3 = 4;
        int mod = 3*4*5;
        QuadraticDiophantineModBit squaresBit1 = new QuadraticDiophantineModBit(mod1, mod2, mod3);
        int bound = 1000000;

        for (int i=0; i< 10; i++) {
            long t1 = System.currentTimeMillis();
            for (int k = 0; k < bound; k++) {
                QuadraticDiophantineModBit squaresBit = new QuadraticDiophantineModBit(mod);
                int nMod = PrimeMath.mod(-k, mod);
                int[] xCandidatesModBit = squaresBit.xArray(nMod);
            }
            long t2 = System.currentTimeMillis();

            System.out.println("Time all : " + (t2-t1+0.0)/1000);

            for (int k = 0; k < bound; k++) {
                int[] xCandidatesModBit1 = squaresBit1.x4ResidueClasses(-k);
            }
            t1 = System.currentTimeMillis();
            System.out.println("Time split: " + (t1-t2+0.0)/1000);
        }
    }
    @Test
    public void mergePerfPow2()
    {
        int mod = 256;
        QuadraticDiophantineModBit squaresBit1 = new QuadraticDiophantineModBit(mod);
        int bound = 100000;

        for (int i=0; i< 10; i++) {
            long t1 = System.currentTimeMillis();
            for (int k = 0; k < bound; k++) {
                QuadraticDiophantineModBit squaresBit = new QuadraticDiophantineModBit(mod);
                int nMod = PrimeMath.mod(-k, mod);
                int[] xCandidatesModBit = squaresBit.xArray(nMod);
            }
            long t2 = System.currentTimeMillis();

            System.out.println("Time bit    : " + (t2-t1+0.0)/1000);

            for (int k = 0; k < bound; k++) {
                SquaresModBitSet squaresBit = new SquaresModBitSet(mod);
                int nMod = PrimeMath.mod(-k, mod);
                int[] xCandidatesModBit = squaresBit.xArray(nMod);
            }
            t1 = System.currentTimeMillis();
            System.out.println("Time bitSet : " + (t1-t2+0.0)/1000);
        }
    }
    @Test
    public void correct() {
        int mod =81;
        for (int k = 0; k < mod; k++) {
            SquaresModBitSet squares = new SquaresModBitSet(mod);

            int nMod =  PrimeMath.mod(-k, mod);
            int[] xCandidatesMod = squares.xArray(nMod);
            List<Integer> list =  Ints.asList(xCandidatesMod);
            int[] expected = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
            Arrays.sort(expected);

            FermatResiduesRecursive squaresBit = new FermatResiduesRecursive(mod);

            squaresBit.initX(nMod);
            int[] xCandidatesModBit = squaresBit.xArray;
            list =  Ints.asList(xCandidatesModBit);
            int[] actual = list.stream().limit(list.indexOf(Integer.valueOf("-1"))).mapToInt(Integer::intValue).toArray();
            Arrays.sort(actual);

            assertTrue(Arrays.equals(expected, actual));
        }

    }
    @Test
    public void perf() {

        System.out.println("start");
        long time = System.currentTimeMillis();
        int mod = 8;
        boolean found = false;
        int limit = 2000000;
//        for (int i = 0; i < limit; i++) {
//            SquaresMod squaresMask = SquaresMod.mod2Pow(mod2Pow);
//
//            long nMod = i % mod2Pow;
//            SquaresMod squaresPlusN = squaresMask.plusRight(-nMod);
//            SquaresMod disjunction = squaresPlusN.disjuctionRight(squaresMask);
//            Collection<Integer> xCandidatesMod = disjunction.values.values();
//
//            if (xCandidatesMod.size() > 0)
//                found = true;
//        }
//        if (found)
//            System.out.println(System.currentTimeMillis() - time);
//        time = System.currentTimeMillis();
//        for (int k = 0; k < limit; k++) {
//            SquaresModArray squaresMask = SquaresModArray.mod2Pow(mod2Pow);
//
//            long nMod = SquaresModArray.mod2Pow(-k, mod2Pow);
//            SquaresModArray squaresPlusN = squaresMask.plusRight(nMod);
//            squaresPlusN.disjuctionRight(squaresMask);
//            int[][] values = squaresPlusN.values;
//            int[] xCandidatesMod = new int[2*mod2Pow];
//            int j = 0;
//            for (int i = 0; i < mod2Pow; i++)
//            {
//                if (values[i][0] >= 0)
//                {
//                    xCandidatesMod[j++] = values[i][0];
//                    // do not add the same value twice; can only happen with 0 and mod2Pow
//                    if (Math.abs(values[i][0] - values[i][1]) != mod2Pow)
//                        xCandidatesMod[j++] = values[i][1];
//                }
//            }
//            xCandidatesMod[j] = -1;
//        }
//        System.out.println(System.currentTimeMillis() - time);

        time = System.currentTimeMillis();
        for (int k = 0; k < limit; k++) {
            SquaresModArray1Dim squares = SquaresModArray1Dim.mod(mod);

            long nMod = SquaresModArray.mod(-k, mod);
            SquaresModArray1Dim squaresPlusN = squares.plusRight(nMod);
            squaresPlusN.disjuctionRight(squares);
            int[] values = squaresPlusN.values;
            int[] xCandidatesMod = new int[2*mod];
            int j = 0;
            for (int i = 0; i < mod; i++)
            {
                if (values[2*i] >= 0)
                {
                    xCandidatesMod[j++] = values[2*i];
                    // do not add the same value twice; can only happen with 0 and mod2Pow
                    if (Math.abs(values[2*i] - values[2*i+1]) != mod)
                        xCandidatesMod[j++] = values[2*i+1];
                }
            }
            xCandidatesMod[j] = -1;
        }
        System.out.println(" array : " + (System.currentTimeMillis() - time));

        time = System.currentTimeMillis();
        int nMod = mod;
        QuadraticDiophantineModBit squares = new QuadraticDiophantineModBit(nMod);

        for (int i = 0; i < limit; i++) {

//            long nMod = SquaresModArray.mod2Pow(-i, mod2Pow);
            int[] xCandidatesMod = squares.xArray(nMod);
            nMod--;
            if (nMod < 0)
                nMod += mod;

//            if (xCandidatesMod.length > 0)
            found = true;
        }
//        if (found)
        System.out.println(" bit   : " + (System.currentTimeMillis() - time));

        time = System.currentTimeMillis();
        nMod = mod;
        SquaresModBitPre squares2 = new SquaresModBitPre(mod);

        for (int i = 0; i < limit; i++) {

//            long nMod = SquaresModArray.mod2Pow(-i, mod2Pow);
            int[] xCandidatesMod = squares2.x(nMod);
            nMod--;
            if (nMod < 0)
                nMod += mod;

//            if (xCandidatesMod.length > 0)
            found = true;
        }
//        if (found)
        System.out.println(" pre   : " + (System.currentTimeMillis() - time));

        time = System.currentTimeMillis();

        int s=1;
        int sum = 0;
        for (int k = 0; k < limit; k++) {
            for (int j =0; j<mod; j++)
            s = k*k -j;
            sum += s;
        }
        if(sum != 0)
            System.out.println(" raw   : " + (System.currentTimeMillis() - time));

    }
}
