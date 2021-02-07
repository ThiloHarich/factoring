package factoring.math;

import factoring.trial.TDiv31Barrett;

import java.math.BigInteger;
import java.util.*;

import static factoring.math.PrimeMath.mod;
import static org.junit.Assert.assertEquals;

public class SquareFinder3 {
    public static final int NO_FACTOR = -1;
    //    List<Pair<Multiset<Integer>, Multiset<Integer>>> factors;
    List<Relation> relations;
    int[] factorBase;
    //    Map<Integer, Integer> factorIndex = new HashMap<>();
    public List<Row> matrix = new ArrayList<>();
    long n;
    boolean print = true;
    boolean check = false;
    int columnIndex = 0;
    static int factorisationWords;


    public SquareFinder3(int[] factorBase, long n) {
        this.factorBase = factorBase;
        factorisationWords = (int) Math.ceil((factorBase.length * .125)); // length * 8 / 64 = length / 8 = length >> 3
        relations = new ArrayList<>();
        this.n = n;
    }

    public Relation addFactors(int lower, int higher, int right, long[] exponentsMod2, int maxFactorIndex) {
        // TODO get rid of Relation, Column, matrix classes
        // TODO delete, since it can not happen
        if (check) {
            long left = ((long)lower) * ((long)higher);
            for (Relation relation :
                    relations) {
                if (((long) relation.lower) * ((long) relation.higher) == left || (relation.lower == lower && relation.higher == higher)) {
                    System.out.println("left = " + left + " already exists. Should not happen");
                    return null;
                }
                if (relation.right == right) {
                    System.out.println("right = " + right + " for left = " + left + " already exists. Should not happen");
                    return null;
                }
            }
        }
        // this is all Overhead
        final Relation rel = getRelation(lower, higher, right, exponentsMod2, maxFactorIndex);
        columnIndex++;
        return rel;
    }

    private Relation getRelation(int lower, int higher, int right, long[] exponentsMod2, int maxFactorIndex) {
        final Relation rel = new Relation(lower, higher, right, maxFactorIndex);
        relations.add(rel);

//        long[] exponents = {exponentsMod2};
        BitSet exponentsSet = BitSet.valueOf(exponentsMod2);
        // column can be an array of longs columIndex, exponents mod 2
        final Row col = Row.of(exponentsSet, columnIndex);
        matrix.add(col);
        return rel;
    }


    // we might have
    // 8 bits per prime for the first 8 primes   -> 2^256
    // 4 bits per prime for the next 16 primes   -> 23^16  > 2^72 = 2^64
    // 2 bits per prime for all other  primes    -> 97^4   > 2^26
//    private long fillExponents(long number, long exponents, long[] factorCounts, long exponentsMod2) {
//        // we iterate over the factors dividing left, deterime how often they divide left
//        int primeIndex = -1;
//        // handle negative right numbers
//        if ((exponents & 1) == 1){
//            number = - number;
//            exponentsMod2 ^= 1;
//            factorCounts[0] += 1;
//            exponents ^= 1;
//        }
//        if (number == 1)
//            return exponentsMod2;
//        do {
//            int primeIndexDiff = Long.numberOfTrailingZeros(exponents);
//            primeIndex += primeIndexDiff + 1;
//            int prime = factorBase[primeIndex];
//            long exponent = 0;
//            long leftDivPrime;
//            // TODO for smaller numbers 52 or 31 bit we do not need the division
//            // TODO we might delay calculation of factorCounts
//            while ((leftDivPrime = (number / prime)) * prime == number) {
//                number = leftDivPrime;
//                exponent++;
//            }
//            if ((exponent & 1l) == 1l )
//                exponentsMod2 ^= (1l << primeIndex);
//            // 8 = 2^3 bits per exponent -> epxonent per prime < 256 -> number > p^256 >= 2^256
//            // TODO dynamic bits ?
//            final int longIndex = primeIndex >> 3;
//            final int indexInWord = (primeIndex - (longIndex << 3)) << 3;
//            factorCounts[longIndex] += (exponent << indexInWord);
//            exponents >>= (primeIndexDiff+1);
//        }while(number > 1);
//        return exponentsMod2;
//    }
//


    public long doGaussElimination(List<Row> matrix){
//        List<Column> matrix = matrixOrig.stream().collect(Collectors.toList());
        int rowCount = matrix.size();
        int pivotCount = 0;
        int colIndex = 0;
        for (int i = factorBase.length-1; i >= 0 ; i--) {
//            for (int i = 0; i < factorBase.length; i++) {
            Row pivot = null;
            int pivotIndex = colIndex-1;
            for (; pivotIndex < matrix.size()-1 && pivot == null;) {
                if (matrix.get(++pivotIndex).entries.get(i)) {
                    pivot = matrix.get(pivotIndex);
                }
            }
            if (pivot != null) {
                colIndex = movePivot(matrix, colIndex, i, pivot, pivotIndex);
            }
            print(matrix);
//            checkForDuplicatedColumn(matrix);
        }
        for (int i = colIndex; i < matrix.size(); i++){
            final long factor = checkForSquare(matrix.get(i));
            if (factor >= 1) {
                final int s = i - colIndex;
//                if(s > 11)
                if (print)
                    System.out.println("factor " + factor +  " found after try " + s);
                return factor;

            }
        }
        if (print) System.out.println("factor not found after try " + (matrix.size() - colIndex));
        return NO_FACTOR;
    }

    private int movePivot(List<Row> matrix, int colIndex, int i, Row pivot, int pivotIndex) {
        // exchange
        matrix.set(pivotIndex, matrix.get(colIndex));
        matrix.set(colIndex, pivot);
        colIndex++;
        xor(matrix, i, pivot, pivotIndex);
        return colIndex;
    }

    private void xor(List<Row> matrix, int i, Row pivot, int pivotIndex) {
        for (int col = pivotIndex +1; col < matrix.size(); col++) {
            if (matrix.get(col).entries.get(i)) {
                matrix.get(col).xor(pivot);
            }
        }
    }

    public List<Row> reduceMatrix(List<Row> matrix) {
        int [] colCounts = new int [factorBase.length];

        List<List<Integer>> rows4Columns = new ArrayList<>();
        for (int l = 0; l < factorBase.length; l++) {
            rows4Columns.add(new ArrayList<>());
        }
        print(matrix, colCounts, rows4Columns);
        print(colCounts);
        HashSet<Integer> deleteRowsSet = new HashSet<>();
        for (int i = 0; i < rows4Columns.size(); i++) {
            List<Integer> rows4Column = rows4Columns.get(i);
            //            for (List<Integer> rows4Column : rows4Columns) {
            if (rows4Column.size() <= 2) {
                if (rows4Column.size() < 2) {
                    if (rows4Column.size() == 1)
                        deleteRowsSet.add(rows4Column.get(0));
                }
                else
//                if (rows4Column.size() == 2)
                {
                    final Integer index = rows4Column.get(1);
                    final Integer index0 = rows4Column.get(0);
                    if (deleteRowsSet.contains(index)) {
                        if (!deleteRowsSet.contains(index0)) {
                            deleteRowsSet.add(index0);
                        }
                    }
                    else{
                        if (!deleteRowsSet.contains(index)){
                            matrix.get(index).xor(matrix.get(index0), columnIndex);
                            deleteRowsSet.add(index0);
                        }
                        else{
                            deleteRowsSet.add(index);
                        }
                    }
                    // If all entries are 0 we are done. Stop early here
//                    checkForSquare(matrix.get(index));
                }
                colCounts[i] = 0;
            }
        }
        if (print) {
            System.out.print("rows to delete : ");
            deleteRowsSet.stream().map(i -> i + ",").forEach(System.out::print);
            System.out.println();
        }
        List<Row> reducedMatrix = new ArrayList<>();
        for (int row = 0; row < matrix.size(); row++) {
            if (!deleteRowsSet.contains(row))
                reducedMatrix.add(matrix.get(row));
        }
        return reducedMatrix;
    }

    private void print(List<Row> matrix, int[] colCounts, List<List<Integer>> rows4Columns) {
        if (!print)
            return;
        for (int j = 0; j < matrix.size(); j++) {
            BitSet row = matrix.get(j).entries;
            String line = String.format("%02d %02d", j, matrix.get(j).id);
            for (int l = 0; l < factorBase.length; l++) {
                final boolean rowContainsPrime = row.get(l);
                line += rowContainsPrime ? "o|" : " |";
                if (rowContainsPrime) {
                    rows4Columns.get(l).add(j);
                    colCounts[l]++;
                }
            }
//            line += matrix.get(j).toString(factors, n);
            System.out.println(line);
        }
    }

    private void print(List<Row> matrix) {
        if (!print)
            return;
        System.out.println();
        for (int j = 0; j < matrix.size(); j++) {
            BitSet row = matrix.get(j).entries;
            String line = String.format("%02d %02d", j, matrix.get(j).id);
            line += String.format(" %02d", matrix.get(j).id);
            for (int l = 0; l < factorBase.length; l++) {
                final boolean rowContainsPrime = row.get(l);
                line += rowContainsPrime ? "o|" : " |";
            }
//            line += matrix.get(j).toString(factors, n);
            System.out.println(line);
        }
    }

    private long checkForSquare(Row col) {
        long [] leftFactorCounts = new long [factorisationWords];
        long [] rightFactorCounts = new long [factorisationWords];
        if (col.entries.cardinality() == 0){
            if (print)
//                if (true)
                System.out.println("candidate found : ");
            long x = 1l;
            long y = 1l;
            BitSet bs = col.columns;
            long oddExponents = 0;
            // multiply together all the factors from the columns mod n
            for (int column  = bs.nextSetBit(0); column >= 0; column = bs.nextSetBit(column+1)) {
                Relation rel = relations.get(column);
                if (print) System.out.println(" rel [" + column + "] : " + rel.toString(factorBase));

                // we need to count all the left factors for each column
                long[] leftFactrorisation = rel.getLeftFactrorisation();
                long[] rightFactrorisation = rel.getRightFactrorisation();
                if (print) {
                    System.out.println(" left exponents  : " + printExponents(leftFactrorisation, factorBase));
                    System.out.println(" right exponents : " + printExponents(rightFactrorisation, factorBase));
                }
                if (check) {
                    if (((long) rel.lower) * rel.higher != multiply(leftFactrorisation, factorBase)) {
                        System.out.println("left wrong");
                    }
                    final long right = multiply(rightFactrorisation, factorBase);
                    if (rel.right != right) {
                        System.out.println("right wrong");
                    }
                }
                for (int i = 0; i < leftFactorCounts.length; i++) {
                    leftFactorCounts[i] += leftFactrorisation[i];
                    rightFactorCounts[i] += rightFactrorisation[i];
                }
            }
            if (print) {
                System.out.println(printExponents(leftFactorCounts, factorBase));
                System.out.println(printExponents(rightFactorCounts, factorBase));
            }
            for (int i = 0; i < factorBase.length; i++) {
                int prime = factorBase[i];
//                int bucket = i >> 3;
//                final int indexInLong = i - bucket;
//                int count = (int) ((leftFactorCounts[bucket] & 255l << indexInLong) >> indexInLong);
                int wordIndex = i >> 3;
                final int indexInWord = (i - (wordIndex << 3)) << 3;
                long countShifted = leftFactorCounts[wordIndex] & (255l << indexInWord);
                int count = (int) (countShifted >> indexInWord);
                long countShiftedRight = rightFactorCounts[wordIndex] & (255l << indexInWord);
                int countRight = (int) (countShiftedRight >> indexInWord);
                // multiply both sides with the prime, if it is odd
                if ((count & 1) == 1) {
//                    if ((countRight & 1) == 0)
//                        System.out.println("left != right");
                    x = PrimeMath.mod(x * prime, n);
                    y = PrimeMath.mod(y * prime, n);
                }
                // TODO calculate prime^count mod n faster
                for (int j = 0; j < ((count+1) >> 1); j++) {
                    x = PrimeMath.mod(x * prime, n);
                }
                for (int j = 0; j < ((countRight+1) >> 1); j++) {
                    y = PrimeMath.mod(y * prime, n);
                }
            }
            BigInteger nBig = BigInteger.valueOf(n);
            if(print){
//                long x = factorsLeft.entrySet().stream().map(e -> sqrt(e)).reduce(1l, (a, b) -> (a * b + n) % n);
//                long y = factorsRight.entrySet().stream().map(e -> sqrt(e)).reduce(1l, (a, b) -> (a * b + n) % n);
                System.out.println("x  :" +x);
                System.out.println("y  :" +y);
                System.out.println("x+y:" +(x+y)   % n);
                System.out.println("x-y:" +(x-y+n) % n);
            }
//            final long sqrtLeft = sqrtMod(factorsLeft, n);
//            final long sqrtRight = sqrtMod(factorsRight, n);
            final long xPlusy = (x + y) % n;
            final long gcd = PrimeMath.gcd(xPlusy, n);
            final long xMiny = (x - y + n ) % n;
            final long diff = mod(xMiny, n);
            final long gcd2 = PrimeMath.gcd(diff, n);


            if (gcd > 1 && gcd < n){
                if (print) System.out.println("factor found : " + gcd);
                return gcd;
            }
            if (gcd2 > 1 && gcd2 < n){
                if (print) System.out.println("factor found : " + gcd2);
                return gcd2;
            }

        }
        return -1;
    }



    private void print(int[] colCounts) {
        if (!print)
            return;
        System.out.print("     ");
        for (int i = 0; i < colCounts.length; i++) {
            System.out.print(String.format("%02d", i));
        }
        System.out.println();
        System.out.print("     ");

        int colEmpty = 0;
        for (int i = 0; i < colCounts.length; i++) {
            if (colCounts[i] <= 2) {
                System.out.print("_" +String.format("%01d", colCounts[i]));
                colEmpty += colCounts[i] == 0 ? 1 : 0;
            }
            else
                System.out.print(String.format("%02d", colCounts[i]));
        }
        System.out.println();
        System.out.println("colums non empty : " + (colCounts.length - colEmpty));
    }

    public static String printExponents(long[] factorCounts, int [] factorBase) {
        final Object apply = apply(factorCounts,
                factorBase,
                "",
                (s, prime, count) -> (String)s + prime + "^" + count + " ");
        return (String) apply;
//        String s = "";
//        if (factorCounts == null)
//            return s;
//        int lastPos = 0;
//        for (int i = 0; i < factorBase.length; i++) {
//            int prime = factorBase[i];
//            int wordIndex = i >> 3;
//            final int indexInWord = (i - (wordIndex << 3)) << 3;
//            long countShifted = factorCounts[wordIndex] & (255l << indexInWord);
//            int count = (int) (countShifted >> indexInWord);
//            if (count != 0)
//                s += prime + "^" + count + " ";
//        }
//        return s;
    }

    public static long multiply(long[] factorCounts, int [] factorBase) {
        final Object apply = apply(factorCounts,
                factorBase,
                1L,
                (s, prime, count) -> (Long)s * ((long)Math.pow(prime, count)));
        return (long) apply;
    }

    public static Object apply(long[] factorCounts, int [] factorBase, Object startVal, TriFunction func) {
        Object s = startVal;
        if (factorCounts == null)
            return s;
        int lastPos = 0;
        for (int i = 0; i < factorBase.length; i++) {
            int prime = factorBase[i];
            int wordIndex = i >> 3;
            final int indexInWord = (i - (wordIndex << 3)) << 3;
            long countShifted = factorCounts[wordIndex] & (255l << indexInWord);
            int count = (int) (countShifted >> indexInWord);
            if (count != 0)
                s = func.apply(s, prime, count);
        }
        return s;
    }


    public static class Relation {
        int lower;
        int higher;
        int right;
        // other then in the classic Quadratic sieve a relation has no sqares on the left side. We can combine non
        // squares on the left side. To calculate the square root, we store the exponents of the factors of the factor
        // base. Each count of a prine is represented by some bits in the long arrray. So we can sum up all
        // exponents by adding the long array, as long as we have reseved enough bits for each prime such that they
        // do not interfer with each other.
        // for one prime factor p we can store an exponents up to 256 = 2^8 -> number > p^256 >= 2^256
        // we do this by shift <<3
        // TODO for bigger primes we do not need 8 bits
        long[] leftFactorCounts;
        long[] rightFactorCounts;
        int maxFactorIndex;

        static TDiv31Barrett smallFactoriser = new TDiv31Barrett();

        public Relation(int lower, int higher, int right, int maxFactorIndex) {
            this.lower = lower;
            this.higher = higher;
            this.right = right;
            this.maxFactorIndex = maxFactorIndex;
        }

        public String toString(int [] factorBase) {
            return "Relation{" +
                    "left=" + ((long)lower) * ((long)higher) +
                    ", right=" + right +
                    ", leftFactorCounts=" + printExponents(leftFactorCounts, factorBase) +
                    ", rightFactorCounts=" + printExponents(rightFactorCounts, factorBase) +
                    '}';
        }

        public long[] getLeftFactrorisation() {
            if (leftFactorCounts == null) {
                leftFactorCounts = new long[factorisationWords];

                // TODO we might divide out small factors on lower,higher,right and as soon as lower*higher*right
                // fit in long, int find the facors there
                smallFactoriser.fillFactorCounts(lower, maxFactorIndex, leftFactorCounts);
                smallFactoriser.fillFactorCounts(higher, maxFactorIndex, leftFactorCounts);
//
//                fillFactorCounts(lower, leftFactorCounts);
//                fillFactorCounts(higher, leftFactorCounts);
            }
            return leftFactorCounts;
        }
        public long[] getRightFactrorisation() {
            if (rightFactorCounts == null) {
                rightFactorCounts = new long[factorisationWords];
                // set sign
                if (right < 0)
                    rightFactorCounts[0] = 1;
                smallFactoriser.fillFactorCounts(Math.abs(right), maxFactorIndex, rightFactorCounts);
            }
            return rightFactorCounts;
        }
//        private void fillFactorCounts(int number, long[] lowerFactorCounts) {
//            // TODO for left side we do not need the max factor, but storing all factorizations uses a lot of memory
//            int numberRemaining = number;
//            int wordIndex = lowerFactorCounts.length -1;
//            // we reduce the number until we have the hole factorisation
//            while (numberRemaining > SmoothNumbersSieve.factorisationBound) {
//                // TODO maybe directy store the prime for Performance
//                int maxFactorIndex = SmoothNumbersSieve.maxFactorIndex[numberRemaining];
//                int factor = SmoothNumbersSieve.primes[maxFactorIndex];
//                // add max factor
//                wordIndex = maxFactorIndex >> 3;
//                // TODO just use maxFactorIndex
//                final int indexInWord = (maxFactorIndex - (wordIndex << 3)) << 3;
//                lowerFactorCounts[wordIndex] += (1l << indexInWord);
//                // TODO use Barret division
//                numberRemaining = numberRemaining / factor;
//            }
//            // copy over the prime factorization ONLY up to max factor only. wordIndex <= lowerFactorCounts.length
//            for (int i = 0; i <= wordIndex; i++) {
//                final int index = numberRemaining * SmoothNumbersSieve.factorisationWords + i;
//                lowerFactorCounts[i] += SmoothNumbersSieve.factorisation[index];
//            }
//        }
//
    }
}
