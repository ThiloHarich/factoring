package factoring.math;

import static factoring.math.PrimeMath.mod;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class SquareFinder2 {
    public static final int NO_FACTOR = -1;
    List<Pair<Multiset<Integer>, Multiset<Integer>>> factors;
    List<Pair<Long,Integer>> relations;
    List<Integer> factorBase;
    Map<Integer, Integer> factorIndex = new HashMap<>();
    List<Row> matrix = new ArrayList<>();
    long n;
    boolean print = false;
    boolean check = false;
    int columnIndex = 0;

    public SquareFinder2(List<Integer> factorBase, long n) {
        this.factorBase = factorBase;
        for (int i = 0; i < factorBase.size(); i++) {
            factorIndex.put(factorBase.get(i), i);
        }
        factors = new ArrayList<>();
        relations = new ArrayList<>();
        this.n = n;
    }

    /**
     * tries to add a relation left and right to the set of relations we like to solve.
     *
     * @param left
     * @param right
     * @return true if the relation does not already exists.
     */
    public boolean addFactors(Multiset<Integer> left, Multiset<Integer> right){
        for (Pair relation: factors) {
            if (relation.getFirst().equals(left)){
                if (check) System.out.println("left = " + left + " already exists");
                return false;
            }
            if (relation.getSecond().equals(right)){
                if (check) System.out.println("right for left " + left + " already exists. Should not happen");
                return false;
            }
        }
        factors.add(Pair.create(left, right));
        return true;
    }
    public boolean addFactors(BitSet leftExponents, int left, BitSet rightExponents, int right) {
        for (Pair relation :
                relations) {
            if (relation.getFirst().equals(left)){
                if (check) System.out.println("left = " + left + " already exists");
                return false;
            }
            if (relation.getSecond().equals(right)){
                if (check) System.out.println("right = " + right + " for left = " + left + " already exists. Should not happen");
                return false;
            }
        }
        relations.add(new Pair(left, right));
        leftExponents.xor(rightExponents);
        final Row col = Row.of(leftExponents, columnIndex);
        matrix.add(col);
        columnIndex++;
        return true;
    }


    public List<Row> initMatrix(){
        for (Pair<Multiset<Integer>, Multiset<Integer>> factorsOfRel: factors) {
            final BitSet factorsMod2 = new BitSet();
            fillMatrix(factorsMod2, factorsOfRel.getFirst());
            fillMatrix(factorsMod2, factorsOfRel.getSecond());
            final Row col = Row.of(factorsMod2, columnIndex);
            matrix.add(col);
            columnIndex++;
            long factor = checkForSquare(col);
            if (factor > 1 && print )
                System.out.println("factor found during init without matrix step");
        }
        return matrix;
    }

    private void fillMatrix(BitSet factorsMod2, Multiset<Integer> factors) {
        for (Multiset.Entry<Integer> entry : factors.entrySet()) {
            int factor = entry.getElement();
            if (entry.getCount() % 2 == 1) {
                final Integer index = factorIndex.get(factor);
                if(index == null)
                    System.out.println();
                factorsMod2.flip(index);
            }
        }
    }

    public long doGaussElimination(List<Row> matrix){
//        List<Column> matrix = matrixOrig.stream().collect(Collectors.toList());
        int rowCount = matrix.size();
        int pivotCount = 0;
        int colIndex = 0;
        for (int i = 0; i < factorBase.size(); i++) {
            Row pivot = null;
            int pivotIndex = colIndex-1;
            for (; pivotIndex < matrix.size()-1 && pivot == null;) {
                if (matrix.get(++pivotIndex).entries.get(i)) {
                    pivot = matrix.get(pivotIndex);
                }
            }
            if (pivot != null) {
                // exchange
                matrix.set(pivotIndex, matrix.get(colIndex));
                matrix.set(colIndex, pivot);
                colIndex++;
                for (int col = pivotIndex+1; col < matrix.size(); col++) {
                    if (matrix.get(col).entries.get(i)) {
                        matrix.get(col).xor(pivot, columnIndex);
                    }
                }
            }
            print(matrix);
//            checkForDuplicatedColumn(matrix);
        }
        for (int i = colIndex; i < matrix.size(); i++){
            final long factor = checkForSquare(matrix.get(i));
            if (factor >= 1) {
                final int s = i - colIndex;
//                if(s > 11)
                if (print) System.out.println("factor found after try " + s);
                return factor;
            }
        }
        if (print) System.out.println("factor not found after try " + (matrix.size() - colIndex));
        return NO_FACTOR;
    }

    public List<Row> reduceMatrix(List<Row> matrix) {
        int [] colCounts = new int [factorBase.size()];

        List<List<Integer>> rows4Columns = new ArrayList<>();
        for (int l = 0; l < factorBase.size(); l++) {
            rows4Columns.add(new ArrayList<>());
        }
        print(matrix, colCounts, rows4Columns);
        print(colCounts);
        HashSet<Integer> deleteRowsSet = new HashSet<>();
        for (List<Integer> rows4Column : rows4Columns) {
            if (rows4Column.size() <= 2) {
                if (rows4Column.size() == 2) {
                    final Integer index = rows4Column.get(1);
                    final Integer index0 = rows4Column.get(0);
                    matrix.get(index).xor(matrix.get(index0), columnIndex);
                    // If all entries are 0 we are done. Stop early here
                    checkForSquare(matrix.get(index));
                }
                if (rows4Column.size() > 0)
                    deleteRowsSet.add(rows4Column.get(0));
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
            for (int l = 0; l < factorBase.size(); l++) {
                final boolean rowContainsPrime = row.get(l);
                line += rowContainsPrime ? "o|" : " |";
                if (rowContainsPrime) {
                    rows4Columns.get(l).add(j);
                    colCounts[l]++;
                }
            }
            line += matrix.get(j).toString(factors, n);
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
            for (int l = 0; l < factorBase.size(); l++) {
                final boolean rowContainsPrime = row.get(l);
                line += rowContainsPrime ? "o|" : " |";
            }
            line += matrix.get(j).toString(factors, n);
            System.out.println(line);
        }
    }

    private long checkForDuplicatedColumn(List<Row> matrix) {
        for (Row col1 : matrix) {
            for (Row col2 : matrix) {
                if (col1.id != col2.id && col1.entries.equals(col2.entries)){
                    System.out.println("two identical lines found");
                    System.out.println(col1.toString(factors, n));
                    System.out.println(col2.toString(factors, n));
                    Row xor = col1.copy();
                    xor.xor(col2, columnIndex);
                    if (xor.gcd(factors, n) > 1)
                        return xor.gcd(factors, n);
                }
            }
        }
        return -1;
    }
    private long checkForSquare(Row col) {
        if (col.entries.cardinality() == 0){
            if (print)
//                if (true)
                System.out.println("candidate found : " + col.toString(factors, n));
            Multiset<Integer> factorsLeft = HashMultiset.create();
            Multiset<Integer> factorsRight = HashMultiset.create();
            BitSet bs = col.columns;
            // multiply together all the factors from the columns mod n
            for (int column  = bs.nextSetBit(0); column >= 0; column = bs.nextSetBit(column+1)) {
                Pair<Multiset<Integer>, Multiset<Integer>> factorsO = factors.get(column);
                factorsLeft.addAll(factorsO.getFirst());
                factorsRight.addAll(factorsO.getSecond());
                if (print) System.out.println("factors " + factorsO.getFirst() + " right " + factorsO.getSecond());
            }
            if(print){
                long x = factorsLeft.entrySet().stream().map(e -> sqrt(e)).reduce(1l, (a, b) -> (a * b + n) % n);
                long y = factorsRight.entrySet().stream().map(e -> sqrt(e)).reduce(1l, (a, b) -> (a * b + n) % n);
                System.out.println("x  :" +x);
                System.out.println("y  :" +y);
                System.out.println("x+y:" +(x+y)   % n);
                System.out.println("x-y:" +(x-y+n) % n);
            }
            final long sqrtLeft = sqrtMod(factorsLeft, n);
            final long sqrtRight = sqrtMod(factorsRight, n);
            final long xPlusy = sqrtLeft + sqrtRight;
            final long gcd = PrimeMath.gcd(xPlusy, n);
            final long xMiny = sqrtLeft - sqrtRight;
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

    private long sqrt(Multiset.Entry<Integer> e) {
        // TODO align with Column.toString() + sqrtMod
        // TODO perfromance Math.pow <-> multiply
        long sqrt =  (long) Math.pow(e.getElement(), e.getCount() / 2) % n;
        sqrt *= e.getCount() % 2 == 1 ? e.getElement() : 1;
        return sqrt;
    }

    private long sqrtMod(Multiset<Integer> factors, long n) {
        long prod = 1;
        for (Multiset.Entry<Integer> entry : factors.entrySet()) {
//            Assert.assertEquals(0, entry.getCount() & 1);
            final int exp = entry.getCount() / 2;
            for (int i = 0; i < exp; i++) {
                prod = mod(prod * entry.getElement(), n);
            }
            // if exponent is odd on one side (but overall exponent) is even, multiply booth sides by the factor and we
            // have a relation
            if (1 == (entry.getCount() & 1)){
                prod = mod(prod * entry.getElement(), n);
            }
        }
        return prod;
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


}
