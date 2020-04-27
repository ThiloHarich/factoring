package factoring.math;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multiset;
import factoring.sieve.triple.TripleLookupSieve;
import org.apache.commons.math3.util.Pair;

import java.util.*;

import static org.junit.Assert.assertEquals;

public class SquareFinder {
    List<Pair<Multiset<Integer>, Multiset<Integer>>> factors;
    List<Integer> factorBase;
    Map<Integer, Integer> factorIndex = new HashMap<>();
    List<Column> matrix = new ArrayList<>();
    long n;

    public SquareFinder(List<Integer> factorBase, long n) {
        this.factorBase = factorBase;
        for (int i = 0; i < factorBase.size(); i++) {
            factorIndex.put(factorBase.get(i), i);
        }
        factors = new ArrayList<>();
        this.n = n;
    }

    public void addFactors(Multiset<Integer> left, Multiset<Integer> right){
        factors.add(Pair.create(left, right));
    }

    public List<Column> initMatrix(){
        int columnIndex = 0;
        for (Pair<Multiset<Integer>, Multiset<Integer>> factorsOfRel: factors) {
            final BitSet factorsMod2 = new BitSet();
            fillMatrix(factorsMod2, factorsOfRel.getFirst(), 0);
            fillMatrix(factorsMod2, factorsOfRel.getSecond(), factorBase.size());
            final Column col = Column.of(factorsMod2, columnIndex);
            matrix.add(col);
            columnIndex++;
            checkForSquare(col);
        }
        return matrix;
    }

    private void fillMatrix(BitSet factorsMod2, Multiset<Integer> factors, int offset) {
        for (Multiset.Entry<Integer> entry : factors.entrySet()) {
            int factor = entry.getElement();
            if (entry.getCount() % 2 == 1) {
                final Integer index = factorIndex.get(factor);
                factorsMod2.set(offset + index);
            }
        }
    }

    public long doGaussElimination(List<Column> matrix){
//        List<Column> matrix = matrixOrig.stream().collect(Collectors.toList());
        int rowCount = matrix.get(0).entries.size();
        int pivotCount = 0;
        for (int i = 0; i < rowCount; i++) {
            Column pivot = null;
            int pivotIndex = i-1;
            for (; pivotIndex < matrix.size()-1 && pivot == null;) {
                if (matrix.get(++pivotIndex).entries.get(i)) {
                    pivot = matrix.get(pivotIndex);
                }
            }
            if (pivot != null) {
                // exchange
                matrix.set(pivotIndex, matrix.get(i));
                matrix.set(i, pivot);
                for (int col = pivotIndex+1; col < matrix.size(); col++) {
                    if (matrix.get(col).entries.get(i)) {
                    	matrix.get(col).xor(pivot);
                    }
                }
            }
            print(matrix);
            checkForDuplicatedColumn(matrix);
        }
        for (int i = pivotCount +1; i < matrix.size(); i++){
            final long factor = checkForSquare(matrix.get(i));
            if (factor >= 1)
                return factor;
        }
        return -1;
    }

    public List<Column> reduceMatrix(List<Column> matrix) {
        int [] colCounts = new int [2 * factorBase.size()];

        List<List<Integer>> rows4Columns = new ArrayList<>();
        for (int l = 0; l < 2 * factorBase.size(); l++) {
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
                    matrix.get(index).xor(matrix.get(index0));
//                    checkForSquare(matrix.get(index));
                }
                if (rows4Column.size() > 0)
                    deleteRowsSet.add(rows4Column.get(0));
            }
        }
        System.out.print("rows to delete : ");
        deleteRowsSet.stream().map(i -> i + ",").forEach(System.out::print);
        System.out.println();
        List<Column> reducedMatrix = new ArrayList<>();
        for (int row = 0; row < matrix.size(); row++) {
            if (!deleteRowsSet.contains(row))
                reducedMatrix.add(matrix.get(row));
        }
        return reducedMatrix;
    }

    private void print(List<Column> matrix, int[] colCounts, List<List<Integer>> rows4Columns) {
        for (int j = 0; j < matrix.size(); j++) {
            BitSet row = matrix.get(j).entries;
            String line = String.format("%02d %02d", j, matrix.get(j).id);
            for (int l = 0; l < 2 * factorBase.size(); l++) {
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

    private void print(List<Column> matrix) {
    	System.out.println();
        for (int j = 0; j < matrix.size(); j++) {
            BitSet row = matrix.get(j).entries;
            String line = String.format("%02d %02d", j, matrix.get(j).id);
            line += String.format(" %02d", matrix.get(j).id);
            for (int l = 0; l < 2 * factorBase.size(); l++) {
                final boolean rowContainsPrime = row.get(l);
                line += rowContainsPrime ? "o|" : " |";
            }
            line += matrix.get(j).toString(factors, n);
            System.out.println(line);
        }
    }

    private long checkForDuplicatedColumn(List<Column> matrix) {
        for (Column col1 : matrix) {
            for (Column col2 : matrix) {
                if (col1.id != col2.id && col1.entries.equals(col2.entries)){
                    System.out.println("two identical lines found");
                    System.out.println(col1.toString(factors, n));
                    System.out.println(col2.toString(factors, n));
                    Column xor = col1.copy();
                    xor.xor(col2);
                    if (xor.gcd(factors, n) > 1)
                    	return xor.gcd(factors, n);
                }
            }
        }
        return -1;
    }
    private long checkForSquare(Column col) {
        if (col.entries.cardinality() == 0){
            System.out.println("candidate found");
            long prodLeft = 1;
            long prodRight = 1;
            BitSet bs = col.columns;
            for (int column  = bs.nextSetBit(0); column >= 0; column = bs.nextSetBit(column+1)) {
                Pair<Multiset<Integer>, Multiset<Integer>> factorsO = factors.get(column);
                prodLeft = TripleLookupSieve.multiplyMod(factorsO.getFirst(), n);
                prodRight = TripleLookupSieve.multiplyMod(factorsO.getSecond(), n);
                final long sqrtLeft = (long) Math.sqrt(prodLeft);
                assertEquals(prodLeft,sqrtLeft*sqrtLeft);
                final long sqrtRight = (long) Math.sqrt(prodRight);
                assertEquals(prodRight,sqrtRight*sqrtRight);
                final long gcd = PrimeMath.gcd(sqrtLeft + sqrtRight, n);
                if (gcd > 1 && PrimeMath.gcd(Math.abs(sqrtLeft - sqrtRight),n) > 0){
                    System.out.println("factor found : " + gcd);
                    return gcd;

                }
                if (column == Integer.MAX_VALUE) {
                    break; // or (column+1) would overflow
                }
            }
        }
        return -1;
    }

    private void print(int[] colCounts) {
        System.out.print(" ");
        for (int i = 0; i < colCounts.length; i++) {
            System.out.print(String.format("%02d", i));
        }
        System.out.println();
        System.out.print(" ");

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
