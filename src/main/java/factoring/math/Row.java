package factoring.math;

import java.math.BigInteger;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.Pair;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableSortedMultiset;
import com.google.common.collect.Multiset;

import factoring.sieve.triple.TripleLookupSieve;

public class Row {
        int id;
        BitSet entries;
        BitSet columns;

        public static Row of(BitSet entries, int columnIndex) {
            Row col = new Row();
            col.entries = entries;
            col.id = columnIndex;
            col.columns = new BitSet();
            col.columns.set(columnIndex);
            return col;
        }
        
        public Row copy (){
        	Row clone = new Row();
        	clone.columns = (BitSet) this.columns.clone();
        	clone.entries = (BitSet) this.entries.clone();
        	clone.id = this.id;
        	return clone;
        }
        
        public void xor(Row col) {
        	this.entries.xor(col.entries);
        	this.columns.xor(col.columns);
        }
        
        public long gcd (List<Pair<Multiset<Integer>, Multiset<Integer>>> factors, long n) {
            Multiset<Integer> factorsLeft = HashMultiset.create();
            Multiset<Integer> factorsRight = HashMultiset.create();
	        for (int column  = columns.nextSetBit(0); column >= 0; column = columns.nextSetBit(column+1)) {
	            Pair<Multiset<Integer>, Multiset<Integer>> newFactors = factors.get(column);
	            addFactors(factorsLeft, newFactors.getFirst());
	            addFactors(factorsRight, newFactors.getSecond());
	            // TODO combine multiply, sqrt, mod
	            if (column == Integer.MAX_VALUE) {
	                break; // or (column+1) would overflow
	            }	        
	        }

            BigInteger x = TripleLookupSieve.multiplyBig(factorsLeft).sqrt();
            final BigInteger nBig = BigInteger.valueOf(n);
			long xModN = x.mod(nBig).longValue();
            long xModN2 = TripleLookupSieve.multiplySqrtMod(factorsLeft, n);
            BigInteger y = TripleLookupSieve.multiplyBig(factorsRight).sqrt();
            long yModN = y.mod(nBig).longValue();
            if (xModN == yModN)
            	return -1;
            final long gcd = PrimeMath.gcd(xModN + yModN, n);
            if (gcd > 1 && gcd != n && (n / gcd) * gcd == n ){
                System.out.println("factor found : " + gcd);
                return gcd;
            }

	        return -1;	        
        }
        
        public String toString(List<Pair<Multiset<Integer>, Multiset<Integer>>> factors, long n) {
            Multiset<Integer> factorsLeft = HashMultiset.create();
            Multiset<Integer> factorsRight = HashMultiset.create();
            BitSet bs = columns;
            for (int column  = bs.nextSetBit(0); column >= 0; column = bs.nextSetBit(column+1)) {
                final Multiset<Integer> leftFactors = factors.get(column).getFirst();
                addFactors(factorsLeft, leftFactors);
                final Multiset<Integer> rightFactors = factors.get(column).getSecond();
                addFactors(factorsRight, rightFactors);
                if (column == Integer.MAX_VALUE) {
                    break; // or (column+1) would overflow
                }
            }

            String factorsL = ImmutableSortedMultiset.copyOf(factorsLeft).entrySet().stream().map(e -> e.getElement() + "^" + e.getCount()).collect(Collectors.joining("*"));
            String factorsR = ImmutableSortedMultiset.copyOf(factorsRight).entrySet().stream().map(e -> e.getElement() + "^" + e.getCount()).collect(Collectors.joining("*"));
            // only works if all exponents are even
            long prodLeft = factorsLeft.entrySet().stream().map(e -> (long)Math.pow(e.getElement(), e.getCount()/2)%n).reduce(1l, (a, b) -> (a * b) % n);
            long prodRight = factorsRight.entrySet().stream().map(e -> (long)Math.pow(e.getElement(), e.getCount()/2)%n).reduce(1l, (a, b) -> (a * b) % n);
//            assertEquals(prodLeft - n, prodRight);
            return "Column{" +
                    "id=" + id +
                    ", entries=" + entries +
                    ", columns=" + columns +
                    ", factors left =" + factorsL +
                    ", factors right =" + factorsR +
                    ", x % n =" + prodLeft +
                    ", y % n =" + prodRight +
                    '}';
        }

		public static void addFactors(Multiset<Integer> result, final Multiset<Integer> newFactors) {
			for (Multiset.Entry<Integer> entry :newFactors.entrySet()) {
			    if (!result.contains(entry.getElement())){
			        result.add(entry.getElement(), entry.getCount());
			    }
			    else{
			        result.setCount(entry.getElement(), entry.getCount() + result.count(entry.getElement()));
			    }
			}
		}
    }