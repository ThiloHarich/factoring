package de.harich.thilo.factoring.hart.util;

import java.util.HashSet;
import java.util.Set;

public class MultiplierFactorizationCount implements Comparable<MultiplierFactorizationCount>{
    public MultiplierFactorizationCount(int multiplier) {
        this.multiplier = multiplier;
        this.count = 0;
    }
    public MultiplierFactorizationCount(int multiplier, int count) {
        this.multiplier = multiplier;
        this.count = count;
    }
    public int multiplier;
    public int count;
    public Set<Long> factorisations = new HashSet<>();

    public void addFactorisation(long semiprime){
        factorisations.add(semiprime);
    }

    public void removeFactorisation(Set<Long> factors) {
        this.factorisations.removeAll(factors);
    }

    @Override
    public int compareTo(MultiplierFactorizationCount o) {
//        return o.factorisations.size() - factorisations.size();
        if (o.count == count)
            return multiplier - o.multiplier;
        return o.count - count;
    }
    @Override
    public String toString() {
//        return "[multiplier=" + multiplier + ", hits=" + factorisations.size() + "]";
        return "[multiplier=" + multiplier + ", hits=" + count + "]";
    }
}