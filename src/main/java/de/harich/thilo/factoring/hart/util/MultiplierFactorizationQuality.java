package de.harich.thilo.factoring.hart.util;

import java.util.HashSet;
import java.util.Set;

public class MultiplierFactorizationQuality implements Comparable<MultiplierFactorizationQuality>{
    public MultiplierFactorizationQuality(int multiplier) {
        this.multiplier = multiplier;
    }
    public MultiplierFactorizationQuality(int multiplier, double quality) {
        this.multiplier = multiplier;
        this.quality = quality;
    }
    public int multiplier;
    public double quality;
    public Set<Long> factorisations = new HashSet<>();

    public void addFactorisation(long semiprime){
        factorisations.add(semiprime);
    }

    public void removeFactorisation(Set<Long> factors) {
        this.factorisations.removeAll(factors);
    }

    @Override
    public int compareTo(MultiplierFactorizationQuality o) {
//        return o.factorisations.size() - factorisations.size();
        return (int) Math.signum(o.quality - quality);
    }
    @Override
    public String toString() {
//        return "[multiplier=" + multiplier + ", hits=" + factorisations.size() + "]";
        return "[multiplier=" + multiplier + ", quality =" + quality + "]";
    }
}