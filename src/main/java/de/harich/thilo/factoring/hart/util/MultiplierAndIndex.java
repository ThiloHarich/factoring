package de.harich.thilo.factoring.hart.util;

import java.util.HashSet;
import java.util.Set;

public class MultiplierAndIndex{
    public MultiplierAndIndex(int multiplier) {
        this.multiplier = multiplier;
    }
    public MultiplierAndIndex(int multiplier, int index) {
        this.multiplier = multiplier;
        this.index = index;
    }
    public int multiplier;
    public int index;

    @Override
    public String toString() {
//        return "[multiplier=" + multiplier + ", hits=" + factorisations.size() + "]";
        return "[multiplier=" + multiplier + ", hits=" + index + "]";
    }
}