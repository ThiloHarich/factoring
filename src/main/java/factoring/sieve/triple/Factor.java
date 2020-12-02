package factoring.sieve.triple;

public class Factor {
    int base;
    int exponent;

    public Factor(int factor) {
        base = factor;
        exponent = 1;
    }

    public void incExp() {
        exponent++;
    }
}
