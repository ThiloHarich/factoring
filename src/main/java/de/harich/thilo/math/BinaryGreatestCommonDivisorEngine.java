package de.harich.thilo.math;

public class BinaryGreatestCommonDivisorEngine {
    /**
     * Faster binary gcd adapted from OpenJdk's MutableBigInteger.binaryGcd(int, int).
     * @param a
     * @param b
     * @return gcd(a, b)
     */
    public long greatestCommonDivisor(long a, long b) {
        a = Math.abs(a);
        if (b == 0) return a;
        b = Math.abs(b);
        if (a == 0) return b;

        // Right shift a & b till their last bits equal to 1.
        final int aZeros = Long.numberOfTrailingZeros(a);
        final int bZeros = Long.numberOfTrailingZeros(b);
        a >>>= aZeros;
        b >>>= bZeros;

        final int t = (Math.min(aZeros, bZeros));

        while (a != b) {
            if ((a+0x8000000000000000L) > (b+0x8000000000000000L)) { // a > b as unsigned
                a -= b;
                a >>>= Long.numberOfTrailingZeros(a);
            } else {
                b -= a;
                b >>>= Long.numberOfTrailingZeros(b);
            }
        }
        return a<<t;
    }

}
