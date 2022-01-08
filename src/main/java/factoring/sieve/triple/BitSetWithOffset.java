package factoring.sieve.triple;

import java.util.BitSet;

public class BitSetWithOffset {

    public BitSetWithOffset(BitSet ysMinyOffset, int yOffset) {
        this.ysMinyOffset = ysMinyOffset;
        this.yOffset = yOffset;
    }

    public BitSet ysMinyOffset;

    public int yOffset;
}
