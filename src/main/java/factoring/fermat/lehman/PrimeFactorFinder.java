package factoring.fermat.lehman;

import com.google.common.collect.TreeMultiset;

/**
 * Created by Thilo Harich on 26.03.2018.
 */
public interface PrimeFactorFinder {
    TreeMultiset<Long> factorLong(long n);
}
