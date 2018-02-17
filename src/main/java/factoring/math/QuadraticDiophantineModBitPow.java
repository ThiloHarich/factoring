package factoring.math;

/**
 * Looks for solutions xArray for xArray^2 - n = y^2 mod m
 *
 * we calculate f(xArray) = xArray^2 mod m and store the values 2^f(xArray) in an long value squaresMask
 * we calculate g(xArray) = xArray^2-n mod m and test if there is a f(xArray)
 * by testing 2^g(x) & squaresMask != 0
 * we
 **
 * @author thiloharich
 *
 */
public class QuadraticDiophantineModBitPow extends QuadraticDiophantineModBit {


}
