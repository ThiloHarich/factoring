package factoring.fermat.lehman;

public class Multiplier implements Comparable<Multiplier>{
	double sqrtInv;
	double sqrt;
	int value;
	public Multiplier(int value, double sqrt, double sqrtInv) {
		this.value = value;
		this.sqrt = sqrt;
		this.sqrtInv = sqrtInv;
	}
	@Override
	public int compareTo(Multiplier mult) {
		return mult.value - value ;
	}
	@Override
	public String toString() {
		return "Multiplier [value=" + value + ", sqrt=" + sqrt + ", sqrtInv=" + sqrtInv + "]";
	}

}
