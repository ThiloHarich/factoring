package factoring.hart.sieve;

import java.util.HashSet;
import java.util.Set;

public class Level2 implements Comparable<Level2>{
	static int smoothCountTrunk;

	
	int level;
	int smoothCountInLevel = 0;
	double smoothCountEst = 0;
	int sieveBeginIndex = 0;
	Set<Integer> primesUsed = new HashSet<>();

	
	
	public Level2(int level, int sieveBeginIndex, int smoothCountInLevel, double smoothCountEst) {
		super();
		this.level = level;
		this.sieveBeginIndex = sieveBeginIndex;
		this.smoothCountInLevel = smoothCountInLevel;
		this.smoothCountEst = smoothCountEst;
	}

	public Level2(int level) {
		super();
		this.level = level;
	}
	
	double rating () {
		return (smoothCountEst + 1) * (1.0 + 1.0/(level * (sieveBeginIndex+1)));
	}

	@Override
	public int compareTo(Level2 level) {
		double diff = level.rating() - this.rating();
		if (diff != 0)
			return (int) Math.signum(diff);
		if (level.level != this.level)
			return level.level - this.level;
		return level.sieveBeginIndex - this.sieveBeginIndex;

	}

	@Override
	public String toString() {
		return "Level [level=" + level + ", smoothCountInLevel=" + smoothCountInLevel + ", smoothCountEst="
				+ smoothCountEst + ", sieveBeginIndex=" + sieveBeginIndex + ", rating=" + rating() + "]";
	}



//	@Override
//	public String toString() {
//		return "Level [level=" + level + ", count=" + smoothCountSum + ", sieveBeginIndex=" + sieveBeginIndex + ", rating=" + rating() + "]";
//	} 
	
}
