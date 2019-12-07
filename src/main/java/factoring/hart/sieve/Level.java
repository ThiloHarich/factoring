package factoring.hart.sieve;

public class Level implements Comparable<Level>{
	int level;
	double count;
	int sieveBeginIndex = 0;
	
	
	public Level(int level, double count, int sieveBeginIndex) {
		super();
		this.level = level;
		this.count = count;
		this.sieveBeginIndex = sieveBeginIndex;
	}


	public Level(int level, double count) {
		super();
		this.level = level;
		this.count = count;
	}
	
	double rating () {
		return count + .5/(Math.cbrt(level) * Math.sqrt(sieveBeginIndex+1));
	}


	@Override
	public int compareTo(Level level) {
		double diff = level.rating() - this.rating();
			return (int) Math.signum(diff);

	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Level other = (Level) obj;
		if (count != other.count)
			return false;
		if (level != other.level)
			return false;
		if (sieveBeginIndex != other.sieveBeginIndex)
			return false;
		return true;
	}


	@Override
	public String toString() {
		return "Level [level=" + level + ", count=" + count + ", sieveBeginIndex=" + sieveBeginIndex + ", rating=" + rating() + "]";
	} 
	
}
