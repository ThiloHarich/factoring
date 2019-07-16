package puzzles;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Points16 {

	public Set<List<Point>> findSolution() {
		final Set<List<Point>> allSolutions = new HashSet<>();
		final Set<Point> hits = new HashSet<>();
		// due to symmetries in x and y we can just take the left upper
		for (int x=0; x <= 5; x++) {
			for (int y=0; y <= 5; y++) {
				final Point point = new Point (x,y);
				final List<Point> polyLine = new ArrayList<>();
				polyLine.add(point);
				final Set<List<Point>> sloutions = findSolution(polyLine, hits, 0);
				allSolutions.addAll(sloutions);
			}
		}
		return allSolutions;
	}

	private Set<List<Point>> findSolution(List<Point> polyLine, Set<Point> hits, int depth) {
		// early aboard
		if (depth >= 6 || hits.size() + (6-depth) * 4 < 16)
			return Collections.emptySet();
		final Set<List<Point>> allSolutions = new HashSet<>();
		final Point point1 = polyLine.get(polyLine.size()-1);
		for (int x=0; x <= 10; x++) {
			for (int y=0; y <= 10; y++) {
				final Point point2 = new Point (x,y);
				final HashSet<Point> newHits = new HashSet<>(hits);
				// we assume each new line must at least add 2 new points
				if (updateHits (point1, point2, newHits) < 2)
					continue;
				final List<Point> solution = new ArrayList<>(polyLine);
				solution.add(point2);
				if (newHits.size() == 16) {
					addSolution(allSolutions, solution);
				}
				allSolutions.addAll(findSolution(solution, newHits, depth + 1));
			}
		}
		return allSolutions;
	}

	public void addSolution(final Set<List<Point>> allSolutions, List<Point> solution) {
		for (final List<Point> existingSol : allSolutions) {
			if (Arrays.equals(solution.toArray(), 1, 5, existingSol.toArray(), 1, 5))
				return;
		}
		System.out.println(solution);
		allSolutions.add(solution);
	}



	private int updateHits(Point point1, Point point2, Set<Point> hits) {
		int hitcount = 0;
		if (point2.x - point1.x == 0) {
			if (point1.x < 3 || point1.x >= 6)
				return 0;
			else {
				final int minY = Math.min(point1.y, point2.y);
				final int maxY = Math.max(point1.y, point2.y);
				for(int y = Math.max(minY, 3); y <= Math.min(maxY,6); y++) {
					final Point hit = new Point(point1.x, y);
					hitcount += hits.add(hit) ? 1 : 0;
				}
				return hitcount;
			}
		}
		final double gradient = (0.0 + point2.y - point1.y) / (0.0 + point2.x - point1.x);
		final int xBegin = Math.max(3,Math.min(point1.x, point2.x));
		final int xEnd = Math.min(6,Math.max(point1.x, point2.x));
		for (int x = xBegin; x <= xEnd; x++) {
			final double y = point1.y + gradient * (x - point1.x);
			if (Math.abs(y- Math.round(y))< .01 && y > 2.9 && y < 6.1) {
				final Point hit = new Point(x, (int) Math.round(y));
				hitcount += hits.add(hit) ? 1 : 0;
			}
		}
		return hitcount;
	}

	public static void main(String[] args) {
		final Points16 points16 = new Points16();
		points16.findSolution().stream().forEach(System.out::println);
	}
}