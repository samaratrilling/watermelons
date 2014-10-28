package watermelon.group6;

import java.util.*;

import watermelon.sim.Pair;
import watermelon.sim.Point;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
	static double distowall = 1.0;
	static double distotree = 2.0;
	static double distoseed = 2.0;
	static final double distBetweenSeeds = Math.sqrt(3)+ 0.0001;

	public void init() {
	}

	static double distance(seed tmp, Pair pair) {
		return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
	}

	static double distance(seed a, seed b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	public ArrayList<seed> checkerboard(ArrayList<Pair> trees, double width, double length, double s) {
		ArrayList<seed> seeds = new ArrayList<seed>();
		int k = 0;
		for (double i = distowall; i <= width - distowall; i = i + distoseed) {
			int t = k;
			for (double j = distowall; j <= length - distowall; j = j + distoseed) {
				seed tmp;
				if (t == 0)
					tmp = new seed(i, j, false);
				else
					tmp = new seed(i, j, true);
				t = 1 - t;
				boolean add = true;
				for (int f = 0; f < trees.size(); f++) {
					if (distance(tmp, trees.get(f)) < distotree) {
						add = false;
						break;
					}
				}
				if (add) {
					seeds.add(tmp);
				}
			}
			k = 1 - k;
		}
		System.out.printf("#seeds = %d\n", seeds.size());
		return seeds;
	}

	public ArrayList<seed> diagonal(ArrayList<Pair> trees, double width,
          double length, double s, boolean packColumns) {
		ArrayList<seed> seeds = new ArrayList<seed>();
		
		// swapped false means squish columns together.
		// swapped true means squish rows together.
		// Used to decide whether or not to flip the last row.
		boolean secondToLastIsTetra = false;
		double lastCoord = 0;
		boolean lastIsTetra = false;
		
		boolean col = false;
		// Alternate coloring every other column (starting on left)
		for (double i = distowall; i <= width - distowall; i += distBetweenSeeds) {
			// Alternate coloring every other row (starting at top)
			boolean row = col;
			// Coord of first seed for this column (top left)
			double first;
			if (col == false)
				first = distowall;
			else
				first = distowall + distoseed / 2;
			for (double j = first; j <= length - distowall; j += distoseed) {
				seed tmp;
				// color it either diploid or tetraploid depending on the row.
				// Save that information
				if (packColumns) {
					// pack columns
					tmp = new seed(j, i, row);
					secondToLastIsTetra = lastIsTetra;
					if (j > lastCoord) {
						lastCoord = j;
					}
					lastIsTetra = col;
				}
				else {
					// pack rows
					tmp = new seed(i, j, row);
					secondToLastIsTetra = lastIsTetra;
					if (j > lastCoord) {
						lastCoord = j;
					}
					lastIsTetra = row;
				}
				
				row = !row;
				boolean add = true;
				for (int f = 0; f < trees.size(); f++) {
					if (distance(tmp, trees.get(f)) < distotree) {
						add = false;
						break;
					}
				}
				if (add) {
					seeds.add(tmp);
				}
			}
			col = !col;
		}
		
		// change the last row to the opposite color
		if (lastIsTetra == secondToLastIsTetra) {
			if (packColumns) {
				// pack columns
				for (seed individ : seeds) {
					if (individ.x == lastCoord) {
						// swap the coloring of the last row
						individ.tetraploid = !individ.tetraploid;
					}
				}
			}
			else {
				//pack rows
				for (seed individ : seeds) {
					if (individ.y == lastCoord) {
						// swap coloring of last column on right
						individ.tetraploid = !individ.tetraploid;
					}
				}
			}
		}
		
		System.out.printf("#seeds = %d\n", seeds.size());
		return seeds;
	}

	public ArrayList<seed> compact (ArrayList<Pair> trees, double width,
          double length, double s) {
      ArrayList<seed> packRows = diagonal(trees, width, length, s, false);
      ArrayList<seed> packColumns = diagonal(trees, length, width, s, true);
      double scoreRows = calculatescore(packRows, s);
      double scoreColumns = calculatescore(packColumns, s);
      if (scoreColumns > scoreRows) {
          return packColumns;
      }
      else {
          return packRows;
      }
  }

	public boolean valid(seed sd, ArrayList<Pair> trees, ArrayList<seed> seeds, double width, double length) {
		final double eps = 0;
		if (sd.x < distowall - eps || width - sd.x < distowall - eps || sd.y < distowall - eps || length - sd.y < distowall - eps)
			return false;
		for (int i = 0; i < trees.size(); ++i) {
			if (distance(sd, trees.get(i)) < distotree - eps)
				return false;
		}
		for (int i = 0; i < seeds.size(); ++i) {
			if (distance(sd, seeds.get(i)) < distoseed - eps)
				return false;
		}
		return true;
	}

	// copied from simulator
	double calculatescore(ArrayList<seed> seedlist, double s) {
		double total = 0;
		
		for (int i = 0; i < seedlist.size(); i++) {
			double score;
			double chance = 0.0;
			double totaldis = 0.0;
			double difdis = 0.0;
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i) {
					totaldis = totaldis+ Math.pow(
									distance(seedlist.get(i),
											seedlist.get(j)), -2);
				}
			}
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i && ((seedlist.get(i).tetraploid && !seedlist.get(j).tetraploid) ||
						(!seedlist.get(i).tetraploid && seedlist.get(j).tetraploid))) {
					difdis = difdis + Math.pow(distance(seedlist.get(i),
							seedlist.get(j)), -2);
				}
			}
			//System.out.println(totaldis);
			//System.out.println(difdis);
			chance = difdis / totaldis;
			score = chance + (1 - chance) * s;
			total = total + score;
		}
		System.out.println("total score: " +total);
		return total;
	}

	public ArrayList<seed> ring(ArrayList<Pair> trees, double width, double length, double s) {
		final double[] xoffset = {-1, 1, 2, 1, -1, -2};
		final double[] yoffset = {-distBetweenSeeds, -distBetweenSeeds, 0, distBetweenSeeds, distBetweenSeeds, 0};
		final double[] cx = {2, 1, -1, -2, -1, 1};
		final double[] cy = {0, distBetweenSeeds, distBetweenSeeds, 0, -distBetweenSeeds, -distBetweenSeeds};

		ArrayList<seed> seeds = new ArrayList<seed>();
		// layer = iteration of the tree placement algorithm
		int layer = 0;
		// count is number of seeds we've been able to place for this particular layer
		int seedsThisLayer = 1;
		// while we can still place seeds on the board
		while (seedsThisLayer > 0) {
			layer += 1;
			seedsThisLayer = 0;
			boolean flagLayer = (layer % 2 == 0);
			for (int i = 0; i < trees.size(); ++i) {
				double xtree = trees.get(i).x;
				double ytree = trees.get(i).y;
				// TODO: within the layer, alternate ploidy between tetra and diploid.
				boolean flag = flagLayer;
				for (int j = 0; j < 6; ++j) {
					double x = xtree + layer * xoffset[j];
					double y = ytree + layer * yoffset[j];
					for (int k = 0; k < layer; ++k) {
						x += cx[j];
						y += cy[j];
						//Random random = new Random();
						//seed tmp = new seed(x, y, (random.nextInt(2) == 0));
						seed tmp = new seed(x, y, flag);
						if (valid(tmp, trees, seeds, width, length)) {
							seeds.add(tmp);
							seedsThisLayer += 1;
						}
						flag = !flag;
					} // end k for
				} // end j for
			} // end i for
		} // end while
		System.out.printf("#seeds = %d\n", seeds.size());
		return seeds;
	}

	public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {

		ArrayList<seed> l1 = compact(treelist, width, length, s);
		ArrayList<seed> l2 = ring(treelist, width, length, s);
		ArrayList<seed> l3 = checkerboard(treelist, width, length, s);
		
		double score1 = calculatescore(l1, s);
		System.out.println("compact: " + score1);
		double score2 = calculatescore(l2, s);
		System.out.println("ring: " + score2);
		double score3 = calculatescore(l3, s);

		double maxScore = Math.max(Math.max(score1, score2), score3);
		if (maxScore == score1) {
			return l1;
		}
		else if (maxScore == score2) {
			return l2;
		}
		else {
			return l3;
		}
	}
}
