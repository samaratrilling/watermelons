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

	//Pack starting top left
	public ArrayList<seed> topLeftPacking(ArrayList<Pair> trees, double width, double length, double s, boolean packColumns){

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
					//Push the last column to the wall
					if(i+distBetweenSeeds > width - distowall)
						i = width - distowall;
					
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
					
					//Push the last column to the wall
					if(i+distBetweenSeeds > width - distowall)
						i = width - distowall;
					
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
						//Push the last column to the wall
						individ.x = distowall;
						// swap coloring of last column on right
						individ.tetraploid = !individ.tetraploid;
					}
				}
			}
		}
		
		System.out.printf("#seeds = %d\n", seeds.size());
		return seeds;	
	}
	
	//Pack starting top right
	public ArrayList<seed> topRightPacking(ArrayList<Pair> trees, double width, double length, double s, boolean packColumns){

		ArrayList<seed> seeds = new ArrayList<seed>();
		
		// swapped false means squish columns together.
		// swapped true means squish rows together.
		// Used to decide whether or not to flip the last row.
		boolean secondToLastIsTetra = false;
		double lastCoord = 0;
		boolean lastIsTetra = false;
		
		boolean col = false;
		
		// Alternate coloring every other column (starting on right)
		for (double i = width - distowall; i >= distowall; i -= distBetweenSeeds) {
			
			//Push the last row /column to the wall
			// Alternate coloring every other row (starting at top)
			boolean row = col;
			// Coord of first seed for this column (top right)
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
					
					//Push the last column to the wall
					if(i-distBetweenSeeds < distowall)
						i = distowall;
					
					// pack columns
					tmp = new seed(j, i, row);
					secondToLastIsTetra = lastIsTetra;
					if (j > lastCoord) {
						lastCoord = j;
					}
					lastIsTetra = col;					
				}
				else {
					//Push the last column to the wall
					if(i-distBetweenSeeds < distowall)
						i = distowall;
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
				System.out.println("Last Row packing columnss...");
				for (seed individ : seeds) {
					if (individ.x == lastCoord) {
						// swap the coloring of the last row
						individ.tetraploid = !individ.tetraploid;
					}
				}
			}
			else {
				//pack rows
				System.out.println("Last Row packing rows...");
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
	
	
	//Choose the packing which gives highest number of seeds
	public ArrayList<seed> diagonal(ArrayList<Pair> trees, double width,
          double length, double s, boolean packColumns) {
		
		//Pack starting from top left
		ArrayList<seed> seedList1 = topLeftPacking(trees, width, length, s, packColumns);
		ArrayList<seed> seedList2 = topRightPacking(trees, width, length, s, packColumns);
		
		//Calculate the score for each packing
		double topLeftScore     = calculatescore(seedList1, s);
		double topRightScore    = calculatescore(seedList2, s);

		
		System.out.println("Top Left Packing Score: " + topLeftScore);
		System.out.println("Top Right Packing Score: " + topRightScore);				
		
		//Return the highest scoring packing
		if(topLeftScore > topRightScore)	
			return seedList1;
		else 
			return seedList2;	
	}

	public ArrayList<seed> compact (ArrayList<Pair> trees, double width,
          double length, double s) {
	  System.out.println("Calculating scores while Packing Rows....");
      ArrayList<seed> packRows = diagonal(trees, width, length, s, false);
      System.out.println("Calculating scores while Packing Columns....");
      ArrayList<seed> packColumns = diagonal(trees, length, width, s, true);
      double scoreRows = calculatescore(packRows, s);
      double scoreColumns = calculatescore(packColumns, s);
      ArrayList<seed> bestPacking;
      double bestScore;
      if (scoreColumns > scoreRows) {
          bestPacking = packColumns;
          bestScore = scoreColumns;
      }
      else {
          bestPacking = packRows;
          bestScore = scoreRows;
      }
      
      // Start off calling recolor with the improvement equal to the score.
      // The improvement metric will get smaller as subsequent recolorings get closer to
      // the optimal.
      ArrayList<seed> recolored = recolor(bestPacking, bestScore, bestScore, s);
      return recolored;
  }
	
	public ArrayList<seed> recolor(ArrayList<seed> bestPacking, double bestScore,
			double scoreImprovement, double sVal) {
		if (scoreImprovement < .00001) {		
		return bestPacking;
		}
		
		/*// Find the best seed to change first.
		seed bestSeedSoFar = new seed();
		int bestSeedIndex = 0;
		double bestIndividScoreSoFar = bestScore;
		for (seed s : bestPacking) {
			seed flipped = new seed(s.x, s.y, !s.tetraploid);
			ArrayList<seed> tempBoard = new ArrayList<>();
			tempBoard.addAll(bestPacking);
			tempBoard.add(bestPacking.indexOf(s), flipped);
			tempBoard.remove(s);
			double newScore = calculatescore(tempBoard, sVal);
			if (newScore > bestIndividScoreSoFar) {
				bestSeedSoFar = flipped;
				bestSeedIndex = tempBoard.indexOf(flipped);
				bestIndividScoreSoFar = newScore;
			}
		}
		
		// Start off bestPacking with the one best flipped seed
		bestPacking.remove(bestSeedIndex);
		bestPacking.add(bestSeedIndex, bestSeedSoFar);
		*/
		
		// Create new board to test flips on. Only keeps the ones that make it better.
		ArrayList<seed> flippedBoard = new ArrayList<>();
		flippedBoard.addAll(bestPacking);
		double bestNewScore = bestScore;

		// Iterate through bestPacking
		for (int i = 0; i< bestPacking.size(); i++) {
			// Flip each seed.
			seed s = bestPacking.get(i);
			seed flipped = new seed(s.x, s.y, !s.tetraploid);
			flippedBoard.add(bestPacking.indexOf(s), flipped);
			flippedBoard.remove(s);
			
			// If flipping the seed gives a better overall score,
			// keep it and update the score.
			double newScore = calculatescore(flippedBoard, sVal);
			if (newScore > bestNewScore) {
				bestNewScore = newScore;
			}
			// Otherwise switch it back. NOTE: this conditions all seed changes
			// on all previous seed changes.
			else {
				flippedBoard.add(bestPacking.indexOf(s), s);
				flippedBoard.remove(flipped);
			}
		}
		scoreImprovement = bestNewScore - bestScore;
		return recolor(flippedBoard, bestNewScore, scoreImprovement, sVal);
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