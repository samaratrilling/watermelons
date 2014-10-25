package watermelon.group6;

import java.util.*;

import watermelon.sim.Pair;
import watermelon.sim.Point;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
	static double distowall = 1.0;
	static double distotree = 2.0;
	static double distoseed = 2.0;

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

	public ArrayList<seed> compact(ArrayList<Pair> trees, double width, double length, double s) {
		ArrayList<seed> seeds = new ArrayList<seed>();
		int k = 0;
		for (double i = distowall; i <= width - distowall; i += Math.sqrt(3)) {
			int t = k;
			double first;
			if (k == 0)
				first = distowall;
			else
				first = distowall + distoseed / 2;
			for (double j = first; j <= length - distowall; j += distoseed) {
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

	public boolean valid(seed sd, ArrayList<Pair> trees, ArrayList<seed> seeds, double width, double length) {
		final double eps = 1e-6;
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
		double total = 0.0;
		for (int i = 0; i < seedlist.size(); i++) {
			double score;
			double chance = 0.0;
			double totaldis = 0.0;
			double difdis = 0.0;
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i) {
					totaldis = totaldis + Math.pow(distance(seedlist.get(i), seedlist.get(j)), -2);
				}
			}
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i && ((seedlist.get(i).tetraploid && !seedlist.get(j).tetraploid) || (!seedlist.get(i).tetraploid && seedlist.get(j).tetraploid))) {
					difdis = difdis	+ Math.pow(distance(seedlist.get(i), seedlist.get(j)), -2);
				}
			}
			chance = difdis / totaldis;
			score = chance + (1 - chance) * s;
			total = total + score;
		}
		return total;
	}

	public ArrayList<seed> ring(ArrayList<Pair> trees, double width, double length, double s) {
		final double[] xoffset = {-1, 1, 2, 1, -1, -2};
		final double[] yoffset = {-Math.sqrt(3), -Math.sqrt(3), 0, Math.sqrt(3), Math.sqrt(3), 0};
		final double[] cx = {2, 1, -1, -2, -1, 1};
		final double[] cy = {0, Math.sqrt(3), Math.sqrt(3), 0, -Math.sqrt(3), -Math.sqrt(3)};

		ArrayList<seed> seeds = new ArrayList<seed>();
		int layer = 0;
		int cnt = 1;
		while (cnt > 0) {
			layer += 1;
			cnt = 0;
			for (int i = 0; i < trees.size(); ++i) {
				double xtree = trees.get(i).x;
				double ytree = trees.get(i).y;
				boolean flag = (layer % 2 == 0);
				for (int j = 0; j < 6; ++j) {
					double x = xtree + layer * xoffset[j];
					double y = ytree + layer * yoffset[j];
					for (int k = 0; k < layer; ++k) {
						x += cx[j];
						y += cy[j];
						Random random = new Random();
						seed tmp = new seed(x, y, (random.nextInt(2) == 0));
						if (valid(tmp, trees, seeds, width, length)) {
							seeds.add(tmp);
							cnt += 1;
						}
						flag = !flag;
					}
				}
			}
		}
		System.out.printf("#seeds = %d\n", seeds.size());
		return seeds;
	}

	public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {
		//return checkerboard(treelist, width, length, s);
		//return compact(treelist, width, length, s);
		//return ring(treelist, width, length, s);

		/*
		   if (treelist.size() == 0)
		   return compact(treelist, width, length, s);
		   else
		   return ring(treelist, width, length, s);
		   */
		ArrayList<seed> l1 = compact(treelist, width, length, s);
		ArrayList<seed> l2 = ring(treelist, width, length, s);

		// relabel the seeds
		
		double score1 = calculatescore(l1, s);
		double score2 = calculatescore(l2, s);

		if (score1 > score2)
			return l1;
		else
			return l2;
	}
}
