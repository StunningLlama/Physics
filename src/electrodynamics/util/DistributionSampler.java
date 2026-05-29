// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

// Generates random samples from a spatial distribution

public class DistributionSampler {

	// Test function
	/*public static void main(String[] args) {
		SpatialDistribution test = new SpatialDistribution();
		test.init(5, 5);
		test.prepare(new double[][] {{0, 0, 0, 0, 0},{0, 0, 0, 0, 0}, {0, 1, 4, 1, 0}, {0,0,0,0,0},{1,1,0,0,0}});
		test.generateNsamples(5000, (c) -> {System.out.println(c.x + "\t" + c.y);});
	}*/


	private double[][] rho;				// Physical density function, unnormalized
	private double rho_tot = 0;			// Total amount
	private double[][] p_dual;			// Probability distribution on dual grid
	private double[][] p_dual_cdf;		// Cumulative probability on dual grid

	private int nx;
	private int ny;

	private int mx;
	private int my;

	FastRandom rand = new FastRandom();

	public void init(int nx, int ny) {
		this.nx = nx;
		this.ny = ny;
		mx = nx-1;
		my = ny-1;

		rho = new double[nx][ny];

		p_dual = new double[nx-1][ny-1];
		p_dual_cdf = new double[nx-1][ny-1];
	}

	public void prepare(double[][] rho_in) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				rho[i][j] = (rho_in[i][j] < 0) ? -rho_in[i][j]: rho_in[i][j];
			}
		}

		rho_tot = 0;
		for (int i = 0; i < mx; i++) {
			for (int j = 0; j < my; j++) {
				p_dual[i][j] = 0.25*(rho[i][j] + rho[i+1][j] + rho[i][j+1] + rho[i+1][j+1]);
				rho_tot += p_dual[i][j];
			}
		}

		for (int i = 0; i < mx; i++) {
			for (int j = 1; j < my; j++) {
				p_dual_cdf[i][j] = p_dual_cdf[i][j-1] + p_dual[i][j]/rho_tot;
			}
			if (i < mx-1) {
				p_dual_cdf[i+1][0] = p_dual_cdf[i][my-1] + p_dual[i+1][0]/rho_tot;
			}
		}
	}

	public double getTotalAmount() {
		return rho_tot;
	}

	private int binarySearch(double p) {
		int M = mx*my;
		int lower = 1;
		int upper = M-1;
		int mid = 0;

		while (lower != upper) {
			mid = (lower+upper)/2;
			if (p <= p_dual_cdf[(mid-1)/my][(mid-1)%my])
			{
				if (upper-lower == 1)
					mid = lower+1;
				upper = mid;
			} else if (p > p_dual_cdf[mid/my][mid%my]){
				lower = mid;
			} else {
				break;
			}
		}

		return mid;
	}

	private Coord generateSample() {
		
		while (true) {
			int index = binarySearch(rand.next());
			int i = index/my;
			int j = index%my;
			double rhomax = Math.max(Math.max(rho[i][j], rho[i+1][j]), Math.max(rho[i][j+1], rho[i+1][j+1]));
			double di = rand.next();
			double dj = rand.next();
			if (rand.next() < bilinearinterp(rho, i+di, j+dj)/rhomax) {
				return new Coord(i+di, j+dj);
			}
		}
	}

	public void generateSamples(double N, GenFunc f) {
		if (!Double.isFinite(rho_tot) || rho_tot <= 0 || N > 100000)
			return;

		int N_floor = (int)Math.floor(N);
		double N_frac = N-N_floor;

		for (int n = 0; n < N_floor; n++) {
			f.generate(generateSample());
		}

		if (rand.next() < N_frac) {
			f.generate(generateSample());
		}
	}

	public double bilinearinterp(double[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;

		if (xfloor < 0) {
			xfloor = 0;
			fx = 0.0;
		} else if (xfloor >= nx - 1) {
			xfloor = nx - 2;
			fx = 1.0;
		}
		if (yfloor < 0) {
			yfloor = 0;
			fy = 0.0;
		} else if (yfloor >= ny - 1) {
			yfloor = ny - 2;
			fy = 1.0;
		}
		double va = array[xfloor][yfloor]*(1.0-fx) + array[xfloor+1][yfloor]*fx;
		double vb = array[xfloor][yfloor+1]*(1.0-fx) + array[xfloor+1][yfloor+1]*fx;

		return va*(1.0-fy) + vb*fy;
	}

	public class Coord {
		public double x;
		public double y;

		public Coord(double x, double y) {
			this.x = x;
			this.y = y;
		}
	}

	public interface GenFunc {
		void generate(Coord c);
	}
}
