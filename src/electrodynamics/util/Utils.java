// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;

public class Utils {
	
	// Compute sech^2(x)
	public static double sech2(double x, double exp2x, boolean useapprox) {
		return useapprox? (1-x*x+0.66666666666666667*x*x*x*x) : 4/(exp2x+2+1/exp2x);
	}

	// Compute tanh(x)
	public static double tanh(double x, double exp2x, boolean useapprox) {
		return useapprox? x*(1-0.33333333333333333*x*x+0.13333333333333333*x*x*x*x) : (exp2x-1)/(exp2x+1);
	}
	
	public static double length(double x, double y) {
		return Math.sqrt(x*x+y*y);
	}
	
	public static double clamp(double val, double min, double max) {
		if ((val != val) || (val < min)) return min;
		if (val > max) return max;
		return val;
	}
	
	public static double max(double x1, double x2, double x3) {
		return Math.max(Math.max(x1, x2), x3);
	}

	public static double max(double x1, double x2, double x3, double x4) {
		return Math.max(Math.max(x1, x2), Math.max(x3, x4));
	}
	
	public static double min(double x1, double x2, double x3) {
		return Math.min(Math.min(x1, x2), x3);
	}

	public static double min(double x1, double x2, double x3, double x4) {
		return Math.min(Math.min(x1, x2), Math.min(x3, x4));
	}
	
	public static double logmean(double x, double y)
	{
		if (x <= 0 || y <= 0)
			return 0;

		//My approximation
		if (Math.abs((x-y)/(x+y)) <  1e-3)
			return (2/3.0)*Math.sqrt(x*y) + (1/6.0)*(x+y);

		return (x-y)/FastLog.log(x/y);
		//return (x-y)/Math.log(x/y);
	}

	// Compute (1/a) tanh(ax)/tanh(x)
	public static double tanhratio(double a, double x, double exp2ax, boolean useapprox) {
		if (useapprox) {
			double a2 = a*a;
			double x2 = x*x;
			return 1 + 0.33333333333333333*(1-a2)*x2 - 0.022222222222222222*(1+5*a2-6*a2*a2)*x2*x2;
		} else {
			double exp2x = FastExp.exp(2*x);
			return (exp2ax-1)*(exp2x+1)/(a*(exp2ax+1)*(exp2x-1));
		}
	}
	
	// Compute x/tanh(x)
	public static double xtanhxm1(double x, double exp2x, boolean useapprox)
	{
		if (useapprox) {
			double x2 = x*x;
			return 1 + 0.33333333333333333*x2-0.022222222222222222*x2*x2;
		} else {
			return x*(exp2x+1)/(exp2x-1);
		}
	}

	public static double bilinearinterp(double[][] array, double x, double y, int nx, int ny) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-6 && Math.abs(y-Math.round(y)) < 1e-6) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= nx) i = nx - 1;
			if (j >= ny) j = ny - 1;
			return array[i][j];
		}

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
	
	public static double bilinearinterp(int[][] array, double x, double y, int nx, int ny) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-6 || Math.abs(y-Math.round(y)) < 1e-6) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= nx) i = nx - 1;
			if (j >= ny) j = ny - 1;
			return array[i][j];
		}

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
	
	public static double bilinearinterp_extrap(double[][] array, double x, double y, int nx, int ny) {
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
		double a = array[xfloor][yfloor];
		double b = array[xfloor+1][yfloor];
		double c = array[xfloor][yfloor+1];
		double d = array[xfloor+1][yfloor+1];

		return a*(1-fx)*(1-fy)+b*fx*(1-fy)+c*(1-fx)*fy+d*fx*fy;
		
		/*double denom = (Double.isFinite(a)? (1-fx)*(1-fy) : 0) + (Double.isFinite(b)? fx*(1-fy) : 0)
			+ (Double.isFinite(c)? (1-fx)*fy : 0) + (Double.isFinite(d)? fx*fy : 0);
		
		double f = (Double.isFinite(a)? a*(1-fx)*(1-fy) : 0) + (Double.isFinite(b)? b*fx*(1-fy) : 0)
			+ (Double.isFinite(c)? c*(1-fx)*fy : 0) + (Double.isFinite(d)? d*fx*fy : 0);
		
		return f/denom;*/
	}
	
	public static double bilinearinterp_extrap(double[][] array, double[][] ref, double x, double y, int nx, int ny) {
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
		double a = array[xfloor][yfloor] + 0*ref[xfloor][yfloor];
		double b = array[xfloor+1][yfloor] + 0*ref[xfloor+1][yfloor];
		double c = array[xfloor][yfloor+1] + 0*ref[xfloor][yfloor+1];
		double d = array[xfloor+1][yfloor+1] + 0*ref[xfloor+1][yfloor+1];

		return a*(1-fx)*(1-fy)+b*fx*(1-fy)+c*(1-fx)*fy+d*fx*fy;
		
		/*double denom = (Double.isFinite(a)? (1-fx)*(1-fy) : 0) + (Double.isFinite(b)? fx*(1-fy) : 0)
			+ (Double.isFinite(c)? (1-fx)*fy : 0) + (Double.isFinite(d)? fx*fy : 0);
		
		double f = (Double.isFinite(a)? a*(1-fx)*(1-fy) : 0) + (Double.isFinite(b)? b*fx*(1-fy) : 0)
			+ (Double.isFinite(c)? c*(1-fx)*fy : 0) + (Double.isFinite(d)? d*fx*fy : 0);
		
		return f/denom;*/
	}
	
	public static double bilinearinterp_geometric_extrap(double[][] array, double x, double y, int nx, int ny) {
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
		
		double a = Math.log(Math.abs(array[xfloor][yfloor]));
		double b = Math.log(Math.abs(array[xfloor+1][yfloor]));
		double c = Math.log(Math.abs(array[xfloor][yfloor+1]));
		double d = Math.log(Math.abs(array[xfloor+1][yfloor+1]));
		
		double denom = (Double.isFinite(a)? (1-fx)*(1-fy) : 0) + (Double.isFinite(b)? fx*(1-fy) : 0)
			+ (Double.isFinite(c)? (1-fx)*fy : 0) + (Double.isFinite(d)? fx*fy : 0);
		
		double f = (Double.isFinite(a)? a*(1-fx)*(1-fy) : 0) + (Double.isFinite(b)? b*fx*(1-fy) : 0)
			+ (Double.isFinite(c)? c*(1-fx)*fy : 0) + (Double.isFinite(d)? d*fx*fy : 0);
		
		return Math.exp(f/denom);
	}
	
	public static double bilinearinterp_length(double[][] Fx, double[][] Fy, double x, double y, int nx, int ny) {
		return length(bilinearinterp(Fx, x-0.5, y, nx, ny), bilinearinterp(Fy, x, y-0.5, nx, ny));
	}
	
	public static <T> List<T> cloneList(List<T> list, UnaryOperator<T> cloner) {
	    List<T> newList = new ArrayList<T>(list.size());
	    for (T element : list) {
	        newList.add(cloner.apply(element));
	    }
	    return newList;
	}
	
	public static <T> List<T> cloneList(List<T> list, UnaryOperator<T> cloner, Predicate<T> criteria) {
	    List<T> newList = new ArrayList<T>(list.size());
	    for (T element : list) {
	    	if (criteria.test(element))
	    		newList.add(cloner.apply(element));
	    }
	    return newList;
	}
	
	private static DecimalFormat df_e = new DecimalFormat("#.########E0");
	private static DecimalFormat df = new DecimalFormat("#.########");
	public static String formatDouble(double d) {
		if ((d != 0 && Math.abs(d) < 1e-3) || Math.abs(d) >= 1e6)  {
			return df_e.format(d);
		} else {
			return df.format(d);
		}
	}
}
