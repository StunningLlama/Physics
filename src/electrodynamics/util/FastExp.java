// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

public class FastExp {

	private static int min_exponent = -200;
	private static int length = 400;

	private static double[] exp_x;

	static {
		exp_x = new double[length];

		for (int i = 0; i < length; i++) {
			double x = (i + min_exponent)/2.0;
			exp_x[i] = Math.exp(x);
		}
	}

	public static double exp(double y) {
		int i = (int)Math.round(2*y);
		double x = y - i/2.0;
		i = i - min_exponent;
		i = (i < 0? 0 : ((i >= length)? length-1 : i));
		
		double x3 = x*x*x;
		return exp_x[i]*(1+x+0.5*x*x+x3*(0.16666666666666667 + x*0.041666666666666667 + x*x*0.0083333333333333333 + x3*0.0013888888888888889));
	}
}
