// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

// Fast accurate approximation of log based on:
// https://math.stackexchange.com/questions/5050163/surprisingly-good-approximation-of-fracx-1-ln-x
public class FastLog {

	//int min_exponent = -1022;
	//int orders = 2045;

	private static int min_exponent = -100;
	private static int orders = 200;
	
	private static int length;
	private static long A = (long)0b1111111111 << 52;

	private static double[] x;
	private static double[] log_x;

	static {
		length = orders*4;
		x = new double[length];
		log_x = new double[length];

		for (int i = 0; i < length; i++) {
			x[i] = Math.pow(2, min_exponent+i/4)*(1+(i%4)/4.0);
			log_x[i] = Math.log(x[i]);
		}
	}

	public static double log(double y) {
		int exp = Math.getExponent(y);
		double mantissa = Double.longBitsToDouble(A|Double.doubleToRawLongBits(y)&(~0x7ff0000000000000l));
		int i = (int)(((exp-min_exponent)<<2) + 4*mantissa - 3.5);
		double w = y/x[i];
		return (w-1)/(0.66666666666666667*Math.sqrt(w) + 0.16666666666666667*(w+1)) + log_x[i];
	}
}
