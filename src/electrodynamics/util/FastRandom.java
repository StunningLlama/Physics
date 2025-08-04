// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

import java.util.Random;

// Very fast cached random number generator

public class FastRandom {
	private int cache_size = 87961; //Big prime
	private int refresh_ratio = 100; //(# calls to next)/(# calls to rand.nextDouble)
	private double[] cache = new double[cache_size];
	private int ptr = 0;
	private int stride = 1;
	private int counter = 0;
	private Random rand = new Random();

	public FastRandom() {
		rand.setSeed(System.currentTimeMillis());
		for (int i = 0; i < cache_size; i++) {
			cache[i] = rand.nextDouble();
		}
	}

	public double next() {
		counter++;
		if (counter >= refresh_ratio) {
			counter = 0;
			cache[rand.nextInt(cache_size)] = rand.nextDouble();
		}

		ptr = (ptr+stride)%cache_size;
		if (ptr == 0) {
			stride = rand.nextInt(cache_size-1)+1;
		}

		return cache[ptr];
	}
}
