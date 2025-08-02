package electrodynamics;

import java.util.Random;

public class FastRandom {
	int cache_size = 87961; //Big prime
	int refresh_ratio = 100; //(# calls to next)/(# calls to rand.nextDouble)
	double[] cache = new double[cache_size];
	int ptr = 0;
	int stride = 1;
	int counter = 0;
	Random rand = new Random();
	
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
