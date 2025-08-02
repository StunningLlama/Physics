// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

// Tool for measuring performance of subroutines
class Timer {
	long tstart = 0;
	String name;
	boolean enabled = true;
	boolean outputavg = false;
	double avgtime = 0;
	double time = 0;
	static boolean allEnabled = false;
	
	public Timer(String name, boolean enabled) {
		this.name = name;
		this.enabled = enabled;
	}

	void start() {
		if (enabled) {
			tstart = System.nanoTime();
		}
	}

	void disableOutput() {
		enabled = false;
	}
	
	void enableOutput() {
		enabled = true;
	}

	void stop() {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.95+time*0.05;
		}
	}

	void stop(String msg) {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.95+time*0.05;
		}
	}
}