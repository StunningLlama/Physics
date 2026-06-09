// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

// Tool for measuring performance of subroutines
public class Timer {

	private long tstart = 0;
	private boolean started = false;
	private String name;
	private boolean enabled = true;
	private double avgtime = 0;
	private double time = 0;
	private double coeff = 1;

	public Timer(String name, int smoothing, boolean enabled) {
		this.name = name;
		this.coeff = 1.0/smoothing;
		this.enabled = enabled;
	}

	public void start() {
		if (enabled) {
			tstart = System.nanoTime();
			started = true;
		}
	}

	public void disableOutput() {
		enabled = false;
	}

	public void enableOutput() {
		enabled = true;
	}

	public double getAverageTime() {
		return avgtime;
	}

	public String getName() {
		return name;
	}

	public void stop() {
		if (enabled && started) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*(1-coeff)+time*coeff;
			started = false;
		}
	}

	public void stop(String msg) {
		if (enabled && started) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*(1-coeff)+time*coeff;
			started = false;
		}
	}
}