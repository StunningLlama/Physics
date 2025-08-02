// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

// Tool for measuring performance of subroutines
public class Timer {

	private long tstart = 0;
	private String name;
	private boolean enabled = true;
	private double avgtime = 0;
	private double time = 0;

	public static boolean allEnabled = false;

	public Timer(String name, boolean enabled) {
		this.name = name;
		this.enabled = enabled;
	}

	public void start() {
		if (enabled) {
			tstart = System.nanoTime();
		}
	}

	public void disableOutput() {
		enabled = false;
	}

	public void enableOutput() {
		enabled = true;
	}

	public double getTime() {
		return time;
	}

	public double getAverageTime() {
		return avgtime;
	}

	public String getName() {
		return name;
	}

	public void stop() {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.95+time*0.05;
		}
	}

	public void stop(String msg) {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.95+time*0.05;
		}
	}
}