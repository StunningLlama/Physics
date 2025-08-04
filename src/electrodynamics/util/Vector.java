// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

// Two-dimensional vector
public class Vector {
	public double x;
	public double y;

	public Vector(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public Vector copy() {
		return new Vector(x, y);
	}

	public void copy(Vector b) {
		this.x = b.x;
		this.y = b.y;
	}

	public void initialize(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public void add(Vector b) {
		x += b.x;
		y += b.y;
	}

	public void scalarmult(double c) {
		x *= c;
		y *= c;
	}

	public void addmult(Vector b, double c) {
		x += b.x * c;
		y += b.y * c;
	}

	public void rotate(double theta) {
		double xf = x*Math.cos(theta) + y*Math.sin(theta);
		double yf = -x*Math.sin(theta) + y*Math.cos(theta);
		x = xf;
		y = yf;
	}

	public void normalize() {
		double magnitude = Math.sqrt(x*x+y*y);
		if (magnitude != 0) {
			x /= magnitude;
			y /= magnitude;
		}
	}

	public double dot(Vector b) {
		return this.x * b.x + this.y * b.y;
	}
}