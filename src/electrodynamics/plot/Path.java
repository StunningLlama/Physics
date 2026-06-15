// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import electrodynamics.Renderer;

public abstract class Path {
	public abstract double getX(double t);
	public abstract double getY(double t);
	public abstract double getArclength();
	public abstract void draw(Renderer r);
	public boolean completed = false;
}
