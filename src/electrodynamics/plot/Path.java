package electrodynamics.plot;

import electrodynamics.Renderer;

public abstract class Path {
	public abstract double getX(double t);
	public abstract double getY(double t);
	public abstract void draw(Renderer r);
}
