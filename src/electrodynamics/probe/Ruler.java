package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class Ruler extends Probe {

	public int x1 = 0;
	public int y1 = 0;
	public int x2 = 0;
	public int y2 = 0;

	public double length = 0;
	
	@Override
	public void reset() {
		length = 0;
	}

	@Override
	public void calculateDefaultLabelCoords() {
		double xa = 0.5*(x1+x2);
		double ya = 0.5*(y1+y2);

		double dx = x2 - x1;
		double dy = y2 - y1;
		double len = Utils.length(dx, dy) + Double.MIN_VALUE;
		dx = dx/len;
		dy = dy/len;
		if (Math.abs(dx) > Math.abs(dy))
		{
			dx = -Math.abs(dx);
		} else {
			dy = -2*Math.abs(dy);
		}
		
		labelcoord.x = (int)(xa-3*dy)-2;
		labelcoord.y = (int)(ya+4*dx);
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		length = e.ds*Utils.length(x2-x1, y2-y1);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return Utils.length(x1-mx, y1-my) < 3 || Utils.length(x2-mx, y2-my) < 3;
	}
	
	@Override
	public Ruler clone() {
		Ruler p = null;
		p = (Ruler) super.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}
	
	@Override
	public void translate(int dx, int dy) {
		x1 += dx;
		y1 += dy;
		x2 += dx;
		y2 += dy;
		labelcoord.translate(dx, dy);
	}
	
	@Override
	public void rotate90(int i_max, int j_max) {
		int y_tmp = y1;
		int x_tmp = x1;
		x1 = j_max-y_tmp;
		y1 = x_tmp;
		y_tmp = y2;
		x_tmp = x2;
		x2 = j_max-y_tmp;
		y2 = x_tmp;
		labelcoord.rotate90(i_max, j_max);
	}

	@Override
	public void flip_h(int i_min, int i_max) {
		x1 = (i_min + i_max) - x1;
		x2 = (i_min + i_max) - x2;
		labelcoord.flip_h(i_min, i_max);
	}

	@Override
	public void flip_v(int j_min, int j_max) {
		y1 = (j_min + j_max) - y1;
		y2 = (j_min + j_max) - y2;
		labelcoord.flip_v(j_min, j_max);
	}
	
	@Override
	public boolean intersects(int xmin, int ymin, int xmax, int ymax) {
		return (x1 >= xmin && x1 <= xmax && y1 >= ymin && y1 <= ymax) || (x2 >= xmin && x2 <= xmax && y2 >= ymin && y2 <= ymax);
	}
	
	@Override
	public boolean checkInBounds(Simulation e) {
		return (x1 >= 0 && x1 < e.nx && y1 >= 0 && y1 < e.ny) && (x2 >= 0 && x2 < e.nx && y2 >= 0 && y2 < e.ny);
	}
}
