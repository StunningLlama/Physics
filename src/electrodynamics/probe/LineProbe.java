// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Renderer;
import electrodynamics.Simulation;
import electrodynamics.Renderer.VectorView;
import electrodynamics.units.Units;
import electrodynamics.util.Utils;

public class LineProbe extends Probe {
	public LineProbe(int mx, int my) {
		super(mx, my);
		x1 = mx;
		x2 = mx;
		y1 = my;
		y2 = my;
		calculateDefaultLabelCoords();
	}

	public int x1 = 0;
	public int y1 = 0;
	public int x2 = 0;
	public int y2 = 0;

	public double[][] vf_x = null;
	public double[][] vf_y = null;
	public VectorView vectorname = VectorView.NONE;
	
	@Override
	public void reset() {
		value = 0;
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
		double[][][] vf = {null, null};
		e.computeVectorField(vf, vectorname);
		vf_x = vf[0];
		vf_y = vf[1];
		
		if (vf_x == null || vf_y == null) return;
		
		int x1_t = x1;
		int y1_t = y1;
		int x2_t = x2;
		int y2_t = y2;
		int dy = y2_t - y1_t;
		int dx = x2_t - x1_t;
		float t = (float) 0.5;
		float J = 0;

		if (Math.abs(dx) > Math.abs(dy)) {
			float m = (float) dy / (float) dx;
			t += y1_t;
			dx = (dx < 0) ? -1 : 1;
			m *= dx;
			while (x1_t != x2_t) {
				int x1_prev = x1_t;
				float t_prev = t;

				x1_t += dx;
				t += m;

				J += accumCurrent(x1_prev, (int)t_prev, x1_t, (int)t, vf_x, vf_y, e.ds);

			}
		} else {
			float m = (float) dx / (float) dy;
			t += x1_t;
			dy = (dy < 0) ? -1 : 1;
			m *= dy;
			while (y1_t != y2_t) {
				int y1_prev = y1_t;
				float t_prev = t;

				y1_t += dy;
				t += m;

				J += accumCurrent((int)t_prev, y1_prev, (int)t, y1_t, vf_x, vf_y, e.ds);
			}
		}

		value = J*e.depth;
		if (savedatapoint) data.addData(value, e.time);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return Utils.length(x1-mx, y1-my) < 3 || Utils.length(x2-mx, y2-my) < 3;
	}
	
	@Override
	public LineProbe clone() {
		LineProbe p = null;
		p = (LineProbe) super.clone();
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}

	public double accumCurrent(int x0, int y0, int x1, int y1, double[][] Jx, double[][] Jy, double ds) {
		int dx = x1-x0;
		int dy = y1-y0;
		if (dx == 1 && dy == 0) {
			return -Jy[x0+1][y0]*ds;
		} else if (dx == -1 && dy == 0) {
			return Jy[x0][y0]*ds;
		} else if (dx == 0 && dy == 1) {
			return Jx[x0][y0+1]*ds;
		} else if (dx == 0 && dy == -1) {
			return -Jx[x0][y0]*ds;
		} else if (dx == 1 && dy == 1) {
			return (Jx[x0][y0+1] - Jy[x0+1][y0+1])*ds;
		} else if (dx == 1 && dy == -1) {
			return (-Jx[x0][y0] - Jy[x0+1][y0-1])*ds;
		} else if (dx == -1 && dy == 1) {
			return (+Jx[x0][y0+1] + Jy[x0][y0+1])*ds;
		} else if (dx == -1 && dy == -1) {
			return (-Jx[x0][y0] + Jy[x0][y0-1])*ds;
		}
		return 0;
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

	@Override
	public void drag(int mx, int my) {
		x2 = mx;
		y2 = my;
		calculateDefaultLabelCoords();
	}

	@Override
	public void draw(Renderer r) {
		r.setalphaBG(1);
		r.setalphaFG(1.0);
		r.setColorFloat(0.5f, 1.0f, 1.0f);
		r.drawPixelRectangle(x1-1, y1-1, 3, 3);
		r.drawPixelRectangle(x2-1, y2-1, 3, 3);

		r.setalphaFG(0.3);
		r.setColorFloat(0.5f, 1.0f, 1.0f);
		r.drawPixelLine(x1, y1, x2, y2);
		double dx_perp = (y2-y1);
		double dy_perp = -(x2-x1);
		double norm = Utils.length(dx_perp, dy_perp);
		if (norm > 0) {
			dx_perp = dx_perp/norm;
			dy_perp = dy_perp/norm;

			int dist = 3;

			r.setPixel((int)Math.round(x1 + dist*dx_perp), (int)Math.round(y1+dist*dy_perp));
			r.setPixel((int)Math.round(x2 + dist*dx_perp), (int)Math.round(y2+dist*dy_perp));
		}
		
		r.setalphaFG(0.1);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelLine((x1 + x2)/2, (y1+y2)/2, labelcoord.x, labelcoord.y);
	}

	@Override
	public String getText(Units units) {
		return shorthand + name + " = " + units.toString_fixedsigfigs(value, quantity);
	}
}