// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Simulation;

public class FluxProbe extends Probe {
	public int x1;
	public int y1;

	public int x2;
	public int y2;

	public double flux = 0;
	
	@Override
	public void reset() {
		flux = 0;
	}
	
	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = (x1+x2)/2;
		labelcoord.y = (y1+y2)/2;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		double phi = 0;

		int n_min = 0;
		int n_max = 0;
		int m_min = 0;
		int m_max = 0;

		n_min = Math.min(x1, x2);
		n_max = Math.max(x1, x2);
		m_min = Math.min(y1, y2);
		m_max = Math.max(y1, y2);

		for (int n = n_min; n <= n_max; n++) {
			for (int m = m_min; m <= m_max; m++) {
				phi += e.Bz[n][m]*(e.ds*e.ds);
			}
		}

		flux = phi;
		if (savedatapoint) data.addData(flux, e.time);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return mx >= Math.min(x1, x2) && mx <= Math.max(x1, x2) && my >= Math.min(y1, y2) && my <= Math.max(y1, y2);
	}
	
	@Override
	public FluxProbe clone() {
		FluxProbe p = null;
		p = (FluxProbe) super.clone();
		p.data = data.clone();
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