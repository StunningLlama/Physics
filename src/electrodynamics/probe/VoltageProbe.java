// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class VoltageProbe extends Probe {
	public int x = 0;
	public int y = 0;

	public double potential = 0;
	
	@Override
	public void reset() {
		potential = 0;
	}

	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = x-2;
		labelcoord.y = y-5;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		potential = e.V_avg[x][y];
		if (e.hasGround() && this != e.getGround())
			potential -= e.getGround().potential;
		
		if (savedatapoint) data.addData(potential, e.time);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return Utils.length(x-mx, y-my) < 3;
	}
	
	@Override
	public VoltageProbe clone() {
		VoltageProbe p = null;
		p = (VoltageProbe) super.clone();
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}

	@Override
	public void translate(int dx, int dy) {
		x += dx;
		y += dy;
		labelcoord.translate(dx, dy);
	}

	@Override
	public void rotate90(int i_max, int j_max) {
		int y_tmp = y;
		int x_tmp = x;
		x = j_max-y_tmp;
		y = x_tmp;
		labelcoord.rotate90(i_max, j_max);
	}

	@Override
	public void flip_h(int i_min, int i_max) {
		x = (i_min + i_max) - x;
		labelcoord.flip_h(i_min, i_max);
	}

	@Override
	public void flip_v(int j_min, int j_max) {
		y = (j_min + j_max) - y;
		labelcoord.flip_v(j_min, j_max);
	}
	
	@Override
	public boolean intersects(int xmin, int ymin, int xmax, int ymax) {
		return (x >= xmin && x <= xmax && y >= ymin && y <= ymax);
	}
	
	@Override
	public boolean checkInBounds(Simulation e) {
		return (x >= 0 && x < e.nx && y >= 0 && y < e.ny);
	}
}