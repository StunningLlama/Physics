// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Renderer;
import electrodynamics.Simulation;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.units.Units;
import electrodynamics.util.Utils;

public class PointProbe extends Probe {
	public PointProbe(int mx, int my) {
		super(mx, my);
		x = mx;
		y = my;
		calculateDefaultLabelCoords();
	}

	public int x = 0;
	public int y = 0;

	public ScalarView scalarname = ScalarView.NONE;
	public double[][] scalarfield = new double[][] {{0}};
	
	@Override
	public void reset() {
		value = 0;
	}

	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = x-2;
		labelcoord.y = y-5;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		e.computeScalarField(scalarfield, x, y, scalarname);
		value = scalarfield[0][0];
		if (savedatapoint) data.addData(value, e.time);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return Utils.length(x-mx, y-my) < 3;
	}
	
	@Override
	public PointProbe clone() {
		PointProbe p = null;
		p = (PointProbe) super.clone();
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

	@Override
	public void drag(int mx, int my) {
		x = mx;
		y = my;
		calculateDefaultLabelCoords();
	}

	@Override
	public void draw(Renderer r) {
		r.setalphaBG(1);
		r.setalphaFG(1.0);
		r.setColorFloat(0.5f, 1.0f, 1.0f);
		r.drawPixelRectangle(x-1, y-1, 3, 3);
		
		r.setalphaFG(0.1);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelLine(x, y, labelcoord.x, labelcoord.y);
	}

	@Override
	public String getText(Units units) {
		return shorthand + name + " = " + units.toString_fixedsigfigs(value, quantity);

	}
}