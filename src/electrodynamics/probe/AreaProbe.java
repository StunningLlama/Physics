// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Renderer;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Simulation;
import electrodynamics.units.Units;

public class AreaProbe extends Probe {
	public AreaProbe(int mx, int my) {
		super(mx, my);
		x1 = mx;
		x2 = mx;
		y1 = my;
		y2 = my;
		calculateDefaultLabelCoords();
	}

	public int x1;
	public int y1;

	public int x2;
	public int y2;

	public QuantityType quantitytype = QuantityType.SCALAR;
	public ScalarView scalarname = ScalarView.NONE;
	public double[][] scalarfield = new double[][] {{0}};
	
	@Override
	public void reset() {
		value = 0;
	}
	
	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = (x1+x2)/2;
		labelcoord.y = (y1+y2)/2;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		double Q = 0;

		int n_min = 0;
		int n_max = 0;
		int m_min = 0;
		int m_max = 0;

		n_min = Math.min(x1, x2);
		n_max = Math.max(x1, x2);
		m_min = Math.min(y1, y2);
		m_max = Math.max(y1, y2);
		
		if (scalarfield.length != n_max-n_min+1 || scalarfield[0].length != m_max-m_min+1) {
			scalarfield = new double[n_max-n_min+1][m_max-m_min+1];
		}
		
		e.computeScalarField(scalarfield, n_min, m_min, scalarname);

		for (int n = n_min; n <= n_max; n++) {
			for (int m = m_min; m <= m_max; m++) {
				Q += scalarfield[n-n_min][m-m_min]*(e.ds*e.ds);
			}
		}

		switch (quantitytype) {
		case DENSITY:
			value = Q*e.depth;
			break;
		case FLUX_DENSITY:
			value = Q;
			break;
		case SCALAR:
			value = Q/((n_max-n_min+1)*(m_max-m_min+1)*e.ds*e.ds);
			break;
		}
		if (savedatapoint) data.addData(value, e.time);
	}
	
	@Override
	public boolean isMouseHovering(int mx, int my) {
		return mx >= Math.min(x1, x2) && mx <= Math.max(x1, x2) && my >= Math.min(y1, y2) && my <= Math.max(y1, y2);
	}
	
	@Override
	public AreaProbe clone() {
		AreaProbe p = null;
		p = (AreaProbe) super.clone();
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

		r.setalphaFG(0.3);
		r.setColorFloat(0.5f, 1.0f, 1.0f);
		r.drawPixelLine(x1, y1, x2, y1);
		r.drawPixelLine(x2, y1, x2, y2);
		r.drawPixelLine(x2, y2, x1, y2);
		r.drawPixelLine(x1, y2, x1, y1);
		
		r.setalphaFG(0.1);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelLine((x1 + x2)/2, (y1+y2)/2, labelcoord.x, labelcoord.y);
	}

	@Override
	public String getText(Units units) {
		if (custom) {
			switch (quantitytype) {
			case DENSITY:
				return shorthand + name + " (total) = " + units.toString_fixedsigfigs(value, quantity);
			case FLUX_DENSITY:
				return shorthand + name + " (flux) = " + units.toString_fixedsigfigs(value, quantity);
			case SCALAR:
				return shorthand + name + " (avg) = " + units.toString_fixedsigfigs(value, quantity);
			default:
				return "";
			}
		} else {
			return shorthand + name + " = " + units.toString_fixedsigfigs(value, quantity);
		}
	}
	
	public enum QuantityType {
		DENSITY, FLUX_DENSITY, SCALAR;
	}
}