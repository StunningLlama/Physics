// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import java.util.ArrayList;

import electrodynamics.Renderer;
import electrodynamics.util.Utils;

public class SegmentedPath extends Path {
	public ArrayList<Integer> vx;
	public ArrayList<Integer> vy;
	public double[] arclength;
	public double total_arclength;
	
	public SegmentedPath() {
		vx = new ArrayList<Integer>();
		vy = new ArrayList<Integer>();
	}
	
	@Override
	public double getX(double t) {
		int i = 0;
		while (t > arclength[i]) i++;
		if (i >= arclength.length) return vx.get(arclength.length-1);
		if (i == 0) return vx.get(0);
		
		double f = (t - arclength[i-1])/(arclength[i] - arclength[i-1]);
		return vx.get(i)*f + vx.get(i-1)*(1-f);
	}

	@Override
	public double getY(double t) {
		int i = 0;
		while (t > arclength[i]) i++;
		if (i >= arclength.length) return vy.get(arclength.length-1);
		if (i == 0) return vy.get(0);
		
		double f = (t - arclength[i-1])/(arclength[i] - arclength[i-1]);
		return vy.get(i)*f + vy.get(i-1)*(1-f);
	}

	@Override
	public void draw(Renderer r) {
		assert(vx.size() == vy.size());
		if (vx.size() > 0) {
			r.setalphaBG(0.5);
			r.setalphaFG(1.0);
			r.setColorFloat(1.0f, 1.0f, 1.0f);
			for (int i = 0; i < vx.size() - 1; i++) {
				r.drawPixelLine((int)vx.get(i), (int)vy.get(i), (int)vx.get(i+1), (int)vy.get(i+1));
			}

			r.setalphaBG(0);
			r.setalphaFG(1.0);
			r.setColorFloat(0.7f, 0.7f, 0.7f);
			r.drawPixelRectangle((int)vx.get(0)-1, (int)vy.get(0)-1, 3, 3);
			r.setColorFloat(1.0f, 1.0f, 1.0f);
			r.drawPixelRectangle((int)vx.get(vx.size()-1)-1, (int)vy.get(vy.size()-1)-1, 3, 3);
		}
	}

	public void addNewJoint(int x, int y) {
		if (vx.size() > 0) {
			if (vx.get(vx.size()-1) == x && vy.get(vy.size()-1) == y) {
				completed = true;
				return;
			}
		}
		vx.add(x);
		vy.add(y);
		computeArc();
	}
	
	public void computeArc() {
		arclength = new double[vx.size()];
		double l_total = 0;
		for (int i = 0; i < vx.size()-1; i++) {
			arclength[i] = l_total;
			l_total += Utils.length(vx.get(i+1)-vx.get(i), vy.get(i+1)-vy.get(i));
		}
		
		arclength[vx.size()-1] = l_total;
		total_arclength = l_total;

		for (int i = 0; i < vx.size(); i++) {
			arclength[i] /= l_total;
		}
	}

	@Override
	public double getArclength() {
		return total_arclength;
	}
}
