// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import electrodynamics.Renderer;
import electrodynamics.util.Utils;

public class LinePath extends Path {
	public int x1;
	public int x2;
	public int y1;
	public int y2;
	
	@Override
	public double getX(double t) {
		return t*(x2 - x1) + x1;
	}

	@Override
	public double getY(double t) {
		return t*(y2 - y1) + y1;
	}

	@Override
	public double getArclength() {
		return Utils.length(x2-x1, y2-y1);
	}	
	

	@Override
	public void draw(Renderer r) {
		r.setalphaBG(0.5);
		r.setalphaFG(1.0);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelLine((int)x1, (int)y1, (int)x2, (int)y2);

		r.setalphaBG(0.0);
		r.setalphaFG(1.0);
		r.setColorFloat(0.7f, 0.7f, 0.7f);
		r.drawPixelRectangle((int)x1-1, (int)y1-1, 3, 3);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelRectangle((int)x2-1, (int)y2-1, 3, 3);
	}
}
