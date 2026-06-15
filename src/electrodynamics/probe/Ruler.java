package electrodynamics.probe;

import electrodynamics.Renderer;
import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.units.Units;
import electrodynamics.util.Utils;

public class Ruler extends LineProbe {
	public Ruler(int mx, int my) {
		super(mx, my);
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		value = e.ds*Utils.length(x2-x1, y2-y1);
	}

	@Override
	public void draw(Renderer r) {
		r.setalphaBG(1);
		r.setalphaFG(1.0);
		r.setColorFloat(1.0f, 0.8f, 0.5f);
		r.drawPixelRectangle(x1-1, y1-1, 3, 3);
		r.drawPixelRectangle(x2-1, y2-1, 3, 3);

		r.setalphaFG(0.3);
		r.setColorFloat(1.0f, 0.8f, 0.5f);
		r.drawPixelLine(x1, y1, x2, y2);
		
		r.setalphaFG(0.1);
		r.setColorFloat(1.0f, 1.0f, 1.0f);
		r.drawPixelLine((x1 + x2)/2, (y1+y2)/2, labelcoord.x, labelcoord.y);
	}

	@Override
	public String getText(Units units) {
		return units.toString_fixedsigfigs(value, Quantity.LENGTH);
	}
}
