// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class ScalarPlot extends Plot {

	public XYSeries data;

	public ScalarPlot() {
		super();
        frame.setTitle("Scalar plot");
	}
	
	@Override
	public void createDataSeries() {
        data = fig.plot("-k", 2.0f, "scalar");
	}
	
	@Override
	public void updatePlot(Simulation e) {
        String title = e.controls.scalarview.getOption().name;
        fig.title(title);
        fig.ylabel(e.controls.scalarview.getOption().unit);
        
		if (frame.isVisible() && e.frame%10 == 0) {

			data.setNotify(false);
			
			data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = t*(x2 - x1) + x1;
				double y = t*(y2 - y1) + y1;

				data.add(t, Utils.bilinearinterp_extrap(e.renderer.scalarfield, x, y, e.nx, e.ny));
			}

			data.setNotify(true);
		}
	}
}
