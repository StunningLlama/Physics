// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.util.Utils;

public class ScalarPlot extends Plot {

	public XYSeries data;

	@Override
	public void initialize() {
		super.initialize();
        frame.setTitle("Scalar plot");
        fig.xlabel("Position");
	}
	
	@Override
	public void createDataSeries() {
        data = fig.plot("-k", 2.0f, "scalar");
	}
	
	@Override
	public void updatePlot(Simulation e) {
        String title = e.controls.scalarview.getOption().name;
        fig.title(title);
        Quantity q = e.controls.scalarview.getOption().unit;
        double unitquantity = e.units.sys.toSI(1, q);
        //String unitname = e.units.toString(unitquantity, e.controls.scalarview.getOption().unit);
        String unitname = e.units.sys.toString(1, q, "%.0f");
        if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
        fig.ylabel(q.name + " (" + unitname + ")");
        
		if (frame.isVisible() && e.frame%10 == 0) {

			data.setNotify(false);
			
			data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = path.getX(t);
				double y = path.getY(t);

				data.add(t, Utils.bilinearinterp_extrap(e.renderer.scalarfield, x, y, e.nx, e.ny)/unitquantity);
			}

			data.setNotify(true);
		}
	}
}
