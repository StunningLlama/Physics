// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.probe.Probe;
import electrodynamics.units.Quantity;

public class XYPlot extends Plot {

	public XYSeries data;
	public Probe x;
	public Probe y;

	@Override
	public void initialize() {
		super.initialize();
        frame.setTitle("XY plot");
	}
	
	@Override
	public void createDataSeries() {
        data = new XYSeries("", false);
        fig.dataset.addSeries(data);
        fig.FindColor("-k", 2.0f);
	}
	
	@Override
	public void updatePlot(Simulation e) {
        
		if (frame.isVisible() && e.frame%10 == 0) {
			if (x != null && y != null) {
		        Quantity q = x.quantity;
		        double unitquantity_x = e.units.sys.toSI(1, q);
		        String unitname_x = e.units.sys.toString(1, q, "%.0f");
		        if (unitname_x.startsWith("1 ")) unitname_x = unitname_x.substring(2);
		        fig.xlabel(q.name + " (" + unitname_x + ")");
		        
		        q = y.quantity;
		        double unitquantity_y = e.units.sys.toSI(1, q);
		        String unitname_y = e.units.sys.toString(1, q, "%.0f");
		        if (unitname_y.startsWith("1 ")) unitname_y = unitname_y.substring(2);
		        fig.ylabel(q.name + " (" + unitname_y + ")");
		        
				data.add(x.value/unitquantity_x, y.value/unitquantity_y);
			}
		}
	}

	@Override
	public void reset() {
		data.clear();
	}
}
