// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.util.Utils;

public class CarrierPlot extends Plot {
	public XYSeries rho_n_data;
	public XYSeries rho_p_data;
	
	public LogarithmicAxis logaxis;
	public NumberAxis linaxis;

	@Override
	public void initialize() {
		super.initialize();
		logaxis = new LogarithmicAxis("");
		logaxis.setLog10TickLabelsFlag(true);
		linaxis = new NumberAxis("");
        fig.chart.getXYPlot().setRangeAxis(logaxis);
        fig.xlabel("Position");
        frame.setTitle("Carrier density plot");
	}
	
	@Override
	public void createDataSeries() {
        rho_n_data = fig.plot("-b", 2.0f, "Electrons");
        rho_p_data = fig.plot("-r", 2.0f, "Holes");
	}
	
	@Override
	public void updatePlot(Simulation e) {
		if (frame.isVisible() && e.frame%10 == 0) {
			if (e.prefs.chkbox_logscale.isSelected() && fig.chart.getXYPlot().getRangeAxis() != logaxis) {
		        fig.chart.getXYPlot().setRangeAxis(logaxis);
			} else if (!e.prefs.chkbox_logscale.isSelected() && fig.chart.getXYPlot().getRangeAxis() != linaxis) {
		        fig.chart.getXYPlot().setRangeAxis(linaxis);
			}
			
	        double unitquantity = e.units.sys.toSI(1, Quantity.NUMBER_DENSITY);
	        String unitname = e.units.sys.toString(1, Quantity.NUMBER_DENSITY, "%.0f");
	        if (unitname.startsWith("1 ")) unitname = "1".concat(unitname.substring(2));
	        fig.ylabel("Density (" + unitname + ")");
	        
			rho_n_data.setNotify(false);
			rho_p_data.setNotify(false);
			
			rho_n_data.clear();
			rho_p_data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = t*(x2 - x1) + x1;
				double y = t*(y2 - y1) + y1;

				double rho_n = Utils.bilinearinterp_geometric_extrap(e.rho_n, x, y, e.nx, e.ny)/e.e_charge/unitquantity;
				if (rho_n > 0) rho_n_data.add(t, rho_n);
				else rho_n_data.add(t, Double.NaN);

				double rho_p = Utils.bilinearinterp_geometric_extrap(e.rho_p, x, y, e.nx, e.ny)/e.e_charge/unitquantity;
				if (rho_p > 0) rho_p_data.add(t, rho_p);
				else rho_p_data.add(t, Double.NaN);
			}

			rho_n_data.setNotify(true);
			rho_p_data.setNotify(true);
		}
	}
}
