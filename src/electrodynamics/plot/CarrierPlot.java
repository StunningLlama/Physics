// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class CarrierPlot extends Plot {
	public XYSeries rho_n_data;
	public XYSeries rho_p_data;

	public CarrierPlot() {
		super();
		LogarithmicAxis yaxis = new LogarithmicAxis("");
		yaxis.setLog10TickLabelsFlag(true);
        fig.chart.getXYPlot().setRangeAxis(yaxis);
        fig.xlabel("Position");
        fig.ylabel("Density (1/m^3)");
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

			rho_n_data.clear();
			rho_p_data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = t*(x2 - x1) + x1;
				double y = t*(y2 - y1) + y1;

				double rho_n = Utils.bilinearinterp_geometric(e.rho_n, x, y, e.nx, e.ny)/e.e_charge;
				if (rho_n > 0) rho_n_data.add(t, rho_n);

				double rho_p = Utils.bilinearinterp_geometric(e.rho_p, x, y, e.nx, e.ny)/e.e_charge;
				if (rho_p > 0) rho_p_data.add(t, rho_p);
			}
		}
	}
}
