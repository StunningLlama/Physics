// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.util.Utils;

public class CarrierPlot extends Plot {
	public XYSeries rho_n_data;
	public XYSeries rho_p_data;
	public static boolean debug_charge = false;

	@Override
	public void initialize() {
		super.initialize();
        fig.chart.getXYPlot().setRangeAxis(logaxis);
        menu_log.setSelected(true);
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
	        setxunitsfixed(e.units, Quantity.LENGTH, 1e-6);
	        setyunits(e.units, Quantity.NUMBER_DENSITY);
	        
			rho_n_data.setNotify(false);
			rho_p_data.setNotify(false);
			
			rho_n_data.clear();
			rho_p_data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = path.getX(t);
				double y = path.getY(t);
				double len = t*path.getArclength()*e.ds/xunitquantity;

				if (!debug_charge) {
					double n_n = Utils.bilinearinterp_charge(e.rho_n, x, y, e.nx, e.ny)/e.e_charge/yunitquantity;
					if (n_n > 0) rho_n_data.add(len, n_n);
					else rho_n_data.add(len, Double.NaN);

					double n_p = Utils.bilinearinterp_charge(e.rho_p, x, y, e.nx, e.ny)/e.e_charge/yunitquantity;
					if (n_p > 0) rho_p_data.add(len, n_p);
					else rho_p_data.add(len, Double.NaN);
				} else {
					double n_n = Utils.bilinearinterp(e.rho_n, x, y, e.nx, e.ny)/e.e_charge/yunitquantity;
					rho_n_data.add(len, n_n);
					double n_p = Utils.bilinearinterp(e.rho_p, x, y, e.nx, e.ny)/e.e_charge/yunitquantity;
					rho_p_data.add(len, n_p);
				}
			}

			rho_n_data.setNotify(true);
			rho_p_data.setNotify(true);
		}
	}

	@Override
	public void reset() {
		rho_n_data.clear();
		rho_p_data.clear();
	}
}
