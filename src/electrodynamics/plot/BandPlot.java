// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.util.Utils;

public class BandPlot extends Plot {
	public XYSeries E_n_data;
	public XYSeries E_p_data;
	public XYSeries F_n_data;
	public XYSeries F_p_data;

	@Override
	public void initialize() {
		super.initialize();
        fig.xlabel("Position");
        frame.setTitle("Band diagram");
	}
	
	@Override
	public void createDataSeries() {
        E_n_data = fig.plot("-b", 2.0f, "E_c");
        E_p_data = fig.plot("-r", 2.0f, "E_v");
        F_n_data = fig.plot(".b", 2.0f, "E_Fc");
        F_p_data = fig.plot(".r", 2.0f, "E_Fv");
	}
	
	@Override
	public void updatePlot(Simulation e) {
		if (frame.isVisible() && e.frame%10 == 0) {
	        fig.ylabel("Energy (eV)");
	        setxunitsfixed(e.units, Quantity.LENGTH, 1e-6);
	        
			E_n_data.setNotify(false);
			E_p_data.setNotify(false);
			F_n_data.setNotify(false);
			F_p_data.setNotify(false);
			
			E_n_data.clear();
			E_p_data.clear();
			F_n_data.clear();
			F_p_data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = path.getX(t);
				double y = path.getY(t);
				double len = t*path.getArclength()*e.ds/xunitquantity;

				// Add chemical energy and electrostatic energy to get band energy
				E_n_data.add(len, -(Utils.bilinearinterp_extrap(e.E0_n, x, y, e.nx, e.ny)/e.q_n+Utils.bilinearinterp_extrap(e.phi, e.E0_n, x, y, e.nx, e.ny)));
				E_p_data.add(len, -(Utils.bilinearinterp_extrap(e.E0_p, x, y, e.nx, e.ny)/e.q_p+Utils.bilinearinterp_extrap(e.phi, e.E0_p, x, y, e.nx, e.ny)));
				F_n_data.add(len, -Utils.bilinearinterp_extrap(e.mu_n, x, y, e.nx, e.ny)/e.q_n);
				F_p_data.add(len, -Utils.bilinearinterp_extrap(e.mu_p, x, y, e.nx, e.ny)/e.q_p);
			}

			E_n_data.setNotify(true);
			E_p_data.setNotify(true);
			F_n_data.setNotify(true);
			F_p_data.setNotify(true);
		}
	}
	
	@Override
	public void reset() {
		E_n_data.clear();
		E_p_data.clear();
		F_n_data.clear();
		F_p_data.clear();
	}
}
