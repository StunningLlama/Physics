// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;

public class BandPlot extends Plot {
	public XYSeries E_n_data;
	public XYSeries E_p_data;
	public XYSeries F_n_data;
	public XYSeries F_p_data;

	public BandPlot() {
		super();
        fig.xlabel("Position");
        fig.ylabel("Energy (eV)");
        frame.setTitle("Band diagram");
	}
	
	@Override
	public void createDataSeries() {
        E_n_data = fig.plot("-b", 2.0f, "E_v");
        E_p_data = fig.plot("-r", 2.0f, "E_c");
        F_n_data = fig.plot(".b", 2.0f, "E_Fv");
        F_p_data = fig.plot(".r", 2.0f, "E_Fc");
	}
	
	@Override
	public void updatePlot(Simulation e) {
		if (frame.isVisible() && e.frame%10 == 0) {

			E_n_data.clear();
			E_p_data.clear();
			F_n_data.clear();
			F_p_data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = t*(x2 - x1) + x1;
				double y = t*(y2 - y1) + y1;

				// Add chemical energy and electrostatic energy to get band energy
				E_n_data.add(t, -(e.renderer.bilinearinterp(e.E0_n, x, y)/e.q_n+e.renderer.bilinearinterp(e.phi, x, y)));
				E_p_data.add(t, -(e.renderer.bilinearinterp(e.E0_p, x, y)/e.q_p+e.renderer.bilinearinterp(e.phi, x, y)));
				F_n_data.add(t, -(e.renderer.bilinearinterp(e.F_n, x, y)/e.q_n+e.renderer.bilinearinterp(e.phi, x, y)));
				F_p_data.add(t, -(e.renderer.bilinearinterp(e.F_p, x, y)/e.q_p+e.renderer.bilinearinterp(e.phi, x, y)));
			}
		}
	}
}
