// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import java.util.function.Predicate;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.probe.Probe;
import electrodynamics.units.Quantity;

public class ProbePlot extends Plot {
	
	public String yaxis;
	public String title;
	public String nameprefix;
	public Quantity quantity;
	public Predicate<Probe> probefilter;
	public double scalefactor;
	public boolean customprobe = false;
	
	public int window_offset = 0;
	
	public ProbePlot(String title, String yaxis, String nameprefix, double scalefactor, Quantity quantity, Predicate<Probe> probefilter, int window_offset) {
		this.title = title;
		this.yaxis = yaxis;
		this.nameprefix = nameprefix;
		this.quantity = quantity;
		this.scalefactor = scalefactor;
		this.probefilter = probefilter;
		this.window_offset = window_offset;
	}

	@Override
	public void initialize() {
		super.initialize();
		frame.setTitle(title);
		fig.chart.getXYPlot().setRangeZeroBaselineVisible(true);
	}
	
	@Override
	public void createDataSeries() {
		fig.plot("-k", 1.0f, "tmp");
	}
	
	public boolean checkProbesExist(Simulation e) {
		for (Probe p : e.probes) {
			if (probefilter.test(p)) {
				return true;
			}
		}
		return false;
	}
	
	@Override
	public void updatePlot(Simulation e) {
		if (frame.isVisible() && e.frame%10 == 0) {
	        setxunitsfixed(e.units, Quantity.TIME, 1e-12);
	        setyunitsfixed(e.units, quantity, scalefactor);
	        
			fig.dataset.removeAllSeries();
			fig.colors.clear();
			fig.strokes.clear();
			
			double yrange = 0;
			double tmax = 0;
			for (Probe p : e.probes) {
				if (probefilter.test(p)) {
					XYSeries dat = fig.plot("-k", 2.0f, nameprefix + p.name);
					dat.setNotify(false);
					dat.clear();

					for (int i = 0; i < p.data.data_size; i++) {
						if (!Double.isNaN(p.data.data[i])) {
							double datapointy = p.data.data[i]/yunitquantity;
							dat.add(p.data.time[i]/xunitquantity, datapointy);
							if (Math.abs(datapointy) > yrange) {
								yrange = Math.abs(datapointy);
							}
						}
					}
					
					if (p.data.time[p.data.data_size-1] > tmax)
						tmax = p.data.time[p.data.data_size-1];
				}
			}
			if (yrange == 0)
				yrange = 1;

			for (Object dat : fig.dataset.getSeries()) {
				((XYSeries)dat).setNotify(true);
			}

			double time_window = Probe.data_size*e.iteration_multiplier*e.controls.plotinterval*e.dt;
			fig.chart.getXYPlot().getDomainAxis().setRange((tmax-time_window)/xunitquantity, tmax/xunitquantity);
			fig.chart.getXYPlot().getRangeAxis().setRange(-1.1*yrange, 1.1*yrange);
		}
	}
	
	@Override
	public void createPlot(Simulation e, Path p) {
		super.createPlot(e, p);
		frame.setLocation(500, 300+window_offset);
	}
}
