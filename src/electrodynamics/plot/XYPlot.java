// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import java.awt.event.ActionEvent;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenuItem;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.probe.Probe;

public class XYPlot extends Plot {

	public XYSeries data;
	public Probe x;
	public Probe y;
	
	JCheckBoxMenuItem menu_paused;
	JMenuItem menu_measure;
	

	@Override
	public void initialize() {
		super.initialize();

		menu_paused = new JCheckBoxMenuItem("Pause continuous data collection");
		menu.add(menu_paused);
		menu_paused.addActionListener(this);

		menu_measure = new JMenuItem("Measure data point");
		menu.add(menu_measure);
		menu_measure.addActionListener(this);
		
        frame.setTitle("XY plot");
	}
	
	@Override
	public void createDataSeries() {
        data = new XYSeries("x/y", false);
        fig.dataset.addSeries(data);
        fig.FindColor("-k", 2.0f);
	}
	
	@SuppressWarnings("deprecation")
	@Override
	public void updatePlot(Simulation e) {
        
		if (frame.isVisible() && e.frame%10 == 0) {
			if (x != null && y != null) {
				this.setxunits(e.units, x.quantity);
				this.setyunits(e.units, y.quantity);
				
				data.setKey(x.toString() + " vs " + y.toString());
		        
		        if (!menu_paused.isSelected())
		        	data.add(x.value/xunitquantity, y.value/yunitquantity);
			}
		}
	}

	@Override
	public void reset() {
		data.clear();
	}
	

	@Override
	public void actionPerformed(ActionEvent ev) {
		super.actionPerformed(ev);
		if (ev.getSource() == menu_measure) {
			data.add(x.value/xunitquantity, y.value/yunitquantity);
		}
	}
}
