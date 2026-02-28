// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;

import electrodynamics.Simulation;

public abstract class Plot {
	public MatlabChart fig;
	public JFrame frame;
	public Font boldfont = new Font(Font.SANS_SERIF, Font.BOLD, 12);
	public Font regularfont = new Font(Font.SANS_SERIF, Font.PLAIN, 12);

	public double x1;
	public double y1;
	public double x2;
	public double y2;

	public Plot() {
		fig = new MatlabChart();
		
		createDataSeries();
		
		fig.RenderPlot();
		fig.title("");
		fig.xlabel("");
		fig.ylabel("");
		fig.grid("on","on");
		fig.font(boldfont);
		fig.legend("northeast", regularfont);

		ChartPanel chartPanel = new ChartPanel(fig.chart);

		frame = new JFrame("");
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.add(chartPanel);
		frame.setSize(600, 400);
		frame.setLocationRelativeTo(null);
		frame.setVisible(false);
	}
	
	public abstract void createDataSeries();
	
	public abstract void updatePlot(Simulation e);
	
	public void createPlot(Simulation e) {
		x1 = e.controls.mx_start;
		y1 = e.controls.my_start;
		x2 = e.controls.mx;
		y2 = e.controls.my;
		frame.setVisible(true);
	}
}
