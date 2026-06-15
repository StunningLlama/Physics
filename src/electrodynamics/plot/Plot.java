// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.plot;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.units.Units;

public abstract class Plot implements ActionListener {
	public MatlabChart fig;
	public JFrame frame;
	public Font boldfont = new Font(Font.SANS_SERIF, Font.BOLD, 12);
	public Font regularfont = new Font(Font.SANS_SERIF, Font.PLAIN, 12);

	public Path path;
	public JMenuBar menuBar;
	public JMenu menu;
	public JMenuItem menu_log;
	public JMenuItem menu_reset;

	public LogarithmicAxis logaxis;
	public NumberAxis linaxis;
	
	public double xunitquantity;
	public double yunitquantity;

	public Plot() {
		fig = new MatlabChart();
	}
	
	public void initialize() {
		createDataSeries();
		
		fig.RenderPlot();
		fig.title("");
		fig.xlabel("");
		fig.ylabel("");
		fig.grid("on","on");
		fig.font(boldfont);
		fig.legend("northeast", regularfont);

		ChartPanel chartPanel = new ChartPanel(fig.chart);

		logaxis = new LogarithmicAxis("");
		logaxis.setLog10TickLabelsFlag(true);
		logaxis.setStrictValuesFlag(false);
		logaxis.setAllowNegativesFlag(true);
		logaxis.setLabelFont(boldfont);
		logaxis.setTickLabelFont(boldfont);
		
		linaxis = new NumberAxis("");
		linaxis.setAutoRangeIncludesZero(false);
		linaxis.setLabelFont(boldfont);
		linaxis.setTickLabelFont(boldfont);

		frame = new JFrame("");
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		

		menuBar = new JMenuBar();
		frame.setJMenuBar(menuBar);

		menu = new JMenu("Plot settings");
		menuBar.add(menu);

		menu_log = new JCheckBoxMenuItem("Logarithmic y axis");
		menu.add(menu_log);

		menu_reset = new JMenuItem("Clear data");
		menu.add(menu_reset);
		
		menu_log.addActionListener(this);
		menu_reset.addActionListener(this);
		
		frame.add(chartPanel);
		frame.setSize(600, 400);
		frame.setLocationRelativeTo(null);
		frame.setVisible(false);
	}
	
	public abstract void createDataSeries();
	
	public abstract void updatePlot(Simulation e);
	
	public void createPlot(Simulation e, Path p) {
		this.path = p;
		frame.setVisible(true);
	}

	public void reset() {};
	
	@Override
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == menu_log) {
			if (menu_log.isSelected() && fig.chart.getXYPlot().getRangeAxis() != logaxis) {
				fig.chart.getXYPlot().setRangeAxis(logaxis);
			} else if (!menu_log.isSelected() && fig.chart.getXYPlot().getRangeAxis() != linaxis) {
				fig.chart.getXYPlot().setRangeAxis(linaxis);
			}
		} else if (ev.getSource() == menu_reset) {
			reset();
		}
	}
	
	public void setxunits(Units units, Quantity quantity) {
		xunitquantity = units.sys.toSI(1, quantity);
	    String unitname = units.sys.toString(1, quantity, "%.0f");
	    if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
	    fig.xlabel(quantity.name + " (" + unitname + ")");
	}
	
	public void setyunits(Units units, Quantity quantity) {
		yunitquantity = units.sys.toSI(1, quantity);
	    String unitname = units.sys.toString(1, quantity, "%.0f");
	    if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
	    fig.ylabel(quantity.name + " (" + unitname + ")");
	}
	

	public void setxunitsfixed(Units units, Quantity quantity, double value) {
		xunitquantity = value;
	    String unitname = units.sys.toStringSI(value, quantity, "%.0f");
	    if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
	    fig.xlabel(quantity.name + " (" + unitname + ")");
	}
	
	public void setyunitsfixed(Units units, Quantity quantity, double value) {
		yunitquantity = value;
	    String unitname = units.sys.toStringSI(value, quantity, "%.0f");
	    if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
	    fig.ylabel(quantity.name + " (" + unitname + ")");
	}
}
