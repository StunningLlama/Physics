// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.Font;

import javax.swing.JFrame;

import org.jfree.chart.ChartPanel;
import org.jfree.data.xy.XYSeries;

public class BandPlot {
	XYSeries E_n_data;
	XYSeries E_p_data;
	XYSeries F_n_data;
	XYSeries F_p_data;
	
	MatlabChart fig;
	JFrame frame;
	
	double x1;
	double y1;
	double x2;
	double y2;
	
	public void createPlot() {
        fig = new MatlabChart();
        
        E_n_data = fig.plot("-b", 2.0f, "E_v"); 
        E_p_data = fig.plot("-r", 2.0f, "E_c");
        F_n_data = fig.plot(".b", 2.0f, "E_Fv");
        F_p_data = fig.plot(".r", 2.0f, "E_Fc");
        
        fig.RenderPlot();
        fig.title("");
        fig.xlabel("Position");
        fig.ylabel("Energy (eV)");
        fig.grid("on","on");
        fig.font("Helvetica",15);
        fig.legend("northeast");
        fig.legend.setItemFont(new Font("Helvetica", Font.PLAIN, 11));
        
        ChartPanel chartPanel = new ChartPanel(fig.chart);

        // Create window
        frame = new JFrame("Band diagram");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.add(chartPanel);
        frame.setSize(600, 400);
        frame.setLocationRelativeTo(null); // center
        frame.setVisible(true);
	}
}
