package electrodynamics.plot;

import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;

import electrodynamics.Simulation;

public abstract class Plot {
	public MatlabChart fig;
	public JFrame frame;

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
		fig.font("Helvetica",15);
		fig.legend("northeast");
		fig.legend.setItemFont(new Font("Helvetica", Font.PLAIN, 11));

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
		x1 = e.controls.mx_start_index;
		y1 = e.controls.my_start_index;
		x2 = e.controls.mx_index;
		y2 = e.controls.my_index;
		frame.setVisible(true);
	}
}
