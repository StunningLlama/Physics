package electrodynamics.plot;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;

public class CarrierPlot extends Plot {
	public XYSeries rho_n_data;
	public XYSeries rho_p_data;

	public CarrierPlot() {
		super();
        fig.chart.getXYPlot().setRangeAxis(new LogarithmicAxis(""));
        fig.xlabel("Position");
        fig.ylabel("log(density)");
        frame.setTitle("Carrier density plot");
	}
	
	@Override
	public void createDataSeries() {
        rho_n_data = fig.plot("-b", 2.0f, "n density");
        rho_p_data = fig.plot("-r", 2.0f, "p density");
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

				double rho_n = -e.renderer.bilinearinterp(e.rho_n, x, y);
				if (rho_n > 0) rho_n_data.add(t, rho_n);

				double rho_p = e.renderer.bilinearinterp(e.rho_p, x, y);
				if (rho_p > 0) rho_p_data.add(t, rho_p);
			}
		}
	}
}
