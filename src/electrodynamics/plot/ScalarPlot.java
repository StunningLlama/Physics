package electrodynamics.plot;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Renderer;
import electrodynamics.Simulation;

public class ScalarPlot extends Plot {

	public XYSeries data;

	public ScalarPlot() {
		super();
        frame.setTitle("Scalar plot");
	}
	
	@Override
	public void createDataSeries() {
        data = fig.plot("-", 2.0f, "scalar");
	}
	
	@Override
	public void updatePlot(Simulation e) {
        String title = ((Renderer.ScalarView)e.opts.gui_view.getSelectedItem()).name.split("\\:\\ ")[1];
        fig.title(title);
        fig.ylabel(((Renderer.ScalarView)e.opts.gui_view.getSelectedItem()).unit);
        
		if (frame.isVisible() && e.frame%10 == 0) {

			data.clear();

			for (int n = 0; n <= 100; n++) {
				double t = n/100.0;
				double x = t*(x2 - x1) + x1;
				double y = t*(y2 - y1) + y1;

				data.add(t, e.renderer.bilinearinterp(e.renderer.scalarfield, x, y));
				updatePlot(e);
			}
		}
	}
}
