package electrodynamics.plot;

import java.util.function.Predicate;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.probe.Probe;

public class ProbePlot extends Plot {
	
	String yaxis;
	String title;
	String nameprefix;
	double scalefactor;
	Predicate<Probe> probefilter;
	
	int window_offset = 0;
	
	public ProbePlot(String title, String yaxis, String nameprefix, double scalefactor, Predicate<Probe> probefilter, int window_offset) {
		this.title = title;
		this.yaxis = yaxis;
		this.nameprefix = nameprefix;
		this.scalefactor = scalefactor;
		this.probefilter = probefilter;
		this.window_offset = window_offset;
	}

	@Override
	public void initialize() {
		super.initialize();
		fig.xlabel("Time [ps]");
		fig.ylabel(yaxis);
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
			int index = 0;
			fig.dataset.removeAllSeries();
			fig.colors.clear();
			fig.strokes.clear();
			
			double yrange = 0;
			for (Probe p : e.probes) {
				if (probefilter.test(p)) {
					XYSeries dat = fig.plot("-k", 2.0f, nameprefix + e.getProbeName(index));
					dat.setNotify(false);
					dat.clear();

					for (int i = 0; i < p.data.data_size; i++) {
						if (!Double.isNaN(p.data.data[i])) {
							double datapointy = scalefactor*p.data.data[i];
							dat.add(1e12*(p.data.time[i] - p.data.time[p.data.data_size-1]), datapointy);
							if (Math.abs(datapointy) > yrange) {
								yrange = Math.abs(datapointy);
							}
						}
					}

					index++;
				}
			}

			for (Object dat : fig.dataset.getSeries()) {
				((XYSeries)dat).setNotify(true);
			}

			fig.chart.getXYPlot().getDomainAxis().setRange(-1e12*100*e.iteration_multiplier*e.controls.plotinterval*e.dt, 0);
			fig.chart.getXYPlot().getRangeAxis().setRange(-1.1*yrange, 1.1*yrange);
		}
	}
	
	@Override
	public void createPlot(Simulation e) {
		super.createPlot(e);
		frame.setLocation(500, 300+window_offset);
	}
}
