package electrodynamics.plot;

import java.util.ArrayList;
import java.util.List;
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
	
	public ProbePlot(String title, String yaxis, String nameprefix, double scalefactor, Predicate<Probe> probefilter) {
		this.title = title;
		this.yaxis = yaxis;
		this.nameprefix = nameprefix;
		this.scalefactor = scalefactor;
		this.probefilter = probefilter;
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
			for (Probe p : e.probes) {
				if (probefilter.test(p)) {
					XYSeries dat = fig.plot("-k", 2.0f, nameprefix + e.getProbeName(index));
					dat.setNotify(false);
					dat.clear();

					for (int i = 0; i < p.data.data_size; i++) {
						if (!Double.isNaN(p.data.data[i]))
							dat.add(1e12*(p.data.time[i] - p.data.time[p.data.data_size-1]), scalefactor*p.data.data[i]);
					}

					index++;
				}
			}

			for (Object dat : fig.dataset.getSeries()) {
				((XYSeries)dat).setNotify(true);
			}
		}
	}
	
	@Override
	public void createPlot(Simulation e) {
		super.createPlot(e);
		frame.setLocation(300+(int)(200*Math.random()), 300+(int)(200*Math.random()));
	}
}
