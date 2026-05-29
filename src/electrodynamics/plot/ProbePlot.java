package electrodynamics.plot;

import java.util.function.Predicate;

import org.jfree.data.xy.XYSeries;

import electrodynamics.Simulation;
import electrodynamics.probe.Probe;
import electrodynamics.units.Quantity;

public class ProbePlot extends Plot {
	
	String yaxis;
	String title;
	String nameprefix;
	double unitquantity;
	double time_quantity = 1e-12;
	Quantity quantity;
	Predicate<Probe> probefilter;
	
	int window_offset = 0;
	
	public ProbePlot(String title, String yaxis, String nameprefix, double scalefactor, Quantity quantity, Predicate<Probe> probefilter, int window_offset) {
		this.title = title;
		this.yaxis = yaxis;
		this.nameprefix = nameprefix;
		this.unitquantity = scalefactor;
		this.quantity = quantity;
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
			int index = 0;
			fig.dataset.removeAllSeries();
			fig.colors.clear();
			fig.strokes.clear();
			
	        String unitname = e.units.sys.toStringSI(unitquantity, quantity, "%.0f");
	        if (unitname.startsWith("1 ")) unitname = unitname.substring(2);
			fig.ylabel(yaxis + " (" + unitname + ")");

	        String unitname_time = e.units.sys.toStringSI(time_quantity, Quantity.TIME, "%.0f");
	        if (unitname_time.startsWith("1 ")) unitname_time = unitname_time.substring(2);
			fig.xlabel("Time (" + unitname_time + ")");
			
			double yrange = 0;
			double tmax = 0;
			for (Probe p : e.probes) {
				if (probefilter.test(p)) {
					XYSeries dat = fig.plot("-k", 2.0f, nameprefix + e.getProbeName(index));
					dat.setNotify(false);
					dat.clear();

					for (int i = 0; i < p.data.data_size; i++) {
						if (!Double.isNaN(p.data.data[i])) {
							double datapointy = p.data.data[i]/unitquantity;
							dat.add(p.data.time[i]/time_quantity, datapointy);
							if (Math.abs(datapointy) > yrange) {
								yrange = Math.abs(datapointy);
							}
						}
					}
					
					if (p.data.time[p.data.data_size-1] > tmax)
						tmax = p.data.time[p.data.data_size-1];

					index++;
				}
			}

			for (Object dat : fig.dataset.getSeries()) {
				((XYSeries)dat).setNotify(true);
			}

			double time_window = Probe.data_size*e.iteration_multiplier*e.controls.plotinterval*e.dt;
			fig.chart.getXYPlot().getDomainAxis().setRange((tmax-time_window)/time_quantity, tmax/time_quantity);
			fig.chart.getXYPlot().getRangeAxis().setRange(-1.1*yrange, 1.1*yrange);
		}
	}
	
	@Override
	public void createPlot(Simulation e) {
		super.createPlot(e);
		frame.setLocation(500, 300+window_offset);
	}
}
