package electrodynamics.probe;

import electrodynamics.Simulation;

public abstract class Probe {
	public LabelCoord labelcoord = new LabelCoord();
	public ProbeData data = new ProbeData();

	public abstract void calculateDefaultLabelCoords();
	public abstract void measure(Simulation e, boolean savedatapoint);
	
	public class LabelCoord {
		public int x = -1;
		public int y = -1;
	}
	
	public class ProbeData {
		public int data_size = 100;
		public double[] data = new double[data_size];
		public double[] time = new double[data_size];
		
		public ProbeData() {
			resetData();
		}
		
		public void resetData() {
			for (int i = 0; i < data_size; i++) {
				data[i] = Double.NaN;
				time[i] = Double.NaN;
			}
		}
		
		public void fixWeirdIssue() {
			if (data == null)
				data = new double[data_size];
			if (time == null)
				time = new double[data_size];
			resetData();
		}
		
		public void addData(double datapoint, double timepoint) {
			for (int i = 0; i < data_size - 1; i++) {
				data[i] = data[i+1];
				time[i] = time[i+1];
			}
			data[data_size-1] = datapoint;
			time[data_size-1] = timepoint;
		}
	}
}
