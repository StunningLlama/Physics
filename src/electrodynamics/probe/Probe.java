package electrodynamics.probe;

import electrodynamics.Material;
import electrodynamics.Simulation;

public abstract class Probe implements Cloneable {
	public LabelCoord labelcoord = new LabelCoord();
	public ProbeData data = new ProbeData();

	public abstract void calculateDefaultLabelCoords();
	public abstract void measure(Simulation e, boolean savedatapoint);
	
	public class LabelCoord implements Cloneable {
		public int x = -1;
		public int y = -1;
		
	    @Override
	    public LabelCoord clone() {
	        try {
				return (LabelCoord) super.clone();
			} catch (CloneNotSupportedException e) {
				return null;
			}
	    }
	}
	
	public class ProbeData implements Cloneable {
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
		
	    @Override
	    public ProbeData clone() {
	        ProbeData dat = null;
			try {
				dat = (ProbeData) super.clone();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace();
			}
	        dat.data = data.clone();
	        dat.time = time.clone();
	        return dat;
	    }
	}
}
