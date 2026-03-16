package electrodynamics.probe;

import electrodynamics.Simulation;

public abstract class Probe implements Cloneable {
	public static int data_size = 100;
	
	public LabelCoord labelcoord = new LabelCoord();
	public ProbeData data = new ProbeData();
	public boolean selected = false;

	public abstract void reset();
	public abstract void calculateDefaultLabelCoords();
	public abstract void measure(Simulation e, boolean savedatapoint);
	public abstract boolean isMouseHovering(int mx, int my);
	public abstract boolean intersects(int xmin, int ymin, int xmax, int ymax);
	public abstract boolean checkInBounds(Simulation e);
	public abstract void translate(int dx, int dy);
	public abstract void rotate90(int i_max, int j_max);
	public abstract void flip_h(int i_min, int i_max);
	public abstract void flip_v(int j_min, int j_max);
	
	@Override
	public Probe clone() {
		try {
			return (Probe) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
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
	    
		public void translate(int dx, int dy) {
			x += dx;
			y += dy;
		}

		public void rotate90(int i_max, int j_max) {
			int y_tmp = y;
			int x_tmp = x;
			x = j_max-y_tmp;
			y = x_tmp;
		}

		public void flip_h(int i_min, int i_max) {
			x = (i_min + i_max) - x;
		}

		public void flip_v(int j_min, int j_max) {
			y = (j_min + j_max) - y;
		}
	}
	
	public class ProbeData implements Cloneable {
		public int data_size = Probe.data_size;
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
