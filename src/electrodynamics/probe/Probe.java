package electrodynamics.probe;

public class Probe {
	public LabelCoord labelcoord = new LabelCoord();
	public ProbeData data = new ProbeData();
	
	public class LabelCoord {
		public int x = -1;
		public int y = -1;
	}
	
	public class ProbeData {
		public int data_size = 100;
		public double[] data = new double[data_size];
		
		public ProbeData() {
			for (int i = 0; i < data_size; i++) {
				data[i] = Double.NaN;
			}
		}
		
		public void addData(double datapoint) {
			for (int i = 0; i < data_size - 1; i++) {
				data[i] = data[i+1];
			}
			data[data_size-1] = datapoint;
		}
	}
}
