package electrodynamics.probe;

import electrodynamics.Simulation;

public class VoltageProbe extends Probe {
	public int x = 0;
	public int y = 0;

	public double potential = 0;
	
	public void calculateDefaultLabelCoords() {
		labelcoord.x = x-2;
		labelcoord.y = y-5;
	}
	
	public void calcVoltage(Simulation e, boolean savedatapoint) {
		potential = e.F[x][y];
		if (this != e.ground && e.ground != null)
			potential -= e.ground.potential;
		
		if (savedatapoint) data.addData(potential, e.time);
	}
}