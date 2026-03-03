package electrodynamics.probe;

import electrodynamics.Simulation;

public class VoltageProbe extends Probe {
	public int x = 0;
	public int y = 0;

	public double potential = 0;

	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = x-2;
		labelcoord.y = y-5;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		potential = e.F[x][y];
		if (this != e.ground && e.ground != null)
			potential -= e.ground.potential;
		
		if (savedatapoint) data.addData(potential, e.time);
	}
	
	@Override
	public VoltageProbe clone() {
		VoltageProbe p = null;
		try {
			p = (VoltageProbe) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}
}