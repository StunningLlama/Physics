package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

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
	public boolean ishovering(int mx, int my) {
		return Utils.length(x-mx, y-my) < 3;
	}
	
	@Override
	public VoltageProbe clone() {
		VoltageProbe p = null;
		p = (VoltageProbe) super.clone();
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}
}