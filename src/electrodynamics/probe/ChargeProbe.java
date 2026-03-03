package electrodynamics.probe;

import electrodynamics.Simulation;

public class ChargeProbe extends Probe {
	public int x1;
	public int y1;

	public int x2;
	public int y2;

	public double charge = 0;
	
	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = (x1+x2)/2;
		labelcoord.y = (y1+y2)/2;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		double Q = 0;

		int n_min = 0;
		int n_max = 0;
		int m_min = 0;
		int m_max = 0;

		n_min = Math.min(x1, x2);
		n_max = Math.max(x1, x2);
		m_min = Math.min(y1, y2);
		m_max = Math.max(y1, y2);

		for (int n = n_min; n <= n_max; n++) {
			for (int m = m_min; m <= m_max; m++) {
				Q += e.rho_free[n][m]*(e.ds*e.ds);
			}
		}

		charge = Q*e.depth;
		if (savedatapoint) data.addData(charge, e.time);
	}
	
	@Override
	public ChargeProbe clone() {
		ChargeProbe p = null;
		try {
			p = (ChargeProbe) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}
}