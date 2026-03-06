package electrodynamics.probe;

import electrodynamics.Simulation;

public class FluxProbe extends Probe {
	public int x1;
	public int y1;

	public int x2;
	public int y2;

	public double flux = 0;
	
	@Override
	public void calculateDefaultLabelCoords() {
		labelcoord.x = (x1+x2)/2;
		labelcoord.y = (y1+y2)/2;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		double phi = 0;

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
				phi += e.Bz[n][m]*(e.ds*e.ds);
			}
		}

		flux = phi;
		if (savedatapoint) data.addData(flux, e.time);
	}
	
	@Override
	public boolean ishovering(int mx, int my) {
		return mx >= Math.min(x1, x2) && mx <= Math.max(x1, x2) && my >= Math.min(y1, y2) && my <= Math.max(y1, y2);
	}
	
	@Override
	public FluxProbe clone() {
		FluxProbe p = null;
		p = (FluxProbe) super.clone();
		p.data = data.clone();
		p.labelcoord = labelcoord.clone();
		return p;
	}
}