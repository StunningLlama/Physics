package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class CurrentProbe extends Probe {
	public int x1 = 0;
	public int y1 = 0;
	public int x2 = 0;
	public int y2 = 0;

	public double current = 0;

	@Override
	public void calculateDefaultLabelCoords() {
		double xa = 0.5*(x1+x2);
		double ya = 0.5*(y1+y2);

		double dx = x2 - x1;
		double dy = y2 - y1;
		double len = Utils.length(dx, dy);
		dx = dx/len;
		dy = dy/len;
		if (Math.abs(dx) > Math.abs(dy))
		{
			dx = -Math.abs(dx);
		} else {
			dy = -2*Math.abs(dy);
		}
		
		labelcoord.x = (int)(xa-3*dy)-2;
		labelcoord.y = (int)(ya+4*dx);
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		int x1_t = x1;
		int y1_t = y1;
		int x2_t = x2;
		int y2_t = y2;
		int dy = y2_t - y1_t;
		int dx = x2_t - x1_t;
		float t = (float) 0.5;
		float J = 0;

		if (Math.abs(dx) > Math.abs(dy)) {
			float m = (float) dy / (float) dx;
			t += y1_t;
			dx = (dx < 0) ? -1 : 1;
			m *= dx;
			while (x1_t != x2_t) {
				int x1_prev = x1_t;
				float t_prev = t;

				x1_t += dx;
				t += m;

				J += accumCurrent(x1_prev, (int)t_prev, x1_t, (int)t, e.Jx_free, e.Jy_free, e.ds);

			}
		} else {
			float m = (float) dx / (float) dy;
			t += x1_t;
			dy = (dy < 0) ? -1 : 1;
			m *= dy;
			while (y1_t != y2_t) {
				int y1_prev = y1_t;
				float t_prev = t;

				y1_t += dy;
				t += m;

				J += accumCurrent((int)t_prev, y1_prev, (int)t, y1_t, e.Jx_free, e.Jy_free, e.ds);
			}
		}

		current = J*e.depth;
		if (savedatapoint) data.addData(current, e.time);
	}

	public double accumCurrent(int x0, int y0, int x1, int y1, double[][] Jx, double[][] Jy, double ds) {
		int dx = x1-x0;
		int dy = y1-y0;
		if (dx == 1 && dy == 0) {
			return -Jy[x0+1][y0]*ds;
		} else if (dx == -1 && dy == 0) {
			return Jy[x0][y0]*ds;
		} else if (dx == 0 && dy == 1) {
			return Jx[x0][y0+1]*ds;
		} else if (dx == 0 && dy == -1) {
			return -Jx[x0][y0]*ds;
		} else if (dx == 1 && dy == 1) {
			return (Jx[x0][y0+1] - Jy[x0+1][y0+1])*ds;
		} else if (dx == 1 && dy == -1) {
			return (-Jx[x0][y0] - Jy[x0+1][y0-1])*ds;
		} else if (dx == -1 && dy == 1) {
			return (+Jx[x0][y0+1] + Jy[x0][y0+1])*ds;
		} else if (dx == -1 && dy == -1) {
			return (-Jx[x0][y0] + Jy[x0][y0-1])*ds;
		}
		return 0;
	}
}