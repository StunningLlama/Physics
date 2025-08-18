package electrodynamics.util;

import java.util.TimerTask;

public abstract class PeriodicTask extends TimerTask {
	double next_start = 0;
	
	public long nextDelay(double frameduration) {
		long end = System.currentTimeMillis();

		next_start = next_start+frameduration;
		long nextDelay = (long)Math.round(next_start-end);
		
		if (nextDelay < 0) {
			next_start = end;
			nextDelay = 0;
		}
		
		return nextDelay;
	}
}
