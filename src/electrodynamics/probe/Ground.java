package electrodynamics.probe;

public class Ground extends VoltageProbe {
	@Override
	public Ground clone() {
		Ground p = null;
		p = (Ground) super.clone();
		return p;
	}
}