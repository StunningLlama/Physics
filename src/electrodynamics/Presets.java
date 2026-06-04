package electrodynamics;

public enum Presets {
	DEFAULT, REALISTIC_SILICON, REALSTIC_GERMANIUM, REALISTIC_GAAS;
	
	public void applyPreset(Simulation e) {
		e.setDefaultParameters();
	}
}
