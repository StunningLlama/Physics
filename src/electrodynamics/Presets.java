package electrodynamics;

public enum Presets {
	DEFAULT("Default"),
	REALISTIC_SILICON("Realistic (Silicon)"),
	REALSTIC_GERMANIUM("Realistic (Germanium)"),
	REALISTIC_GAAS("Realistic (GaAs)");
	
	String name;

	Presets(String name) {
		this.name = name;
	}
	
	public void applyPreset(Simulation e) {
		e.setDefaultParameters();
		switch(this) {
		case REALISTIC_SILICON:
			break;
		case REALSTIC_GERMANIUM:
			break;
		case REALISTIC_GAAS:
			break;
		case DEFAULT:
			break;
		}
	}
	
	@Override
	public String toString() {
		return name;
	}
}
