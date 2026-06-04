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
		case DEFAULT:
			break;
		case REALISTIC_GAAS:
			break;
		case REALISTIC_SILICON:
			break;
		case REALSTIC_GERMANIUM:
			break;
		}
	}
	
	@Override
	public String toString() {
		return name;
	}
}
