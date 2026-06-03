package electrodynamics;

public class Material implements Cloneable {
	MaterialType type = MaterialType.VACUUM;

	boolean modified = false;
	boolean auto_placed = true;
	int activated = 1;
	int conducting = 0;
	int semiconducting = 0;
	double emf = 0.0;			// EMF strength
	double emf_direction = 0.0;	// EMF direction
	double eps_r = 1.0;			// Permittivity
	double mu_r = 1.0;			// Permeability
	double rho_back = 0.0;		// Background charge density
	double ni = 0;				// Equilibrium carrier density
	double W = 0;				// Work function
	double Eb = 0;				// Bandgap
	double r_rel = 1;			// Relative recombination rate
	double absorptivity = 0.0;

	public void erase() {
		type = MaterialType.VACUUM;
		modified = false;
		auto_placed = true;
		activated = 1;
		conducting = 0;
		semiconducting = 0;
		emf = 0.0;
		emf_direction = 0.0;
		eps_r = 1.0;
		mu_r = 1.0;
		rho_back = 0.0;
		ni = 0;
		W = 0;
		Eb = 0;
		r_rel = 1;
		absorptivity = 0;
	}

    @Override
    public Material clone() {
        try {
			return (Material) super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
    }
}