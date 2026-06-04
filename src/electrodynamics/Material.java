package electrodynamics;

public class Material {
	public MaterialType type;
	public String name; // For custom materials
	public int cust_id;

	public boolean modified;
	public boolean auto_placed;
	
	public int activated;
	public int conducting;
	public int semiconducting;
	
	public double emf;				// EMF strength
	public double emf_direction;	// EMF direction in radians
	
	public double eps_r;			// Permittivity
	public double mu_r;				// Permeability
	public double rho_back;			// Background charge density
	
	public double ni;				// Equilibrium carrier density
	public double W;				// Work function
	public double Eb;				// Bandgap

	public double D_n;				// Electron diffusion
	public double D_p;				// Hole diffusion
	
	public double v_sat_n;			// Velocity at which carrier velocity saturates
	public double v_sat_p;			// Velocity at which carrier velocity saturates
	
	public double k_rad;			// Radiative recombination rate constant
	public double k_SRH_n;			// Shockley-Read-Hall recombination rate
	public double k_SRH_p;
	public double k_aug_n;			// Auger recombination rate
	public double k_aug_p;
	
	public double absorptivity;
	
	public Material() {
		initialize();
	}
	
	public Material(Material mat) {
		copyFrom(mat);
	}
	
	public void initialize() {
		type = MaterialType.VACUUM;
		name = null;
		cust_id = -1;

		modified = false;
		auto_placed = true;
		
		activated = 1;
		
		setDefaultParameters();
	}
	
	public void setDefaultParameters() {
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

		D_n = 0;
		D_p = 0;
		
		v_sat_n = 0;
		v_sat_p = 0;
		
		k_rad = 0;
		k_aug_n = 0;
		k_aug_p = 0;
		k_SRH_n = 0;
		k_SRH_p = 0;

		absorptivity = 0.0;
	}
    
    public void copyFrom(Material mat) {
    	type = mat.type;
    	name = mat.name; 
    	cust_id = mat.cust_id;

    	modified = mat.modified;
    	auto_placed = mat.auto_placed;
    	
    	activated = mat.activated;
    	conducting = mat.conducting;
    	semiconducting = mat.semiconducting;
    	
    	emf = mat.emf;				
    	emf_direction = mat.emf_direction;	
    	
    	eps_r = mat.eps_r;			
    	mu_r = mat.mu_r;				
    	rho_back = mat.rho_back;			
    	
    	ni = mat.ni;				
    	W = mat.W;				
    	Eb = mat.Eb;				

    	D_n = mat.D_n;				
    	D_p = mat.D_p;				
    	
    	v_sat_n = mat.v_sat_n;			
    	v_sat_p = mat.v_sat_p;			
    	
    	k_rad = mat.k_rad;			
    	k_SRH_n = mat.k_SRH_n;			
    	k_SRH_p = mat.k_SRH_p;
    	k_aug_n = mat.k_aug_n;			
    	k_aug_p = mat.k_aug_p;
    	
    	absorptivity = mat.absorptivity;
    }
    
    public boolean isEmpty() {
    	return type == MaterialType.VACUUM && cust_id == -1;
    }
    
    @Override
    public String toString() {
    	if (cust_id == -1)
    		return type.name;
    	else if (name != null)
    		return name + " [c]";
    	else
    		return "?";
    }
}