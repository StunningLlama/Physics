// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

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
	public double rho_back;			// Background or doping charge density
	
	public double Ec;				// Conduction band energy rel. vacuum
	public double Ev;				// Valence band energy rel. vacuum
	public double gc;				// Effective conduction band density of states @ sim temperature
	public double gv;				// Effective valence band density of states @ sim temperature

	public double D_n;				// Electron diffusivity
	public double D_p;				// Hole diffusivity
	
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
		
		Ec = 0;
		Ev = 0;
		gc = 0;
		gv = 0;

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

		Ec = mat.Ec;
		Ev = mat.Ev;
		gc = mat.gc;
		gv = mat.gv;			

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
    
	public void computeBandstructurePhi(double Phi, double Eg, double ni, double kT)
	{
		double K = ni*ni;

		this.Ec = -Phi + 0.5*Eg;
		this.Ev = -Phi - 0.5*Eg;

		double TS_n = (Math.log(K)*kT + this.Ec - this.Ev)/2.0;
		double TS_p = TS_n;

		this.gc = Math.exp(TS_n/kT);
		this.gv = Math.exp(TS_p/kT);
	}
	
	public void computeBandstructureChi(double chi, double Eg, double ni, double kT)
	{
		double W = chi + 0.5*Eg;
		double K = ni*ni;

		this.Ec = -W + 0.5*Eg;
		this.Ev = -W - 0.5*Eg;

		double TS_n = (Math.log(K)*kT + this.Ec - this.Ev)/2.0;
		double TS_p = TS_n;

		this.gc = Math.exp(TS_n/kT);
		this.gv = Math.exp(TS_p/kT);
	}
	
	public void computeBandstructureChiG(double chi, double Eg, double gc, double gv, double kT)
	{
		this.Ec = -chi;
		this.Ev = -(chi+Eg);

		this.gc = gc;
		this.gv = gv;
	}
	
	public void computeBandstructurePhiG(double Phi, double Eg, double gc, double gv, double kT)
	{
		double Esum = kT*(Math.log(gc) - Math.log(gv)) - 2*Phi;
		this.Ec = 0.5*(Esum + Eg);
		this.Ev = 0.5*(Esum - Eg);

		this.gc = gc;
		this.gv = gv;
	}

	public double calcPhi(double kT)
	{
		double Esum = Ec+Ev;
		return -(0.5*Esum - 0.5*(Math.log(gc) - Math.log(gv))*kT);
	}
	
	public double calc_ni(double kT)
	{
		return Math.sqrt(Math.exp(-(this.Ec-this.Ev)/kT)*this.gc*this.gv);
	}
	
	public double calc_ni2(double kT)
	{
		return Math.exp(-(this.Ec-this.Ev)/kT)*this.gc*this.gv;
	}

	public double calcEquilibriumElectronCharge(double e_charge, double kT) {
		double K = calc_ni2(kT);
		double B = rho_back/e_charge;
		return -e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}

	public double calcEquilibriumHoleCharge(double e_charge, double kT) {
		double K = calc_ni2(kT);
		double B = -rho_back/e_charge;
		return e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}
	
	public double calcConductivity(double e_charge, double kT) {
		double rho_n = calcEquilibriumElectronCharge(e_charge, kT);
		double rho_p = calcEquilibriumHoleCharge(e_charge, kT);
		return e_charge*(-D_n*rho_n/kT + D_p*rho_p/kT);
	}
    
    public boolean isEmpty() {
    	return type == MaterialType.VACUUM && cust_id == -1;
    }
    
    @Override
    public String toString() {
    	if (cust_id == -1)
    		return type.name;
    	else if (name != null)
    		return name;
    	else
    		return "";
    }

    public String getDisplayName() {
    	if (cust_id == -1)
    		return type.name;
    	else if (name != null)
    		return name + "*";
    	else
    		return "?";
    }
}