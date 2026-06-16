// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.util.HashMap;

public enum Preset {
	DEFAULT("Default"),
	SILICON_300K("Silicon @ 300K, clean"),
	GERMANIUM_300K("Germanium @ 300K, clean"),
	GAAS_300K("GaAs @ 300K, clean"),
	GAN_300K("GaN @ 300K, clean"),
	VERSION_1("SemiSim v1.x");
	
	String name;

	Preset(String name) {
		this.name = name;
	}
	
	public void applyPreset(Simulation e) {
		double eVtoJ = 1.6e-19;
		e.setDefaultParameters();
		
		String semi_name = "";
		
		switch(this) {
		case SILICON_300K:
			e.ds = 2e-8;
			e.mu0 = 1e5*1.257e-6;
			e.T = 300;
			
			e.Eg_semi = 1.12*eVtoJ;
			e.chi_semi = 4.05*eVtoJ;
			e.gc_semi = 3.2e19*1e6;
			e.gv_semi = 1.8e19*1e6;
			e.mu_electron_semi = 1400/1e4;
			e.mu_hole_semi = 450/1e4;
			e.v_sat_n_semi = 1e7*1e-2; //Estimate
			e.v_sat_p_semi = 1e7*1e-2; //Estimate
			e.eps_r_semi = 11.7;
			e.d_crit_n = 1e17*1e6; //Estimate
			e.d_crit_p = 7.5e17*1e6; //Conflicting values
			
			e.k_rad_semi = 1.1e-14 * 1e-6;
			e.k_SRH_n_semi = 1e3; //Estimate
			e.k_SRH_p_semi = 1e3; //Estimate
			e.k_aug_n_semi = 1.1e-30 * 1e-12;
			e.k_aug_p_semi = 0.3e-30 * 1e-12;

			setStandardDoping(e);
			setStandardMetalConstants(e);

			e.max_EMF = 5e6;
			e.default_AC_freq = 1e12;
			e.junction_size = 1;
			
			semi_name = "silicon";
			break;
		case GERMANIUM_300K:
			e.ds = 2e-8;
			e.mu0 = 2e4*1.257e-6;
			e.T = 300;

			e.Eg_semi = 0.661*eVtoJ;
			e.chi_semi = 4.0*eVtoJ;
			e.gc_semi = 1.0e19*1e6;
			e.gv_semi = 5.0e18*1e6;
			e.mu_electron_semi = 3900/1e4;
			e.mu_hole_semi = 1900/1e4;
			e.v_sat_n_semi = 7e6*1e-2; //Estimate
			e.v_sat_p_semi = 7e6*1e-2; //Guess
			e.eps_r_semi = 16.2;
			e.d_crit_n = 1e17*1e6; //Estimate
			e.d_crit_p = 1e17*1e6; //Estimate
			
			e.k_rad_semi = 6.41e-14 * 1e-6;
			e.k_SRH_n_semi = 1e3;
			e.k_SRH_p_semi = 1e3;
			e.k_aug_n_semi = 8e-32 * 1e-12; //(Dieter)
			e.k_aug_p_semi = 2.8e-31 * 1e-12; //(Dieter)

			setStandardDoping(e);
			setStandardMetalConstants(e);

			e.max_EMF = 5e6;
			e.default_AC_freq = 1e12;
			e.junction_size = 1;
			
			semi_name = "germanium";
			break;
		case GAAS_300K:
			e.ds = 2e-8;
			e.mu0 = 0.4e4*1.257e-6;
			e.T = 300;
			
			e.Eg_semi = 1.424*eVtoJ;
			e.chi_semi = 4.07*eVtoJ;
			e.gc_semi = 4.7e17*1e6;
			e.gv_semi = 9.0e18*1e6;
			e.mu_electron_semi = 8500/1e4;
			e.mu_hole_semi = 400/1e4;
			e.v_sat_n_semi = 1e7*1e-2; //Estimate
			e.v_sat_p_semi = 1e7*1e-2; //Estimate
			e.eps_r_semi = 12.9;
			e.d_crit_n = 2e17*1e6; //Estimate
			e.d_crit_p = 3e17*1e6; //Estimate
			
			e.k_rad_semi = 7.2e-10 * 1e-6; //Conflicting values
			e.k_SRH_n_semi = 1/(5e-9);
			e.k_SRH_p_semi = 1/(3e-6);
			e.k_aug_n_semi = 5e-30 * 1e-12; //Conflicting values
			e.k_aug_p_semi = 2e-30 * 1e-12; //Conflicting values

			setStandardDoping(e);
			setStandardMetalConstants(e);

			e.g_metal = 2.5e26*1e6;
			e.g_metal_high = 2*e.g_metal;
			e.g_metal_low = 0.25*e.g_metal;
			
			e.W_metal_default = 4.7*eVtoJ;
			e.W_metal_high = 5.4*eVtoJ;
			e.W_metal_low = 4.1*eVtoJ;

			e.max_EMF = 5e6;
			e.default_AC_freq = 1e12;
			e.junction_size = 2;
			
			semi_name = "GaAs";
			break;
		case GAN_300K:
			e.ds = 2e-8;
			e.mu0 = 0.5e5*1.257e-6;
			e.T = 300;

			e.Eg_semi = 3.39*eVtoJ;
			e.chi_semi = 4.1*eVtoJ;
			e.gc_semi = 1.2e18*1e6;
			e.gv_semi = 4.1e19*1e6;
			e.mu_electron_semi = 1000/1e4;
			e.mu_hole_semi = 200/1e4;
			e.v_sat_n_semi = 2.5e7*1e-2; //Estimate
			e.v_sat_p_semi = 2.5e7*1e-2; //Guess
			e.eps_r_semi = 8.9;
			e.d_crit_n = 7.5e17*1e6; //Estimate
			e.d_crit_p = 2e17*1e6; //Estimate
			
			e.k_rad_semi = 1.1e-8 * 1e-6;
			e.k_SRH_n_semi = 1/(1e-3); //Guess
			e.k_SRH_p_semi = 1/(1e-3); //Guess
			e.k_aug_n_semi = 0.5e-30 * 1e-12; //Guess
			e.k_aug_p_semi = 0.5e-30 * 1e-12; //Guess

			setStandardDoping(e);
			setStandardMetalConstants(e);
			
			e.W_metal_default = 5.5*eVtoJ;
			e.W_metal_high = 7*eVtoJ;
			e.W_metal_low = 4.1*eVtoJ;

			e.max_EMF = 5e6;
			e.default_AC_freq = 1e12;
			e.junction_size = 1;
			
			semi_name = "GaN";
			break;
		case DEFAULT:
			semi_name = "semiconductor";
			break;
		case VERSION_1:
			e.junction_size = 3;
			semi_name = "semiconductor";
			break;
		default:
			break;
		}

		e.modified_names = new HashMap<MaterialType, String>(e.default_names);
		for (MaterialType type : MaterialType.values()) {
			String name = e.modified_names.get(type);
			name = name.replace("semiconductor", semi_name);
			e.modified_names.put(type, name);
		}

		e.calculateDependentConstants();
	}
	
	public void setStandardDoping(Simulation e) {
		e.n_light_doping_concentration = 1e15*1e6;
		e.p_light_doping_concentration = 1e15*1e6;
		e.n_default_doping_concentration = 2e16*1e6;
		e.p_default_doping_concentration = 2e16*1e6;
		e.n_heavy_doping_concentration = 4e17*1e6;
		e.p_heavy_doping_concentration = 4e17*1e6;
	}
	
	public void setStandardMetalConstants(Simulation e) {
		double eVtoJ = 1.6e-19;
		
		e.Eg_metal = 1*eVtoJ;
		e.g_metal = 1.25e26*1e6;
		e.g_metal_high = 2*e.g_metal;
		e.g_metal_low = 0.25*e.g_metal;
		e.W_metal_default = 4.6*eVtoJ;
		e.W_metal_high = 5.1*eVtoJ;
		e.W_metal_low = 4.1*eVtoJ;
		e.k_rad_metal = 1.1e-6 * 1e-6;
		e.mu_electron_metal = 200/1e4;
		e.mu_hole_metal = 200/1e4;
		e.v_sat_n_metal = 1e7*1e-2;
		e.v_sat_p_metal = 1e7*1e-2;
		e.eps_r_metal = 5;
	}
	
	@Override
	public String toString() {
		return name;
	}
}
