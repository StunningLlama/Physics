package electrodynamics;

import java.util.HashMap;

public enum Preset {
	DEFAULT("Default"),
	SILICON_300K("Silicon @ 300K"),
	GERMANIUM_300K("Germanium @ 300K"),
	GAAS_300K("GaAs @ 300K");
	
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
			e.default_width = 2*2.56e-6;
			e.mu0 = 5*1.257e-2;
			e.T = 300;
			
			e.E_b_semi = 1.12*eVtoJ;
			e.W_semi = 4.61*eVtoJ;
			e.ni_semi = 9.65e9*1e6;
			e.mu_electron_semi = 1450/1e4;
			e.mu_hole_semi = 500/1e4;
			e.v_sat_n_semi = 1e7*1e-2;
			e.v_sat_p_semi = 1e7*1e-2;
			e.eps_r_semi = 11.7;
			e.a_factor_n = Math.pow(10, -16.9);
			e.a_factor_p = Math.pow(10, -17.9);
			
			e.k_rad_semi = 1.1e-14 * 1e-6;
			e.k_SRH_n_semi = 1e3;
			e.k_SRH_p_semi = 1e3;
			e.k_aug_n_semi = 1.1e-30 * 1e-12;
			e.k_aug_p_semi = 0.3e-30 * 1e-12;

			e.n_light_doping_concentration = 1e15*1e6;
			e.p_light_doping_concentration = 1e15*1e6;
			e.n_default_doping_concentration = 2e16*1e6;
			e.p_default_doping_concentration = 2e16*1e6;
			e.n_heavy_doping_concentration = 4e17*1e6;
			e.p_heavy_doping_concentration = 4e17*1e6;
			
			e.ni_metal = 2.5e17*1e6;
			e.ni_metal_high = 2*e.ni_metal;
			e.ni_metal_low = 0.25*e.ni_metal;
			e.W_metal_default = e.W_semi;
			e.W_metal_high = e.W_semi + 0.3*eVtoJ;
			e.W_metal_low = e.W_semi - 0.3*eVtoJ;
			e.E_b_metal = e.E_b_semi-0.2*eVtoJ;
			e.k_rad_metal = 1.1e-10 * 1e-6;

			e.max_EMF = 5e6;
			e.default_AC_freq = 1e12;
			
			semi_name = "silicon";
			break;
		case GERMANIUM_300K:
			semi_name = "germanium";
			break;
		case GAAS_300K:
			semi_name = "GaAs";
			break;
		case DEFAULT:
			semi_name = "semiconductor";
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
	
	@Override
	public String toString() {
		return name;
	}
}
