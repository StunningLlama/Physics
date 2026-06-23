// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public class SIUnits extends UnitSystem {
	
	public SIUnits() {
		time_SI = 1;
		len_SI = 1;
		mass_SI = 1;
		charge_SI = 1;
		temp_SI = 1;
	}

	@Override
	public String getSymbol(Quantity q) {
		switch (q) {

		case DIMENSIONLESS: return "";
		case TIME: return "s";
		case LENGTH: return "m";
		case MASS: return "g";
		case CHARGE: return "C";
		case TEMPERATURE: return "K";
		case INFORMATION: return "B";

		case FREQUENCY: return "Hz";
		case ENERGY: return "J";
		case FORCE: return "N";
		case VELOCITY: return "m/s";
		
		case ELECTRIC_POTENTIAL: return "V";
		case ELECTRIC_CURRENT: return "A";
		case ELECTRIC_FIELD: return "V/m";
		case ELECTRIC_FLUX_DENSITY: return "C/m^2";
		case MAGNETIC_FIELD_STRENGTH: return "A/m";
		case MAGNETIC_FLUX_DENSITY: return "T";
		case MAGNETIC_FLUX: return "Wb";
		
		case CHARGE_DENSITY: return "C/m^3";
		case CURRENT_DENSITY: return "A/m^2";
		case NUMBER_DENSITY: return "/m^3";
		case RATE_DENSITY: return "/(m^3 s)";
		case ENERGY_DENSITY: return "J/m^3";
		case POWER_DENSITY: return "W/m^3";
		case ENTROPY_DENSITY_RATE: return "J/(K s m^3)";
		case INTENSITY: return "W/m^2";
		case CONDUCTIVITY: return "S/m";
		case RESISTIVITY: return "Ωm";
		case DIFFUSIVITY: return "(m^2/s)";
		case ELECTRIC_MOBILITY: return "(m^2/(V s))";
		case ENTROPY_RATE: return "J/(K s)";
		case POWER: return "W";
		case RATE: return "/s";
		case RATE_DENSITY_SQUARED: return "/(m^6 s)";
		}

		return "(" + q.name + ")";
	}
	

	@Override
	public double getDefaultMagnitude(Quantity q) {
		if (q == Quantity.MASS)
			return 1e3;
		else
			return 1;
	}

	@Override
	public String getName() {
		return "SI";
	}

}
