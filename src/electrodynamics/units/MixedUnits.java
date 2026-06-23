// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public class MixedUnits extends UnitSystem {
	
	public MixedUnits() {
		time_SI = 1;
		len_SI = 0.01;
		mass_SI = 10000;
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
		case FORCE: return "J/cm";
		case VELOCITY: return "cm/s";
		
		case ELECTRIC_POTENTIAL: return "V";
		case ELECTRIC_CURRENT: return "A";
		case ELECTRIC_FIELD: return "V/cm";
		case ELECTRIC_FLUX_DENSITY: return "C/cm^2";
		case MAGNETIC_FIELD_STRENGTH: return "A/cm";
		case MAGNETIC_FLUX_DENSITY: return "Wb/cm^2";
		case MAGNETIC_FLUX: return "Wb";
		
		case CHARGE_DENSITY: return "C/cm^3";
		case CURRENT_DENSITY: return "A/cm^2";
		case NUMBER_DENSITY: return "/cm^3";
		case RATE_DENSITY: return "/(cm^3 s)";
		case ENERGY_DENSITY: return "J/cm^3";
		case POWER_DENSITY: return "W/cm^3";
		case ENTROPY_DENSITY_RATE: return "J/(K s cm^3)";
		case INTENSITY: return "W/cm^2";
		case CONDUCTIVITY: return "A/(cm V)";
		case RESISTIVITY: return "(cm V)/A";
		case DIFFUSIVITY: return "(cm^2/s)";
		case ELECTRIC_MOBILITY: return "(cm^2/(V s))";

		case ENTROPY_RATE: return "J/(K s)";
		case POWER: return "W";
		case RATE: return "/s";
		case RATE_DENSITY_SQUARED: return "/(cm^6 s)";
		}

		return "(" + q.name + ")";
	}
	

	@Override
	public double getDefaultMagnitude(Quantity q) {
		if (q == Quantity.TIME)
			return time_SI;
		else if (q == Quantity.LENGTH)
			return len_SI;
		else if (q == Quantity.MASS)
			return mass_SI*1e3;
		else if (q == Quantity.CHARGE)
			return charge_SI;
		else if (q == Quantity.TEMPERATURE)
			return temp_SI;
		else
			return 1;
	}

	@Override
	public String getName() {
		return "Engineering";
	}

}
