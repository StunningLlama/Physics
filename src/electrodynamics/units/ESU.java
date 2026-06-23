// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public class ESU extends UnitSystem {
	
	public ESU() {
		time_SI = 1;
		len_SI = 0.01;
		mass_SI = 0.001;
		charge_SI = 1/(0.1*29979245800.0);
		temp_SI = 1;
	}

	@Override
	public String getSymbol(Quantity q) {
		switch (q) {
		case DIMENSIONLESS: return "";
		case TIME: return "s";
		case LENGTH: return "m";
		case MASS: return "g";
		case CHARGE: return "statC";
		case TEMPERATURE: return "K";
		case INFORMATION: return "B";
		

		case FREQUENCY: return "Hz";
		case ENERGY: return "erg";
		case FORCE: return "dyn";
		case VELOCITY: return "cm/s";
		
		case ELECTRIC_POTENTIAL: return "statV";
		case ELECTRIC_CURRENT: return "statA";
		case ELECTRIC_FIELD: return "statV/cm";
		case ELECTRIC_FLUX_DENSITY: return "statC/cm^2";
		case MAGNETIC_FIELD_STRENGTH: return "statA/cm";
		case MAGNETIC_FLUX_DENSITY: return "statT";
		case MAGNETIC_FLUX: return "statWb";
		
		case CHARGE_DENSITY: return "statC/cm^3";
		case CURRENT_DENSITY: return "statA/cm^2";
		case NUMBER_DENSITY: return "/cm^3";
		case RATE_DENSITY: return "/(cm^3 s)";
		case ENERGY_DENSITY: return "erg/cm^3";
		case POWER_DENSITY: return "erg/(cm^3 s)";
		case ENTROPY_DENSITY_RATE: return "erg/(K s cm^3)";
		case INTENSITY: return "erg/(cm^2 s)";
		case RATE_DENSITY_SQUARED: return "/(cm^6 s)";

		case CONDUCTIVITY: return "statA/(statV cm)";
		case RESISTIVITY: return "(statV cm)*statA";
		case DIFFUSIVITY: return "(cm^2/s)";
		case ELECTRIC_MOBILITY: return "(cm^2/(statV s))";
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
			return 1;
		else if (q == Quantity.TEMPERATURE)
			return temp_SI;
		else if (q == Quantity.MAGNETIC_FIELD_STRENGTH || q == Quantity.ELECTRIC_FLUX_DENSITY)
			return (4*Math.PI);
		else
			return 1;
	}

	@Override
	public String getName() {
		return "CGS-ESU";
	}

}
