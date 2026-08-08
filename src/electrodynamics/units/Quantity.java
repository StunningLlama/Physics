// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public enum Quantity {
											//  s, m, kg, C, K
	DIMENSIONLESS("", "N", 						0, 0, 0, 0, 0),
	TIME("Time", "t",							1, 0, 0, 0, 0),
	LENGTH("Length", "l",						0, 1, 0, 0, 0),
	MASS("Mass", "m",							0, 0, 1, 0, 0),
	CHARGE("Charge", "Q",						0, 0, 0, 1, 0),
	TEMPERATURE("Temperature", "T",				0, 0, 0, 0, 1),
	INFORMATION("Information", "I",				0, 0, 0, 0, 0),
	
	FREQUENCY("Frequency", "f",					-1, 0, 0, 0, 0),
	RATE("Rate", "R",							-1, 0, 0, 0, 0),
	ENERGY("Energy", "E",						-2, 2, 1, 0, 0),
	POWER("Power", "P",							-3, 2, 1, 0, 0),
	FORCE("Force", "F",							-2, 1, 1, 0, 0),
	VELOCITY("Velocity", "v",					-1, 1, 0, 0, 0),
	DIFFUSIVITY("Diffusivity", "D",				-1, 2, 0, 0, 0),
	ENTROPY_RATE("Entropy rate", "s",			-3, 2, 1, 0, -1),
	
	ELECTRIC_POTENTIAL("Electric potential", "V",		-2, 2, 1, -1, 0),
	ELECTRIC_CURRENT("Electric current", "I",			-1, 0, 0, 1, 0),
	ELECTRIC_FIELD("Electric field", "E",				-2, 1, 1, -1, 0),
	ELECTRIC_FLUX_DENSITY("Electric flux density", "D",	0, -2, 0, 1, 0),
	ELECTRIC_MOBILITY("Mobility", "μ",					1, 0, -1, 1, 0),
	MAGNETIC_FIELD_STRENGTH("Magnetic field", "Η",		-1, -1, 0, 1, 0),
	MAGNETIC_FLUX_DENSITY("Magnetic flux density", "Β", -1, 0, 1, -1, 0),
	MAGNETIC_FLUX("Magnetic flux", "Φ",					-1, 2, 1, -1, 0),
	CONDUCTIVITY("Conductivity", "σ",					1, -3, -1, 2, 0),
	RESISTIVITY("Resistivity", "ρ",						-1, 3, 1, -2, 0),
	
	CHARGE_DENSITY("Charge density", "ρ",				0, -3, 0, 1, 0),
	CURRENT_DENSITY("Current density", "J",				-1, -2, 0, 1, 0),
	NUMBER_DENSITY("Number density", "n",				0, -3, 0, 0, 0, true),
	RATE_DENSITY("Rate density", "r",					-1, -3, 0, 0, 0, true),
	RATE_DENSITY_SQUARED("Rate coefficient", "r",		-1, -6, 0, 0, 0, true),
	ENERGY_DENSITY("Energy density", "ρ", 				-2, -1, 1, 0, 0),
	POWER_DENSITY("Power density", "ρ",					-3, -1, 1, 0, 0),
	ENTROPY_DENSITY_RATE("Entropy density rate", "σ", 	-3, -1, 1, 0, -1),
	INTENSITY("Intensity", "I",							-3, 0, 1, 0, 0);
	
	Quantity(String name, String shorthand, int time, int len, int mass, int charge, int temp) {
		this.name = name;
		this.shorthand = shorthand;
		this.time = time;
		this.len = len;
		this.mass = mass;
		this.charge = charge;
		this.temp = temp;
		this.no_prefixes = false;
	}
	
	Quantity(String name, String shorthand, int time, int len, int mass, int charge, int temp, boolean no_prefixes) {
		this.name = name;
		this.shorthand = shorthand;
		this.time = time;
		this.len = len;
		this.mass = mass;
		this.charge = charge;
		this.temp = temp;
		this.no_prefixes = no_prefixes;
	}
	
	public Quantity multiplyVolume() {
		switch (this) {
		case CHARGE_DENSITY: return CHARGE;
		case NUMBER_DENSITY: return DIMENSIONLESS;
		case RATE_DENSITY: return RATE;
		case ENERGY_DENSITY: return ENERGY;
		case POWER_DENSITY: return POWER;
		case ENTROPY_DENSITY_RATE: return ENTROPY_RATE;
		default: return null;
		}
	}
	
	public Quantity multiplyArea() {
		switch (this) {
		case CURRENT_DENSITY: return ELECTRIC_CURRENT;
		case INTENSITY: return POWER;
		case MAGNETIC_FLUX_DENSITY: return MAGNETIC_FLUX;
		case ELECTRIC_FLUX_DENSITY: return CHARGE;
		default: return null;
		}
	}
	

	public Quantity divideArea() {
		switch (this) {
		case ELECTRIC_CURRENT: return CURRENT_DENSITY;
		case POWER: return INTENSITY;
		case MAGNETIC_FLUX: return MAGNETIC_FLUX_DENSITY;
		case CHARGE: return ELECTRIC_FLUX_DENSITY;
		default: return null;
		}
	}
	
	public String name;
	public String shorthand;
	int time;
	int len;
	int mass;
	int charge;
	int temp;
	boolean no_prefixes;
}