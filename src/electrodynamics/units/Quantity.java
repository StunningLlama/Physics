// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public enum Quantity {
										//  s, m, kg, C, K
	DIMENSIONLESS("", 						0, 0, 0, 0, 0),
	TIME("(Time)", 							1, 0, 0, 0, 0),
	LENGTH("(Length)", 						0, 1, 0, 0, 0),
	MASS("(Mass)", 							0, 0, 1, 0, 0),
	CHARGE("(Charge)", 						0, 0, 0, 1, 0),
	TEMPERATURE("(Temperature)", 				0, 0, 0, 0, 1),
	INFORMATION("(Information)", 				0, 0, 0, 0, 0),
	
	FREQUENCY("(Frequency)", 					-1, 0, 0, 0, 0),
	ENERGY("(Energy)", 							-2, 2, 1, 0, 0),
	FORCE("(Force)", 							-2, 1, 1, 0, 0),
	VELOCITY("(Velocity)", 						-1, 1, 0, 0, 0),
	DIFFUSIVITY("(Diffusivity)", 				-1, 2, 0, 0, 0),
	
	ELECTRIC_POTENTIAL("(Electric potential)",			-2, 2, 1, -1, 0),
	ELECTRIC_CURRENT("(Electric current)", 				-1, 0, 0, 1, 0),
	ELECTRIC_FIELD("(Electric field)", 					-2, 1, 1, -1, 0),
	ELECTRIC_FLUX_DENSITY("(Electric flux density)", 	0, -2, 0, 1, 0),
	ELECTRIC_MOBILITY("(Mobility)", 					1, 0, -1, 1, 0),
	MAGNETIC_FIELD_STRENGTH("(Magnetic field)", 		-1, -1, 0, 1, 0),
	MAGNETIC_FLUX_DENSITY("(Magnetic flux density)", 	-1, 0, 1, -1, 0),
	MAGNETIC_FLUX("(Magnetic flux)", 					-1, 2, 1, -1, 0),
	CONDUCTIVITY("(Conductivity)", 						1, -3, -1, 2, 0),
	RESISTIVITY("(Resistivity)", 						-1, 3, 1, -2, 0),
	
	CHARGE_DENSITY("(Charge density)", 					0, -3, 0, 1, 0),
	CURRENT_DENSITY("(Current density)", 				-1, -2, 0, 1, 0),
	NUMBER_DENSITY("(Number density)", 					0, -3, 0, 0, 0, true),
	RATE_DENSITY("(Rate density)", 						-1, -3, 0, 0, 0, true),
	ENERGY_DENSITY("(Energy density)", 					-2, -1, 1, 0, 0),
	POWER_DENSITY("(Power density)", 					-3, -1, 1, 0, 0),
	ENTROPY_DENSITY_RATE("(Entropy density rate)",		-3, -1, 1, 0, -1),
	INTENSITY("(Intensity)", 							-3, 0, 1, 0, 0);
	
	Quantity(String name, int time, int len, int mass, int charge, int temp) {
		this.name = name;
		this.time = time;
		this.len = len;
		this.mass = mass;
		this.charge = charge;
		this.temp = temp;
		this.no_prefixes = false;
	}
	
	Quantity(String name, int time, int len, int mass, int charge, int temp, boolean no_prefixes) {
		this.name = name;
		this.time = time;
		this.len = len;
		this.mass = mass;
		this.charge = charge;
		this.temp = temp;
		this.no_prefixes = no_prefixes;
	}
	
	public String name;
	int time;
	int len;
	int mass;
	int charge;
	int temp;
	boolean no_prefixes;
}