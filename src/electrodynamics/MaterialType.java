// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

public enum MaterialType
{
	EMF					("Voltage source",						230, 216, 46, 230),
	AC_EMF				("AC voltage source",					230, 150, 216, 230),
	CURRENT				("Current source",						100, 216, 216, 230),
	SWITCH				("Switch",								194, 194, 194, 120),
	METAL				("Metal",								153, 153, 153, 130),
	METAL_HIGH_C		("Conductive metal",					191, 191, 191, 130),
	METAL_LOW_C			("Resistive metal",						94, 94, 94, 130),
	METAL_HIGH_W		("High workfunction metal",				163, 116, 116, 130),
	METAL_LOW_W			("Low workfunction metal",				116, 121, 163, 130),
	SEMI				("Intrinsic semiconductor",				207,  161, 212, 110),
	SEMI_P_TYPE			("P-type semiconductor",				191,  74,  34, 110),
	SEMI_N_TYPE			("N-type semiconductor",				 84, 123, 191, 110),
	SEMI_HEAVY_P_TYPE	("Heavily doped P-type semiconductor",	204,  41,  41, 110),
	SEMI_HEAVY_N_TYPE	("Heavily doped N-type semiconductor",	 39,  52, 194, 110),
	SEMI_LIGHT_P_TYPE	("Lightly doped P-type semiconductor",	201, 131,  73, 110),
	SEMI_LIGHT_N_TYPE	("Lightly doped N-type semiconductor",	137, 188, 204, 110),
	DIELECTRIC			("Dielectric",							 81, 171,  51, 80),
	FERROMAGNET			("Ferromagnet",							116, 50, 117, 80),
	POS_CHARGE			("Positive static charge",				116, 50, 50, 80),
	NEG_CHARGE			("Negative static charge",				50, 50, 117, 80),
	DECO				("Decoration",							255, 255, 255, 255),
	ABSORBER			("Absorber",							 50,  50,  50),
	VACUUM				("Vacuum",								 20,  20,  20),
	CUSTOM				("Custom material",						 120,  120, 120, 120);

	String name;
	int color_r;
	int color_g;
	int color_b;
	int color_grayscale;

	MaterialType(String name, int r, int g, int b) {
		this.name = name;
		color_r = r;
		color_g = g;
		color_b = b;
		color_grayscale = (int)(0.7*Math.max(Math.max(r, g), b));
	}

	MaterialType(String name, int r, int g, int b, int grayscale_brightness) {
		this.name = name;
		color_r = r;
		color_g = g;
		color_b = b;
		color_grayscale = grayscale_brightness;
	}

	public boolean isConducting() {
		return (this == MaterialType.EMF
				|| this == MaterialType.AC_EMF
				|| this == MaterialType.CURRENT
				|| this == MaterialType.SWITCH
				|| this == MaterialType.METAL
				|| this == MaterialType.METAL_HIGH_W
				|| this == MaterialType.METAL_LOW_W
				|| this == MaterialType.METAL_HIGH_C
				|| this == MaterialType.METAL_LOW_C);
	}

	public boolean isSemiconducting() {
		return (this == MaterialType.SEMI_P_TYPE
				|| this == MaterialType.SEMI_N_TYPE
				|| this == MaterialType.SEMI
				|| this == MaterialType.SEMI_HEAVY_P_TYPE
				|| this == MaterialType.SEMI_HEAVY_N_TYPE
				|| this == MaterialType.SEMI_LIGHT_P_TYPE
				|| this == MaterialType.SEMI_LIGHT_N_TYPE);
	}

	public boolean isInteractable() {
		return (this == MaterialType.EMF
				|| this == MaterialType.AC_EMF
				|| this == MaterialType.CURRENT
				|| this == MaterialType.SWITCH);
	}
	
	public boolean hasEMF() {
		return (this == MaterialType.EMF
				|| this == MaterialType.AC_EMF
				|| this == MaterialType.CURRENT);
	}
	
	public String getName() {
		return name;
	}
}