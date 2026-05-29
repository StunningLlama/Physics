// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public enum Units {
	SI(new SIUnits()),
	ESU(new ESU()),
	EMU(new EMU()),
	GAUSS(new GaussianUnits()),
	MIXED(new MixedUnits());
	
	public UnitSystem sys;
	Units(UnitSystem s) {
		this.sys = s;
	}
	
	public String toString(double value, Quantity q, double lowerbound) {
		return toString(Math.abs(value) < lowerbound? 0 : value, q);
	}

	public String toString(double value, Quantity q) {
		return sys.toStringSI(value, q, "%.2f");
	}

	public String toString_fixedsigfigs(double value, Quantity q, double lowerbound) {
		return toString_fixedsigfigs(Math.abs(value) < lowerbound? 0 : value, q);
	}
	
	public String toString_fixedsigfigs(double value, Quantity q) {
		return sys.toStringSI(value, q, ((value > 0)? " " : "") +"%.4g");
	}
	
	@Override
	public String toString() {
		return sys.getName();
	}
	
	private static void test() {
		for (Units u : Units.values()) {
			System.out.println("Testing " + u.toString() + ".");
			
			for (Quantity q : Quantity.values()) {
				System.out.println("1 " + q.name + " (SI) = " + u.toString(1, q));
			}
			
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		Units.test();
		return;
	}
}
