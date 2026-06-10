// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.units;

public abstract class UnitSystem {
	public double time_SI = 0;
	public double len_SI = 0;
	public double mass_SI = 0;
	public double charge_SI = 0;
	public double temp_SI = 0;
	
	public abstract String getName();
	
	public double toSI(double val, Quantity q) {
		return val*Math.pow(time_SI, q.time)*Math.pow(len_SI, q.len)*Math.pow(mass_SI, q.mass)*Math.pow(charge_SI, q.charge)*Math.pow(temp_SI, q.temp);
	}
	
	public double fromSI(double val, Quantity q) {
		return val*Math.pow(1/time_SI, q.time)*Math.pow(1/len_SI, q.len)*Math.pow(1/mass_SI, q.mass)*Math.pow(1/charge_SI, q.charge)*Math.pow(1/temp_SI, q.temp);
	}
	
	public abstract String getSymbol(Quantity q);
	
	public double getDefaultMagnitude(Quantity q) {return 1;}
	
	
	public String toString(double val, Quantity q, String format) {
		return attachPrefix(val, q, format);
	}
	
	public String toStringSI(double val_SI, Quantity q, String format) {
		return attachPrefix(fromSI(val_SI, q), q, format);
	}
	
	String attachPrefix(double value, Quantity q, String format) {
		String unit = getSymbol(q);
		
		if (q.no_prefixes) {
			return String.format("%.4g", value) + " " + unit;
		}
		
		double default_magnitude = getDefaultMagnitude(q);
		value = value*default_magnitude;
		
		if (!Double.isFinite(value))
			return Double.toString(value) + " " + unit;

		double mag = Math.abs(value);
		
		String precision = format;
		if (mag < 1E-27)
			return "0 " + unit;
		else if (mag < 1E-21)
			return String.format(precision, value*1e24) + " y" + unit;
		else if (mag < 1E-18)
			return String.format(precision, value*1e21) + " z" + unit;
		else if (mag < 1E-15)
			return String.format(precision, value*1e18) + " a" + unit;
		else if (mag < 1E-12)
			return String.format(precision, value*1e15) + " f" + unit;
		else if (mag < 1E-9)
			return String.format(precision, value*1e12) + " p" + unit;
		else if (mag < 1E-6)
			return String.format(precision, value*1e9) + " n" + unit;
		else if (mag < 1E-3)
			return String.format(precision, value*1e6) + " \u00b5" + unit;
		else if (mag < 1)
			return String.format(precision, value*1e3) + " m" + unit;
		else if (mag < 1E3)
			return String.format(precision, value) + " " + unit;
		else if (mag < 1E6)
			return String.format(precision, value*1e-3) + " k" + unit;
		else if (mag < 1E9)
			return String.format(precision, value*1e-6) + " M" + unit;
		else if (mag < 1E12)
			return String.format(precision, value*1e-9) + " G" + unit;
		else if (mag < 1E15)
			return String.format(precision, value*1e-12) + " T" + unit;
		else if (mag < 1E18)
			return String.format(precision, value*1e-15) + " P" + unit;
		else if (mag < 1E21)
			return String.format(precision, value*1e-18) + " E" + unit;
		else if (mag < 1E27)
			return String.format(precision, value*1e-21) + " Z" + unit;
		else if (mag < 1E30)
			return String.format(precision, value*1e-23) + " Y" + unit;
		else
			return "infinity " + unit;
	}
}