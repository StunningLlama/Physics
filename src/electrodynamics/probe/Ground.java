// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.units.Quantity;
import electrodynamics.units.Units;

public class Ground extends VoltageProbe {
	public Ground(int mx, int my) {
		super(mx, my);
	}

	@Override
	public Ground clone() {
		Ground p = null;
		p = (Ground) super.clone();
		return p;
	}
	
	@Override
	public String getText(Units units) {
		return "Ground = " + units.toString_fixedsigfigs(value - value, Quantity.ELECTRIC_POTENTIAL, 1e-6);
	}
}