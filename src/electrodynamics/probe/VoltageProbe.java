// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.units.Quantity;

public class VoltageProbe extends PointProbe {
	public VoltageProbe(int mx, int my) {
		super(mx, my);
		shorthand = "V";
		quantity = Quantity.ELECTRIC_POTENTIAL;
	}
	

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		value = e.V_avg[x][y];
		if (e.hasGround() && this != e.getGround())
			value -= e.getGround().value;
		
		if (savedatapoint) data.addData(value, e.time);
	}
}