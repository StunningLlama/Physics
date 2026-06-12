// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

import electrodynamics.Simulation;
import electrodynamics.Renderer.VectorView;
import electrodynamics.units.Quantity;

public class CurrentProbe extends LineProbe {
	public CurrentProbe(int mx, int my) {
		super(mx, my);
		shorthand = "I";
		quantity = Quantity.ELECTRIC_CURRENT;
		vectorname = VectorView.TOTAL_CURRENT;
	}

	@Override
	public void measure(Simulation e, boolean savedatapoint) {
		super.measure(e, savedatapoint);
	}
}