// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.probe;

public class Ground extends VoltageProbe {
	@Override
	public Ground clone() {
		Ground p = null;
		p = (Ground) super.clone();
		return p;
	}
}