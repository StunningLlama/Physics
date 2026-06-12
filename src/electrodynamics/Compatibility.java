// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import electrodynamics.Renderer.ScalarMode;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorMode;
import electrodynamics.Renderer.VectorView;

public class Compatibility {
	enum VectorMode_v1 {
		ARROWS,
		LINES,
		DOTS,
		CONTOUR,
		SPECIES;
		
		public void applySetting(Simulation e) {
			switch(this) {
			case ARROWS:
				e.controls.vectormode.setOption(VectorMode.ARROWS);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case CONTOUR:
				e.controls.vectormode.setOption(VectorMode.NONE);
				e.controls.scalarmode.setOption(ScalarMode.CONTOUR_COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case DOTS:
				e.controls.vectormode.setOption(VectorMode.DOTS);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case LINES:
				e.controls.vectormode.setOption(VectorMode.LINES);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case SPECIES:
				e.controls.vectormode.setOption(VectorMode.NONE);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(true);
				e.opts.gui_carrier_density.setValue(e.opts.gui_brightness_vec.getValue());
				break;
			default:
				break;
			}
		}
	}
	
	public static ScalarView[] scalar_order_v1 = {
			ScalarView.NONE,
			ScalarView.E_FIELD,
			ScalarView.B_FIELD,
			ScalarView.CHARGE,
			ScalarView.CURRENT,
			ScalarView.H_FIELD,
			ScalarView.POTENTIAL,
			ScalarView.ENERGY,
			ScalarView.ELECTRON_CHARGE,
			ScalarView.HOLE_CHARGE,
			ScalarView.COMBINED_CHARGE,
			ScalarView.BACKGROUND_CHARGE,
			ScalarView.HEAT,
			ScalarView.ENTROPY,
			ScalarView.ELECTRON_POTENTIAL,
			ScalarView.HOLE_POTENTIAL,
			ScalarView.AVERAGE_POTENTIAL,
			ScalarView.GENERATION,
			ScalarView.RECOMBINATION,
			ScalarView.LIGHT
	};

	public static VectorView[] vector_order_v1 = {
			VectorView.NONE,
			VectorView.E_FIELD,
			VectorView.D_FIELD,
			VectorView.ELECTRON_CURRENT,
			VectorView.HOLE_CURRENT,
			VectorView.TOTAL_CURRENT,
			VectorView.EMF,
			VectorView.POYNTING
	};
}
