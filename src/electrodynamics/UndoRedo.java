// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.Adjustable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.swing.AbstractButton;

import electrodynamics.Renderer.ScalarMode;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorMode;
import electrodynamics.Renderer.VectorView;
import electrodynamics.Simulation.BoundaryCondition;
import electrodynamics.probe.Probe;
import electrodynamics.util.Utils;

public class UndoRedo {
	
	public List<Snapshot> prev_states = new ArrayList<Snapshot>();
	public int undoredo_pointer = 0;
	public int history_size = 0;
	
	public UndoRedo(int history_size) {
		this.history_size = history_size;
	}

	public void resetUndoHistory(Simulation e) {
		undoredo_pointer = 0;
		prev_states.clear();

		Snapshot state = new Snapshot();
		state.store(e);
		prev_states.add(state);
	}

	public void undo(Simulation e) {
		undoredo_pointer--;
		if (undoredo_pointer < 0) undoredo_pointer = 0;

		if (undoredo_pointer >= 0 && undoredo_pointer < prev_states.size()) {
			prev_states.get(undoredo_pointer).load(e);
		}
	}

	public void redo(Simulation e) {
		undoredo_pointer++;
		if (undoredo_pointer >= prev_states.size()) undoredo_pointer = prev_states.size()-1;

		if (undoredo_pointer >= 0 && undoredo_pointer < prev_states.size())
			prev_states.get(undoredo_pointer).load(e);
	}
	
	public boolean canRedo() {
		return undoredo_pointer < prev_states.size()-1;
	}
	
	public boolean canUndo() {
		return undoredo_pointer > 0;
	}
	
	public void captureState(Simulation e) {
		int i = undoredo_pointer+1;
		
		while (i < prev_states.size()) {
			prev_states.remove(i);
		}
		
		undoredo_pointer++;
		Snapshot state = new Snapshot();
		state.store(e);
		prev_states.add(state);

		while (prev_states.size() > history_size)
		{
			prev_states.remove(0);
			undoredo_pointer--;
		}
	}
}

class Snapshot {
	int resolution;
	double width;
	double time;
	
	HashMap<AbstractButton, Boolean> boolean_values;
	HashMap<Adjustable, Integer> integer_values;
	
	boolean gui_paused;
	boolean gui_tooltip;
	boolean gui_text_bg;
	boolean gui_elem_colors;
	boolean gui_interface;
	boolean gui_carriers;
	boolean gui_hide_carriers_metal;
	boolean gui_carrier_diffusion;
	boolean gui_gen_recomb;
	boolean gui_border;
	int gui_simspeed;
	int gui_simspeed_2;
	int gui_brightness;
	int gui_brightness_vec;
	int gui_carrier_number;
	
	String description;
	ScalarView gui_view;
	VectorView gui_view_vec;
	ScalarMode gui_scalar_mode;
	VectorMode gui_view_vec_mode;
	BoundaryCondition gui_bc;
	int gui_parameter1;
	
	double[][] ex;
	double[][] ey;
	double[][] hz;
	double[][] rho_c;
	double[][] rho_n;
	double[][] rho_p;
	double[][] rho_back;
	double[][] rho_free;
	double[][] jx_c;
	double[][] jy_c;
	double[][] jx_n;
	double[][] jy_n;
	double[][] jx_p;
	double[][] jy_p;
	Material[][] materials;

	List<Probe> probes;

	public void store(Simulation e) {
		time = e.time;

		gui_view = e.controls.scalarview.getOption();
		gui_view_vec = e.controls.vectorview.getOption();
		gui_scalar_mode = e.controls.scalarmode.getOption();
		gui_view_vec_mode = e.controls.vectormode.getOption();
		gui_bc = (BoundaryCondition) e.opts.gui_bc.getSelectedItem();
		description = e.description;

		boolean_values = new HashMap<AbstractButton, Boolean>();
		integer_values = new HashMap<Adjustable, Integer>();

		for (AbstractButton item : e.opts.boolean_names.values()) {
			boolean_values.put(item, item.isSelected());
		}

		for (Adjustable item : e.opts.integer_names.values()) {
			integer_values.put(item, item.getValue());
		}

		ex = copy(e.Ex);
		ey = copy(e.Ey);
		hz = copy(e.Hz);
		rho_c = copy(e.rho_abs);
		rho_n = copy(e.rho_n);
		rho_p = copy(e.rho_p);
		rho_back = copy(e.rho_back);
		rho_free = copy(e.rho_free);
		jx_c = copy(e.Jx_abs);
		jy_c = copy(e.Jy_abs);
		jx_n = copy(e.Jx_n);
		jy_n = copy(e.Jy_n);
		jx_p = copy(e.Jx_p);
		jy_p = copy(e.Jy_p);
		materials = copy(e.materials);
		probes = Utils.cloneList(e.probes, Probe::clone);
	}


	public void load(Simulation e) {

		e.rwLock.writeLock().lock();
		try {
			e.time = time;

			if (e.prefs.chkbox_undo.isSelected()) {
				e.description = description;
				e.controls.scalarview.setOption(gui_view);
				e.controls.vectorview.setOption(gui_view_vec);
				e.controls.scalarmode.setOption(gui_scalar_mode);
				e.controls.vectormode.setOption(gui_view_vec_mode);
				e.opts.gui_bc.setSelectedItem(gui_bc);

				for (AbstractButton item : boolean_values.keySet()) {
					item.setSelected(boolean_values.get(item));
				}

				for (Adjustable item : integer_values.keySet()) {
					item.setValue(integer_values.get(item));
				}

				e.opts.setRedundantOptions();
			}

			e.Ex = copy(ex); 
			e.Ey = copy(ey);
			e.Hz = copy(hz);

			e.rho_abs = copy(rho_c);
			e.rho_n = copy(rho_n);
			e.rho_p = copy(rho_p);
			e.rho_back = copy(rho_back);
			e.rho_free = copy(rho_free);

			e.Jx_abs = copy(jx_c);
			e.Jy_abs = copy(jy_c);
			e.Jx_n = copy(jx_n);
			e.Jy_n = copy(jy_n);
			e.Jx_p = copy(jx_p);
			e.Jy_p = copy(jy_p);

			e.materials = copy(materials);
			e.probes = Utils.cloneList(probes, Probe::clone);

			e.opts.textPane.setText(e.description);
			e.opts.textPane.setEditable(false);
			e.opts.textPane.setCaretPosition(0);
			e.updateAllMaterials(false);
			e.calcMiscFields(true);
		}
		finally {
			e.rwLock.writeLock().unlock();
		}
	}
	
	public static double[][] copy(double[][] arr) {
	    if (arr == null) {
	        return null;
	    }

	    final double[][] result = new double[arr.length][];
	    for (int i = 0; i < arr.length; i++) {
	        result[i] = Arrays.copyOf(arr[i], arr[i].length);
	    }
	    return result;
	}
	
	public static Material[][] copy(Material[][] arr) {
	    if (arr == null) {
	        return null;
	    }

	    final Material[][] result = new Material[arr.length][];
	    for (int i = 0; i < arr.length; i++) {
	        result[i] = new Material[arr[i].length];
		    for (int j = 0; j < arr[i].length; j++) {
		    	result[i][j] = new Material(arr[i][j]);
		    }
	    }
	    return result;
	}
}
