// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ButtonGroup;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTextArea;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorView;
import electrodynamics.plot.Plot;
import electrodynamics.util.Font7x5;
import electrodynamics.util.Utils;
import electrodynamics.util.Vector;

public class Controls implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener, ItemListener {
	Simulation e;
	
	/* Keyboard controls */

	public boolean advanceframe = false;
	public boolean clear = false;
	public boolean reset = false;
	public boolean save = false;
	public boolean load = false;
	public boolean debugging = false;
	public boolean cut = false;
	public boolean copy = false;
	public boolean paste = false;
	public boolean undo = false;
	public boolean redo = false;
	public boolean delete = false;
    public boolean shift_down = false;
    public boolean ctrl_down = false;
    public boolean alt_down = false;
    public boolean logdata = false;
    public boolean rotate_selection = false;
    public boolean flip_h_selection = false;
    public boolean flip_v_selection = false;


	/* Mouse controls */

	PointerInfo pointerinfo = MouseInfo.getPointerInfo();
	public boolean mouse_pressed = false;
	public boolean mouse_pressed_prev = false;
	public boolean moving_selection = false;
	public boolean dragging_selection = false;
	public boolean brush_changed = false;

	public int mousebutton = 0;
	public int mx = 0;
	public int my = 0;
	public int mx_start = 0;
	public int my_start = 0;

	public int mx_index = 0;
	public int my_index = 0;
	public int mx_start_index = 0;
	public int my_start_index = 0;

	public double mx_realspace = 0;
	public double my_realspace = 0;
	public double mxp_realspace = 0;
	public double myp_realspace = 0;
	public double mx_start_realspace = 0;
	public double my_start_realspace = 0;

	public int delta_mx_index = 0;
	public int delta_my_index = 0;

	public boolean EMF_selected = false;
	public double max_EMF = 5e5;

	Brush prev_brush;
	public double brushsize = 0;
	public int prev_EMF_setting = 0;
	BoundaryCondition prev_boundary = null;

	public boolean[][] under_brush;
	public boolean[][] selected;
	public boolean[][] selected_EMF;

	public ClipboardMaterial[][] selection;
	public ClipboardMaterial[][] clipboard;

	public int text_x = 0;
	public int text_y = 0;
	public boolean texting = false;

	public double flashlight_strength = 1e31;

	Cursor HAND_CURSOR = new Cursor(Cursor.HAND_CURSOR);
	Cursor DEFAULT_CURSOR = new Cursor(Cursor.DEFAULT_CURSOR);

	public HashMap<Brush, JRadioButtonMenuItem> brushbuttonmap;
	public ButtonGroup buttongroup;
	
	public List<Snapshot> prev_states = new ArrayList<Snapshot>();
	public int undoredo_pointer = 0;
	
	public Controls(Simulation e) {
		this.e = e;
	}
	
	void setResolution(int resolution) {
		selection = new ClipboardMaterial[e.nx][e.ny];
		clipboard = new ClipboardMaterial[e.nx][e.ny];

		under_brush = new boolean[e.nx][e.ny];
		selected = new boolean[e.nx][e.ny];
		selected_EMF = new boolean[e.nx][e.ny];
		
		text_x = 0;
		text_y = 0;
		endTextInput();
	}
	
	public void handleMouseInput() {
		Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();
		BrushShape brushshape = (BrushShape) e.opts.gui_brush_1.getSelectedItem();

		boolean pressing = false;
		boolean releasing = false;

		if (mouse_pressed) {
			if (!mouse_pressed_prev) {
				pressing = true;
				e.canvas.requestFocus();
				if (Brush.isMaterialModifyingBrush(brush) || brush == Brush.SELECT)
					e.opts.gui_paused.setSelected(true);
			}
		} else {
			if (mouse_pressed_prev) {
				releasing = true;
			}
		}
		mouse_pressed_prev = mouse_pressed;

		mx_realspace = Math.round((mx-1)/e.renderer.scalefactor_real - 0.5)*e.ds;
		my_realspace = Math.round((my-1)/e.renderer.scalefactor_real - 0.5)*e.ds;

		mx_start_realspace =  Math.round((mx_start-1)/e.renderer.scalefactor_real - 0.5)*e.ds;
		my_start_realspace = Math.round((my_start-1)/e.renderer.scalefactor_real - 0.5)*e.ds;

		mx_index = (int)Math.round((mx-1)/e.renderer.scalefactor_real - 0.5);
		my_index = (int)Math.round((my-1)/e.renderer.scalefactor_real - 0.5);

		mx_start_index = (int)Math.round((mx_start-1)/e.renderer.scalefactor_real - 0.5);
		my_start_index = (int)Math.round((my_start-1)/e.renderer.scalefactor_real - 0.5);

		if (mx_index < 0) mx_index = 0;
		if (my_index < 0) my_index = 0;
		if (mx_index >= e.nx) mx_index = e.nx-1;
		if (my_index >= e.ny) my_index = e.ny-1;

		if (mx_start_index < 0) mx_start_index = 0;
		if (my_start_index < 0) my_start_index = 0;
		if (mx_start_index >= e.nx) mx_start_index = e.nx-1;
		if (my_start_index >= e.ny) my_start_index = e.ny-1;

		e.opts.gui_stepsizelbl.setText("Timestep: " + Utils.getSI(e.dt, "s"));
		e.opts.gui_stepslbl.setText("Sim steps/frame: " + e.opts.gui_simspeed_2.getValue());

		brushsize = e.ds*(Math.pow(10.0, 2*e.opts.gui_brushsize.getValue()/(50.0*10.0) - 0.75) + e.opts.gui_brushsize.getValue()/10.0 + 0.5);
		e.opts.lblBrushSize.setText("Brush size: " + (int)Math.ceil(brushsize/e.ds));

		if (!Brush.isMaterialModifyingBrush(brush))
		{
			e.opts.gui_parameter2.setVisible(false);
			e.opts.gui_parameter2_text.setVisible(false);
			e.opts.gui_parameter2_text.setText("");
		}


		if (Brush.isBrushShapeImportant(brush)) {
			e.opts.gui_brush_1.setVisible(true);
			e.opts.gui_brush_highlight.setVisible(true);
			e.opts.gui_brushsize.setVisible(true);
			e.opts.lblBrushSize.setVisible(true);
		} else {
			e.opts.gui_brush_1.setVisible(false);
			e.opts.gui_brush_highlight.setVisible(false);
			e.opts.gui_brushsize.setVisible(false);
			e.opts.lblBrushSize.setVisible(false);
		}


		if (Brush.isMaterialModifyingBrush(brush) && brush != Brush.ERASE) {
			e.opts.gui_material.setVisible(true);
		} else {
			e.opts.gui_material.setVisible(false);
		}

		if (brush_changed) {
			if (!(brush == Brush.SELECT || brush == Brush.FLOODSELECT)) {
				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						selected[i][j] = false;
					}
				}
			}

			if (brush != Brush.INTERACT) {
				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						selected_EMF[i][j] = false;
					}
				}
				EMF_selected = false;
			}

			brush_changed = false;
		}
		

		boolean update = false;

		if (cut || copy) {
			int i_min = e.nx-1;
			int j_min = e.ny-1;
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					clipboard[i][j].erase();
					if (selected[i][j] && e.materials[i][j].type != MaterialType.VACUUM) {
						if (i < i_min) i_min = i;
						if (j < j_min) j_min = j;
					}
				}
			}
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					if (selected[i][j] && e.materials[i][j].type != MaterialType.VACUUM) {
						
						clipboard[i-i_min][j-j_min] = new ClipboardMaterial(e, i, j);
						
						if (cut) {
							e.eraseMaterial(i, j);
						}
					}
					if (cut) {
						selected[i][j] = false;
					}
				}
			}

			if (cut) update = true;

			cut = false;
			copy = false;
		}

		if (paste) {
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					selected[i][j] = false;
				}
			}

			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					selection[i][j] = clipboard[i][j].clone();
				}
			}
			moving_selection = true;
			dragging_selection = false;
			paste = false;
			update = true;
		}

		if (delete) {
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					if (selected[i][j] && e.materials[i][j].type != MaterialType.VACUUM) {
						e.eraseMaterial(i, j);
					}
				}
			}
			delete = false;
			update = true;
		}
		
		if (rotate_selection || flip_h_selection || flip_v_selection) {
			int i_max = 0;
			int j_max = 0;
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					if (selection[i][j].m.type != MaterialType.VACUUM) {
						if (i > i_max) i_max = i;
						if (j > j_max) j_max = j;
					}
				}
			}
			
			if (rotate_selection) {
				ClipboardMaterial[][] new_selection = new ClipboardMaterial[j_max+1][i_max+1];
				for (int i = 0; i <= i_max; i++)
				{
					for (int j = 0; j <= j_max; j++)
					{
						new_selection[j_max-j][i] = selection[i][j].clone();
						new_selection[j_max-j][i].m.emf_direction += Math.PI/2.0;
						selection[i][j].erase();
					}
				}
				
				for (int i = 0; i <= j_max; i++)
				{
					for (int j = 0; j <= i_max; j++)
					{
						selection[i][j] = new_selection[i][j];
					}
				}
				rotate_selection = false;
			}
			
			if (flip_h_selection) {
				ClipboardMaterial[][] new_selection = new ClipboardMaterial[i_max+1][j_max+1];
				for (int i = 0; i <= i_max; i++)
				{
					for (int j = 0; j <= j_max; j++)
					{
						new_selection[i_max-i][j] = selection[i][j].clone();
						new_selection[i_max-i][j].m.emf_direction = Math.PI - new_selection[i_max-i][j].m.emf_direction;
						selection[i][j].erase();
					}
				}
				
				for (int i = 0; i <= i_max; i++)
				{
					for (int j = 0; j <= j_max; j++)
					{
						selection[i][j] = new_selection[i][j];
					}
				}
				flip_h_selection = false;
			}
			
			if (flip_v_selection) {
				ClipboardMaterial[][] new_selection = new ClipboardMaterial[i_max+1][j_max+1];
				for (int i = 0; i <= i_max; i++)
				{
					for (int j = 0; j <= j_max; j++)
					{
						new_selection[i][j_max - j] = selection[i][j].clone();
						new_selection[i][j_max - j].m.emf_direction = -new_selection[i][j_max - j].m.emf_direction;
						selection[i][j].erase();
					}
				}
				
				for (int i = 0; i <= i_max; i++)
				{
					for (int j = 0; j <= j_max; j++)
					{
						selection[i][j] = new_selection[i][j];
					}
				}
				flip_v_selection = false;
			}
		}
		


		switch(brush) {
		case DRAW:
		case LINE:
		case REPLACE:
		case ERASE:
		case FILL:
		case LIGHT:

			if ((mousebutton == MouseEvent.BUTTON2 || alt_down) && pressing) {
				e.opts.gui_material.setSelectedItem(e.materials[mx_index][my_index].type);
			}

			double angle = 0;
			MaterialType mat = (MaterialType) e.opts.gui_material.getSelectedItem();

			if (mat == MaterialType.EMF || mat == MaterialType.AC_EMF) {
				e.opts.gui_parameter2.setVisible(true);
				e.opts.gui_parameter2_text.setVisible(true);

				int directionval = e.opts.gui_parameter2.getValue()/6;
				if (directionval == 0) {
					e.opts.gui_parameter2_text.setText("EMF direction: Up");
					angle = -Math.PI/2;
				}
				if (directionval == 1) {
					e.opts.gui_parameter2_text.setText("EMF direction: Right");
					angle = 0;
				}
				if (directionval == 2) {
					e.opts.gui_parameter2_text.setText("EMF direction: Down");
					angle = Math.PI/2;
				}
				if (directionval == 3) {
					e.opts.gui_parameter2_text.setText("EMF direction: Left");
					angle = Math.PI;
				}
				//e.opts.gui_parameter2_text.setText("Brush orientation: " + directionval*(360/24) + " deg");
				//angle = Math.PI * directionval/12.0;
			} else {
				e.opts.gui_parameter2.setVisible(false);
				e.opts.gui_parameter2_text.setVisible(false);
				e.opts.gui_parameter2_text.setText("");
			}


			if (mousebutton == MouseEvent.BUTTON3 || brush == Brush.ERASE)
				mat = MaterialType.VACUUM;

			if (!(mousebutton == MouseEvent.BUTTON2 || alt_down)) {
				if (brush == Brush.LINE) {
					if (releasing) {
						drawMaterialLine(mx_start_realspace, my_start_realspace, mx_realspace, my_realspace, brush, brushshape, mat, brushsize, angle);
					}
				} else if (brush == Brush.FILL) {
					if (pressing) {
						MaterialType old_mat = e.materials[mx_index][my_index].type;
						MaterialType new_mat = mat;
						double new_angle = angle;
						if (new_mat != old_mat) {
							this.floodFill(mx_index, my_index, new FloodFillFunc() {
								@Override
								public boolean isValid(int i, int j) {
									return e.materials[i][j].type == old_mat;
								}

								@Override
								public void fill(int i, int j) {
									e.eraseMaterial(i, j);
									e.initializeMaterial(i, j, new_mat);
									if (new_mat == MaterialType.EMF || new_mat == MaterialType.AC_EMF) e.materials[i][j].emf_direction = new_angle;
								}
							});
						}
					}
				} else if (brush == Brush.LIGHT) {
					if (mouse_pressed) {
						for (int i = 0; i < e.nx; i++)
						{
							for (int j = 0; j < e.ny; j++)
							{
								double cx = 0;
								double cy = 0;
								cx = i*e.ds;
								cy = j*e.ds;

								double px = (cx-mx_realspace);
								double py = (cy-my_realspace);
								double r = 0;

								if (brushshape == BrushShape.CIRCLE)
									r = Math.sqrt(px*px+py*py);
								else if (brushshape == BrushShape.SQUARE)
									r = Math.max(Math.abs(px), Math.abs(py));
								e.L[i][j] = (r <= brushsize)? flashlight_strength : 0;
							}
						}
					} else if (releasing) {

						for (int i = 0; i < e.nx; i++)
						{
							for (int j = 0; j < e.ny; j++)
							{
								e.L[i][j] = 0;
							}
						}
					}
				}
				else if (mouse_pressed) {
					drawMaterialLine(mxp_realspace, myp_realspace, mx_realspace, my_realspace, brush, brushshape, mat, brushsize, angle);
				}
			}

			if (Brush.isBrushShapeImportant(brush)) {
				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						double cx = 0;
						double cy = 0;
						cx = i*e.ds;
						cy = j*e.ds;

						double px = (cx-mx_realspace);
						double py = (cy-my_realspace);
						double r = 0;

						if (brushshape == BrushShape.CIRCLE)
							r = Math.sqrt(px*px+py*py);
						else if (brushshape == BrushShape.SQUARE)
							r = Math.max(Math.abs(px), Math.abs(py));
						under_brush[i][j] = (r <= brushsize);
					}
				}
			}

			break;

		case INTERACT:
			if (e.materials[mx_index][my_index].type == MaterialType.EMF || e.materials[mx_index][my_index].type == MaterialType.AC_EMF || e.materials[mx_index][my_index].type == MaterialType.SWITCH)
				e.canvas.setCursor(HAND_CURSOR);
			else
				e.canvas.setCursor(DEFAULT_CURSOR);

			if (pressing) {
				boolean turn_on_EMF = !selected_EMF[mx_index][my_index];

				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						selected_EMF[i][j] = false;
					}
				}
				EMF_selected = false;

				if ((e.materials[mx_index][my_index].type == MaterialType.EMF || e.materials[mx_index][my_index].type == MaterialType.AC_EMF) && turn_on_EMF) {
					this.floodFill(mx_index, my_index, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == e.materials[mx_index][my_index].type && selected_EMF[i][j] != true;
						}

						@Override
						public void fill(int i, int j) {
							selected_EMF[i][j] = true;
						}
					});
					
					EMF_selected = true;
					int setting = (int)(Math.round(50*e.materials[mx_index][my_index].emf/max_EMF));
					e.opts.gui_parameter3.setValue(setting);
					prev_EMF_setting = setting;
				}

				if (e.materials[mx_index][my_index].type == MaterialType.SWITCH) {
					int active = 1-e.materials[mx_index][my_index].activated;
					
					this.floodFill(mx_index, my_index, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == MaterialType.SWITCH && e.materials[i][j].activated != active;
						}

						@Override
						public void fill(int i, int j) {
							e.materials[i][j].activated = active;
						}
					});
				
					update = true;
				}
			}
			break;
		case FLOODSELECT:
		case SELECT:
			if (pressing) {
				if (brush == Brush.FLOODSELECT && !moving_selection) {

					this.floodFill(mx_index, my_index, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == e.materials[mx_index][my_index].type && selected[i][j] != !selected[mx_index][my_index];
						}

						@Override
						public void fill(int i, int j) {
							selected[i][j] = !selected[mx_index][my_index];
						}
					});
				}
				else if (moving_selection && !dragging_selection) {
					for (int i = 0; i < e.nx; i++)
					{
						for (int j = 0; j < e.ny; j++)
						{
							int si = i-delta_mx_index;
							int sj = j-delta_my_index;
							if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && selection[si][sj].m.type != MaterialType.VACUUM) {
								e.eraseMaterial(i, j);
								selection[si][sj].paste(e, i, j);
								selected[i][j] = true;
							}
						}
					}
					update = true;
					moving_selection = false;
					dragging_selection = false;
				} else if (!moving_selection && selected[mx_index][my_index]) {
					for (int i = 0; i < e.nx; i++)
					{
						for (int j = 0; j < e.ny; j++)
						{
							selection[i][j].erase();
							if (selected[i][j]) {
								selection[i][j] = new ClipboardMaterial(e, i, j);
								selected[i][j] = false;
								e.eraseMaterial(i, j);
							}
						}
					}
					update = true;
					moving_selection = true;
					dragging_selection = true;
					delta_mx_index = 0;
					delta_my_index = 0;
				}
			} else if (mouse_pressed) {
				if (brush != Brush.FLOODSELECT) {
					if (dragging_selection) {
						delta_mx_index = mx_index - mx_start_index;
						delta_my_index = my_index - my_start_index;
					} else {
						int mx0 = Math.min(mx_start_index, mx_index);
						int my0 = Math.min(my_start_index, my_index);
						int mx1 = Math.max(mx_start_index, mx_index);
						int my1 = Math.max(my_start_index, my_index);

						for (int i = 0; i < e.nx; i++)
						{
							for (int j = 0; j < e.ny; j++)
							{
								if (i >= mx0 && i <= mx1 && j >= my0 && j <= my1)
									selected[i][j] = true;
								else
									selected[i][j] = false;
							}
						}
					}
				}
			} else if (releasing) {
				if (brush != Brush.FLOODSELECT) {
					delta_mx_index = mx_index - mx_start_index;
					delta_my_index = my_index - my_start_index;
					if (dragging_selection) {
						for (int i = 0; i < e.nx; i++)
						{
							for (int j = 0; j < e.ny; j++)
							{
								int si = i-delta_mx_index;
								int sj = j-delta_my_index;
								if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && selection[si][sj].m.type != MaterialType.VACUUM) {
									e.eraseMaterial(i, j);
									selection[si][sj].paste(e, i, j);
									selected[i][j] = true;
								}
							}
						}
						moving_selection = false;
						dragging_selection = false;
						update = true;
					} else {
						if (delta_mx_index == 0 && delta_my_index == 0) {
							for (int i = 0; i < e.nx; i++)
							{
								for (int j = 0; j < e.ny; j++)
								{
									selected[i][j] = false;
								}
							}
						}
					}
				}
			} else {
				delta_mx_index = mx_index;
				delta_my_index = my_index;
			}
			break;
		case CURRENT:
			if (pressing) {
				CurrentProbe p = new CurrentProbe();
				p.x1 = mx_start_index;
				p.y1 = my_start_index;
				p.x2 = mx_index;
				p.y2 = my_index;
				e.currentprobes.add(p);
			} else if (mouse_pressed) {
				e.currentprobes.get(e.currentprobes.size()-1).x2 = mx_index;
				e.currentprobes.get(e.currentprobes.size()-1).y2 = my_index;
			}
			break;
		case VOLTAGE:
			if (pressing) {
				VoltageProbe p = new VoltageProbe();
				p.x = mx_start_index;
				p.y = my_start_index;
				e.voltageprobes.add(p);
			} else if (mouse_pressed) {
				e.voltageprobes.get(e.voltageprobes.size()-1).x = mx_index;
				e.voltageprobes.get(e.voltageprobes.size()-1).y = my_index;
			}
			break;
		case CHARGE:
				if (pressing) {
					ChargeProbe p = new ChargeProbe();
					p.x1 = mx_start_index;
					p.y1 = my_start_index;
					p.x2 = mx_index;
					p.y2 = my_index;
					e.chargeprobes.add(p);
				} else if (mouse_pressed) {
					e.chargeprobes.get(e.chargeprobes.size()-1).x2 = mx_index;
					e.chargeprobes.get(e.chargeprobes.size()-1).y2 = my_index;
				}
			break;
		case DELETEPROBE:
			e.canvas.setCursor(DEFAULT_CURSOR);
			int i = 0;
			while(i < e.voltageprobes.size()) {
				VoltageProbe p = e.voltageprobes.get(i);
				if (Utils.length(p.x-mx_index, p.y-my_index) < 3) {
					e.canvas.setCursor(HAND_CURSOR);
					if (pressing) {
						e.voltageprobes.remove(i);
						i--;
					}
				}
				i++;
			}

			i = 0;
			while(i < e.currentprobes.size()) {
				CurrentProbe p = e.currentprobes.get(i);
				if (Utils.length(p.x1-mx_index, p.y1-my_index) < 3 || Utils.length(p.x2-mx_index, p.y2-my_index) < 3) {
					e.canvas.setCursor(HAND_CURSOR);
					if (pressing) {
						e.currentprobes.remove(i);
						i--;
					}
				}
				i++;
			}

			i = 0;
			while(i < e.chargeprobes.size()) {
				ChargeProbe p = e.chargeprobes.get(i);
				if (mx_index >= Math.min(p.x1, p.x2) && mx_index <= Math.max(p.x1, p.x2) && my_index >= Math.min(p.y1, p.y2) && my_index <= Math.max(p.y1, p.y2)) {
					e.canvas.setCursor(HAND_CURSOR);
					if (pressing) {
						e.chargeprobes.remove(i);
						i--;
					}
				}
				i++;
			}
			
			if (e.ground != null && Utils.length(e.ground.x-mx_index, e.ground.y-my_index) < 3) {
				e.canvas.setCursor(HAND_CURSOR);
				if (pressing) {
					e.ground = null;
				}
			}

			break;
		case TEXT:
			e.canvas.setCursor(HAND_CURSOR);
			if (mouse_pressed) {
				startTextInput();
				text_x = mx_index;
				text_y = my_index;
			}
			break;
		case GROUND:
			if (pressing) {
				if (e.ground == null)
					e.ground = new VoltageProbe();
				e.ground.x = mx_start_index;
				e.ground.y = my_start_index;
			} else if (mouse_pressed) {
				e.ground.x = mx_index;
				e.ground.y = my_index;
			}
			break;
		case BANDS:
			if (releasing) e.bandplot.createPlot(e);
			break;
		case SCALARPLOT:
			if (releasing) e.scalarplot.createPlot(e);
			break;
		case CARRIERPLOT:
			if (releasing) e.carrierplot.createPlot(e);
			break;
		default:
			break;
		}

		if (brush != Brush.TEXT)
		{
			endTextInput();
		}
		
		SwingUtilities.invokeLater(() -> { //TODO
			for (Plot p : e.plots) {
				p.updatePlot(e);
			}
		});


		setEMFs();

		if ((releasing && Brush.isMaterialModifyingBrush(brush)) || (BoundaryCondition)e.opts.gui_bc.getSelectedItem() != prev_boundary || update) {

			//e.resetFields(false);
			e.updateAllMaterials(false);
			e.multigridSolve(true, false);
			
			captureState();
		}

		prev_boundary = (BoundaryCondition)e.opts.gui_bc.getSelectedItem();

		mxp_realspace = mx_realspace;
		myp_realspace = my_realspace;
	}
	
	public void startTextInput() {
		if (!texting) {
			removeKeyBinds(e.canvas);
			removeKeyBinds(e.opts.panel);
			texting = true;
		}
	}

	public void endTextInput() {
		if (texting) {
			addKeyBinds(e.canvas);
			addKeyBinds(e.opts.panel);
			texting = false;
		}
	}
	
	public void resetUndoHistory() {
		undoredo_pointer = 0;
		prev_states.clear();

		Snapshot state = new Snapshot();
		state.store(e);
		prev_states.add(state);
	}
	
	public void handleUndoRedo() {
		if (undo) {
			undoredo_pointer--;
			if (undoredo_pointer < 0) undoredo_pointer = 0;
			
			if (undoredo_pointer >= 0 && undoredo_pointer < prev_states.size()) {
				prev_states.get(undoredo_pointer).load(e);
			}
			
			undo = false;
		}
		
		if (redo) {
			undoredo_pointer++;
			if (undoredo_pointer >= prev_states.size()) undoredo_pointer = prev_states.size()-1;
			
			if (undoredo_pointer >= 0 && undoredo_pointer < prev_states.size())
				prev_states.get(undoredo_pointer).load(e);

			redo = false;
		}
	}
	
	public void captureState() {
		int i = undoredo_pointer+1;
		
		while (i < prev_states.size()) {
			prev_states.remove(i);
		}
		
		undoredo_pointer++;
		Snapshot state = new Snapshot();
		state.store(e);
		prev_states.add(state);

		while (prev_states.size() > 4)
		{
			prev_states.remove(0);
			undoredo_pointer--;
		}
	}

	public void drawMaterialLine(double x1, double y1, double x2, double y2, Brush brush, BrushShape brushshape, MaterialType mat, double brushsize, double EMF_angle) {
		Vector a = new Vector(0, 0);
		Vector b = new Vector(0, 0);
		Vector p = new Vector(0, 0);
		Vector ab = new Vector(0, 0);
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				double cx = 0;
				double cy = 0;
				cx = i*e.ds;
				cy = j*e.ds;

				a.initialize(x1, y1);
				b.initialize(x2, y2);
				p.initialize(cx, cy);
				p.addmult(a, -1);
				ab.copy(b);
				ab.addmult(a, -1);

				double l2 = ab.dot(ab);
				if (l2 == 0)
					l2 = 1;
				double t = Utils.clamp(p.dot(ab)/l2, 0, 1);
				ab.scalarmult(t);
				p.addmult(ab, -1);
				double r = 0;
				if (brushshape == BrushShape.CIRCLE)
					r = Math.sqrt(p.dot(p));
				else if (brushshape == BrushShape.SQUARE)
					r = Math.max(Math.abs(p.x), Math.abs(p.y));
				if (r <= brushsize) {
					if (mat == MaterialType.VACUUM) {
						e.eraseMaterial(i, j);
					} else if (e.materials[i][j].type == MaterialType.VACUUM || brush == Brush.REPLACE) {
						e.eraseMaterial(i, j);
						e.initializeMaterial(i, j, mat);
						if (mat == MaterialType.EMF || mat == MaterialType.AC_EMF) e.materials[i][j].emf_direction = EMF_angle;
					}
				}
			}
		}
	}
	
	public void setEMFs() {
		int EMF_setting = e.opts.gui_parameter3.getValue();
		double new_EMF = max_EMF*EMF_setting/50.0;

		if (e.opts.gui_brush.getSelectedItem() == Brush.INTERACT && EMF_selected) {
			e.opts.gui_parameter3.setVisible(true);
			e.opts.gui_parameter3_text.setVisible(true);
			e.opts.gui_parameter3_text.setText("EMF: " + Utils.getSI(new_EMF, "V/m"));
		} else {
			e.opts.gui_parameter3.setVisible(false);
			e.opts.gui_parameter3_text.setVisible(false);
			e.opts.gui_parameter3_text.setText("");
		}

		if (EMF_setting != prev_EMF_setting && EMF_selected) {
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					if ((e.materials[i][j].type == MaterialType.EMF || e.materials[i][j].type == MaterialType.AC_EMF) && selected_EMF[i][j]) {
						e.materials[i][j].emf = new_EMF;
					}
				}
			}

			e.updateJustEMFs();
		}

		prev_EMF_setting = EMF_setting;
		

		e.AC_freq = 1e13*Math.pow(10, e.opts.gui_parameter1.getValue()/10.0);
		if (e.AC_source_exists) {
			e.opts.gui_parameter1.setEnabled(true);
			e.opts.gui_parameter1.setVisible(true);
			e.opts.gui_parameter1_text.setVisible(true);
		} else {
			e.opts.gui_parameter1.setEnabled(false);
			e.opts.gui_parameter1.setVisible(false);
			e.opts.gui_parameter1_text.setVisible(false);
		}
		
		e.opts.gui_parameter1_text.setText("AC Freq: " + Utils.getSI(e.AC_freq, "Hz"));

	}

	public void floodFill(int i, int j, FloodFillFunc f) {
		Queue<FloodFillCoordinate> queue = new LinkedList<>();
		queue.add(new FloodFillCoordinate(i, j));

		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < e.nx && coord.j >= 0 && coord.j < e.ny && f.isValid(coord.i, coord.j)) {
				f.fill(coord.i, coord.j);
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}

	@Override
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == e.opts.gui_reset)
			clear = true;
		else if (ev.getSource() == e.opts.menu_new)
			reset = true;
		else if (ev.getSource() == e.opts.menu_save)
			save = true;
		else if (ev.getSource() == e.opts.menu_open)
			load = true;
		else if (ev.getSource() == e.opts.menu_help)
			try {
				File helpfile = new File("README.html");
				java.awt.Desktop.getDesktop().browse(helpfile.toURI());
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		else if (ev.getSource() == e.opts.menu_editdesc) {
			SwingUtilities.invokeLater(() -> {
				String text = e.opts.textPane.getText();
				JFrame frame = new JFrame();
				JTextArea area = new JTextArea();
				JButton b = new JButton("Save");
				frame.setSize(500, 500);
				area.setText(text);
				area.setEditable(true);
				area.setLineWrap(true);
				area.setWrapStyleWord(true);
				frame.add(area);
				frame.add(b, BorderLayout.SOUTH);
				b.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent ev) {
						e.opts.textPane.setText(area.getText());
						e.opts.textPane.setEditable(false);
						frame.dispose();
					}
				});
				frame.setVisible(true);
				
			});
		} else if (ev.getSource() == e.opts.gui_view) {
			e.updateMiscFields = true;
		} else if (ev.getSource() == e.opts.gui_view_vec) {
			e.updateMiscFields = true;
		} else if (ev.getSource() == e.opts.gui_brush) {
			e.controls.brush_changed = true;
		} else if (ev.getSource() == e.opts.gui_adv_settings) {
			e.savemanager.writeAdvancedSettings();
			e.adv_opts.setVisible(true);
		} else if (ev.getSource() == e.adv_opts.btn_apply) {
			e.savemanager.readAdvancedSettings();
			e.adv_opts.setVisible(false);
		} else if (ev.getSource() == e.adv_opts.btn_cancel) {
			e.adv_opts.setVisible(false);
		} else if (ev.getSource() == e.opts.menu_cut) {
			cut = true;
		} else if (ev.getSource() == e.opts.menu_copy) {
			copy = true;
		} else if (ev.getSource() == e.opts.menu_paste) {
			paste = true;
		} else if (ev.getSource() == e.opts.menu_rotate) {
			rotate_selection = true;
		} else if (ev.getSource() == e.opts.menu_flip_h) {
			flip_h_selection = true;
		} else if (ev.getSource() == e.opts.menu_flip_v) {
			flip_v_selection = true;
		}else if (ev.getSource() == e.opts.menu_undo) {
			undo = true;
		} else if (ev.getSource() == e.opts.menu_redo) {
			redo = true;
		} else if (ev.getSource() == e.opts.menu_about) {
			JOptionPane.showConfirmDialog(e.opts, "Brandon's Semiconductor Simulator (SemiSim).\n (c) 2026 Brandon Li", "About", JOptionPane.OK_OPTION);
		} else if (ev.getSource() instanceof JRadioButtonMenuItem) {
			for (Brush b : brushbuttonmap.keySet())
				if (ev.getSource() == brushbuttonmap.get(b))
					e.opts.gui_brush.setSelectedItem(b);
		}
	}
	


	@Override
	public void itemStateChanged(ItemEvent ev) {
		if (ev.getSource() == e.opts.gui_brush) {
			buttongroup.setSelected(brushbuttonmap.get((Brush) e.opts.gui_brush.getSelectedItem()).getModel(), true);
		}
	}
	
	@Override
	public void mouseClicked(MouseEvent arg0) {}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent e) {
		mouse_pressed = true;
		mousebutton = e.getButton();
		mx = e.getX();
		my = e.getY();
		mx_start = e.getX();
		my_start = e.getY();
	}
	@Override
	public void mouseReleased(MouseEvent e) {
		mouse_pressed = false;
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		mx = e.getX();
		my = e.getY();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		mx = arg0.getX();
		my = arg0.getY();
	}

	private Action key_pause = new AbstractAction(null) {
		@Override
		public void actionPerformed(ActionEvent ev) {
			e.opts.gui_paused.setSelected(!e.opts.gui_paused.isSelected());
		}
	};

    private Action key_frame = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			advanceframe = true;
        }
    };

    private Action key_dbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			debugging = !debugging;
			//Timer.allEnabled = debugging;
        }
    };

    private Action key_changebrush = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush_1.setSelectedIndex((e.opts.gui_brush_1.getSelectedIndex()+1)%2);
        }
    };
    
    boolean shift_draw_override = false;

    private Action key_shift = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		shift_down = true;
    		
    		if ((Brush) e.opts.gui_brush.getSelectedItem() == Brush.DRAW) {
    			shift_draw_override = true;
    			e.opts.gui_brush.setSelectedItem(Brush.LINE);
    		}
        }
    };
    private Action key_shift_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		shift_down = false;
    		
    		if (shift_draw_override) {
    			shift_draw_override = false;
        		if ((Brush) e.opts.gui_brush.getSelectedItem() == Brush.LINE) {
        			e.opts.gui_brush.setSelectedItem(Brush.DRAW);
        		}
    		}
        }
    };
    private Action key_ctrl = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		ctrl_down = true;
        }
    };
    private Action key_ctrl_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		ctrl_down = false;
        }
    };

    private Action key_cut = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		cut = true;
        }
    };

    private Action key_copy = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		copy = true;
        }
    };

    private Action key_undo = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		undo = true;
        }
    };
    
    private Action key_redo = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		redo = true;
        }
    };

    private Action key_paste = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		paste = true;
        }
    };

    private Action key_delete = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		delete = true;
        }
    };

    private Action key_color = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_elem_colors.setSelected(!e.opts.gui_elem_colors.isSelected());
        }
    };

    ScalarView prev_scalar_view = ScalarView.NONE;
    VectorView prev_vector_view = VectorView.NONE;

    private Action key_scalar_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		if (e.opts.gui_view.getSelectedItem() == ScalarView.NONE)
    			e.opts.gui_view.setSelectedItem(prev_scalar_view);
    		else
    		{
    			prev_scalar_view = (ScalarView) e.opts.gui_view.getSelectedItem();
    			e.opts.gui_view.setSelectedItem(ScalarView.NONE);
    		}
        }
    };

    private Action key_vector_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		if (e.opts.gui_view_vec.getSelectedItem() == VectorView.NONE)
    			e.opts.gui_view_vec.setSelectedItem(prev_vector_view);
    		else
    		{
    			prev_vector_view = (VectorView) e.opts.gui_view_vec.getSelectedItem();
    			e.opts.gui_view_vec.setSelectedItem(VectorView.NONE);
    		}
        }
    };

    private Action key_tooltip = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_tooltip.setSelected(!e.opts.gui_tooltip.isSelected());
        }
    };

    private Action key_textbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_text_bg.setSelected(!e.opts.gui_text_bg.isSelected());
        }
    };

    private Action key_alt = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		alt_down = true;
        }
    };

    private Action key_alt_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			alt_down = false;
        }
    };

    private Action key_logdata = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			logdata = true;
        }
    };

    private Action key_rendertext = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			e.opts.gui_interface.setSelected(!e.opts.gui_interface.isSelected());
        }
    };
    
    private Action key_save = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		save = true;
        }
    };
    
    private Action key_open = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		load = true;
        }
    };
    
    private Action key_new = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		reset = true;
        }
    };
    
    private Action key_rotate = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		rotate_selection = true;
        }
    };
    
    private Action key_flip_v = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		flip_v_selection = true;
        }
    };
    
    private Action key_flip_h = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		flip_h_selection = true;
        }
    };
    
    private Action key_1 = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush.setSelectedItem(Brush.INTERACT);
        }
    };
    
    private Action key_2 = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush.setSelectedItem(Brush.DRAW);
        }
    };
    
    private Action key_3 = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush.setSelectedItem(Brush.LINE);
        }
    };
    
    private Action key_4 = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush.setSelectedItem(Brush.FILL);
        }
    };
    
    private Action key_5 = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.gui_brush.setSelectedItem(Brush.SELECT);
        }
    };
    
    private Action key_prevtool = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (e.opts.gui_brush.getSelectedIndex() > 0)
				e.opts.gui_brush.setSelectedIndex(e.opts.gui_brush.getSelectedIndex()-1);
        }
    };
    
    private Action key_nexttool = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (e.opts.gui_brush.getSelectedIndex() < e.opts.gui_brush.getItemCount()-1)
				e.opts.gui_brush.setSelectedIndex(e.opts.gui_brush.getSelectedIndex()+1);
        }
    };


    public void addKeyBinds(JPanel contentPane) {
    	InputMap map = contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    	InputMap map2 = contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    	
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_P, 0), key_pause);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0), key_pause);
    	contentPane.getActionMap().put(key_pause, key_pause);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_F, 0), key_frame);
    	contentPane.getActionMap().put(key_frame, key_frame);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, 0), key_dbg);
    	contentPane.getActionMap().put(key_dbg, key_dbg);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Q, 0), key_changebrush);
    	contentPane.getActionMap().put(key_changebrush, key_changebrush);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, InputEvent.SHIFT_DOWN_MASK), key_shift);
    	contentPane.getActionMap().put(key_shift, key_shift);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, 0, true), key_shift_up);
    	contentPane.getActionMap().put(key_shift_up, key_shift_up);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, InputEvent.CTRL_DOWN_MASK), key_ctrl);
    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, InputEvent.META_DOWN_MASK), key_ctrl);
    	contentPane.getActionMap().put(key_ctrl, key_ctrl);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, 0, true), key_ctrl_up);
    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, 0, true), key_ctrl_up);
    	contentPane.getActionMap().put(key_ctrl_up, key_ctrl_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_X, InputEvent.CTRL_DOWN_MASK), key_cut);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_X, InputEvent.META_DOWN_MASK), key_cut);
    	contentPane.getActionMap().put(key_cut, key_cut);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, InputEvent.CTRL_DOWN_MASK), key_copy);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, InputEvent.META_DOWN_MASK), key_copy);
    	contentPane.getActionMap().put(key_copy, key_copy);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, InputEvent.CTRL_DOWN_MASK), key_paste);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, InputEvent.META_DOWN_MASK), key_paste);
    	contentPane.getActionMap().put(key_paste, key_paste);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK), key_undo);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.META_DOWN_MASK), key_undo);
    	contentPane.getActionMap().put(key_undo, key_undo);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK), key_redo);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.META_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK), key_redo);
    	contentPane.getActionMap().put(key_redo, key_redo);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK), key_save);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.META_DOWN_MASK), key_save);
    	contentPane.getActionMap().put(key_save, key_save);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_DOWN_MASK), key_open);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.META_DOWN_MASK), key_open);
    	contentPane.getActionMap().put(key_open, key_open);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.CTRL_DOWN_MASK), key_new);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.META_DOWN_MASK), key_new);
    	contentPane.getActionMap().put(key_new, key_new);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_R, InputEvent.CTRL_DOWN_MASK), key_rotate);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_R, InputEvent.META_DOWN_MASK), key_rotate);
    	contentPane.getActionMap().put(key_rotate, key_rotate);
    	
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_G, InputEvent.CTRL_DOWN_MASK), key_flip_v);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_G, InputEvent.META_DOWN_MASK), key_flip_v);
    	contentPane.getActionMap().put(key_flip_v, key_flip_v);
    	
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_F, InputEvent.CTRL_DOWN_MASK), key_flip_h);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_F, InputEvent.META_DOWN_MASK), key_flip_h);
    	contentPane.getActionMap().put(key_flip_h, key_flip_h);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE, 0), key_delete);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0), key_delete);
    	contentPane.getActionMap().put(key_delete, key_delete);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, 0), key_color);
    	contentPane.getActionMap().put(key_color, key_color);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, 0), key_scalar_view);
    	contentPane.getActionMap().put(key_scalar_view, key_scalar_view);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, 0), key_vector_view);
    	contentPane.getActionMap().put(key_vector_view, key_vector_view);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_T, 0), key_tooltip);
    	contentPane.getActionMap().put(key_tooltip, key_tooltip);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_G, 0), key_textbg);
    	contentPane.getActionMap().put(key_textbg, key_textbg);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, InputEvent.ALT_DOWN_MASK), key_alt);
    	contentPane.getActionMap().put(key_alt, key_alt);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, 0, true), key_alt_up);
    	contentPane.getActionMap().put(key_alt_up, key_alt_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_R, 0), key_logdata);
    	contentPane.getActionMap().put(key_logdata, key_logdata);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_H, 0), key_rendertext);
    	contentPane.getActionMap().put(key_rendertext, key_rendertext);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_1, 0), key_1);
    	contentPane.getActionMap().put(key_1, key_1);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_2, 0), key_2);
    	contentPane.getActionMap().put(key_2, key_2);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_3, 0), key_3);
    	contentPane.getActionMap().put(key_3, key_3);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_4, 0), key_4);
    	contentPane.getActionMap().put(key_4, key_4);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_5, 0), key_5);
    	contentPane.getActionMap().put(key_5, key_5);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_OPEN_BRACKET, 0), key_prevtool);
    	contentPane.getActionMap().put(key_prevtool, key_prevtool);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_CLOSE_BRACKET, 0), key_nexttool);
    	contentPane.getActionMap().put(key_nexttool, key_nexttool);

    }
    

    public void removeKeyBinds(JPanel contentPane) {
    	InputMap map = contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    	InputMap map2 = contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
    	
    	map.clear();
    	map2.clear();
    }

	@Override
	public void mouseWheelMoved(MouseWheelEvent ev) {
		e.opts.gui_brushsize.setValue(e.opts.gui_brushsize.getValue() - (int)(10*ev.getPreciseWheelRotation()));
		System.out.println(ev.getPreciseWheelRotation());
	}
	

	@Override
	public void keyTyped(KeyEvent ev) {
		if (texting) {
			if (Font7x5.getCharacter(ev.getKeyChar()) != null) {
				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 7; j++) {
						if (text_x+i+1 >= 0 && text_x+i+1 < e.nx && text_y+j >= 0 && text_y+j < e.ny
								&& Font7x5.getPixel(ev.getKeyChar(), 4-i, j) == 1 && e.materials[text_x+i+1][text_y+j].type == MaterialType.VACUUM) {
							e.initializeMaterial(text_x+i+1, text_y+j, MaterialType.DECO);
						}
					}
				}
				text_x += 6;
			}
			e.updateAllMaterials(false);
		}
	}

	@Override
	public void keyPressed(KeyEvent ev) {
		if (texting) {
			if (ev.getKeyCode() == KeyEvent.VK_BACK_SPACE || ev.getKeyCode() == KeyEvent.VK_DELETE) {
				text_x -= 6;
				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 7; j++) {
						if (text_x+i+1 >= 0 && text_x+i+1 < e.nx && text_y+j >= 0 && text_y+j < e.ny
								&& e.materials[text_x+i+1][text_y+j].type == MaterialType.DECO) {
							e.eraseMaterial(text_x+i+1,text_y+j);
						}
					}
				}
			}
			e.updateAllMaterials(false);
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {}

	public double bilinearinterp(double[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-2 && Math.abs(y-Math.round(y)) < 1e-2) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= e.nx) i = e.nx - 1;
			if (j >= e.ny) j = e.ny - 1;
			return array[i][j];
		}

		if (xfloor < 0) {
			xfloor = 0;
			fx = 0.0;
		} else if (xfloor >= e.nx - 1) {
			xfloor = e.nx - 2;
			fx = 1.0;
		}
		if (yfloor < 0) {
			yfloor = 0;
			fy = 0.0;
		} else if (yfloor >= e.ny - 1) {
			yfloor = e.ny - 2;
			fy = 1.0;
		}
		double va = array[xfloor][yfloor]*(1.0-fx) + array[xfloor+1][yfloor]*fx;
		double vb = array[xfloor][yfloor+1]*(1.0-fx) + array[xfloor+1][yfloor+1]*fx;

		return va*(1.0-fy) + vb*fy;
	}
	
	class FloodFillCoordinate {
		int i;
		int j;

		public FloodFillCoordinate(int i, int j) {
			this.i = i;
			this.j = j;
		}
	}

	public enum BrushShape {
		CIRCLE("Circle brush"),
		SQUARE("Square brush");
	
		public String name;
		BrushShape(String name)
		{
			this.name = name;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}

	public enum Brush {
		INTERACT("Interact"),
		LIGHT("Flashlight"),
		DRAW("Draw"),
		REPLACE("Replace"),
		LINE("Line"),
		FILL("Fill"),
		ERASE("Eraser"),
		SELECT("Select and move"),
		FLOODSELECT("Flood select"),
		TEXT("Text"),
		VOLTAGE("Voltage probe"),
		CURRENT("Current probe"),
		CHARGE("Charge probe"),
		GROUND("Ground"),
		DELETEPROBE("Delete probe"),
		BANDS("Plot bands"),
		SCALARPLOT("Plot scalar field"),
		CARRIERPLOT("Plot carriers");
	
		public String name;
		Brush(String name)
		{
			this.name = name;
		}
	
		@Override
		public String toString() {
			return name;
		}
	
		public static boolean isMaterialModifyingBrush(Brush brush) {
			return (brush == Brush.DRAW
					|| brush == Brush.LINE
					|| brush == Brush.REPLACE
					|| brush == Brush.ERASE
					|| brush == Brush.FILL);
		}
	
		public static boolean isBrushShapeImportant(Brush brush) {
			return (brush == Brush.DRAW
					|| brush == Brush.LINE
					|| brush == Brush.REPLACE
					|| brush == Brush.ERASE
					|| brush == Brush.LIGHT);
		}
		
		public static boolean drawLine(Brush brush) {
			return (brush == Brush.LINE || brush == Brush.BANDS || brush == Brush.SCALARPLOT || brush == Brush.CARRIERPLOT);
		}
	}
}

interface FloodFillFunc {
	boolean isValid(int i, int j);
	void fill(int i, int j);
}
