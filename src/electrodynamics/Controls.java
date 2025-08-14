// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.Cursor;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.LinkedList;
import java.util.Queue;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.KeyStroke;

import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorView;
import electrodynamics.Simulation.BoundaryCondition;
import electrodynamics.plot.Plot;
import electrodynamics.util.Font7x5;
import electrodynamics.util.Timer;
import electrodynamics.util.Utils;
import electrodynamics.util.Vector;

public class Controls implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener {
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
	public boolean delete = false;
    public boolean shift_down = false;
    public boolean ctrl_down = false;
    public boolean alt_down = false;
    public boolean logdata = false;


	/* Mouse controls */

	PointerInfo pointerinfo = MouseInfo.getPointerInfo();
	public boolean mouse_pressed = false;
	public boolean mouse_pressed_prev = false;
	public boolean modifier_pressed = false;
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
		texting = false;
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


		if (!mouse_pressed && !releasing) {

			if (ctrl_down || shift_down) {
				if (!modifier_pressed && Brush.isMaterialModifyingBrush(brush)) {
					modifier_pressed = true;
					prev_brush = brush;

					if (ctrl_down) {
						e.opts.gui_brush.setSelectedItem(Brush.FILL);
						brush = (Brush) e.opts.gui_brush.getSelectedItem();
					}
					else if (shift_down) {
						e.opts.gui_brush.setSelectedItem(Brush.LINE);
						brush = (Brush) e.opts.gui_brush.getSelectedItem();
					}
					//r.requestFocus();
				}
			} else {
				if (modifier_pressed) {
					modifier_pressed = false;
					e.opts.gui_brush.setSelectedItem(Brush.DRAW);
					brush = (Brush) e.opts.gui_brush.getSelectedItem();
					//r.requestFocus();
				}
			}
		}

		mx_realspace = Math.round((mx-1)/(double)e.renderer.scalefactor - 0.5)*e.ds;
		my_realspace = Math.round((my-1)/(double)e.renderer.scalefactor - 0.5)*e.ds;

		mx_start_realspace =  Math.round((mx_start-1)/(double)e.renderer.scalefactor - 0.5)*e.ds;
		my_start_realspace = Math.round((my_start-1)/(double)e.renderer.scalefactor - 0.5)*e.ds;

		mx_index = (int)Math.round((mx-1)/(double)e.renderer.scalefactor - 0.5);
		my_index = (int)Math.round((my-1)/(double)e.renderer.scalefactor - 0.5);

		mx_start_index = (int)Math.round((mx_start-1)/(double)e.renderer.scalefactor - 0.5);
		my_start_index = (int)Math.round((my_start-1)/(double)e.renderer.scalefactor - 0.5);

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

		brushsize = ((e.ds*e.nx)/500)*(Math.pow(10.0, e.opts.gui_brushsize.getValue()/500.0) + e.opts.gui_brushsize.getValue()/100.0);
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

			if (mat == MaterialType.EMF) {
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
						floodFillSet(mx_index, my_index, e.materials[mx_index][my_index].type, mat, angle);
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
						if (e.opts.gui_brush_highlight.isSelected()) {
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
						else {
							under_brush[i][j] = false;
						}
					}
				}
			}

			break;

		case INTERACT:
			if (e.materials[mx_index][my_index].type == MaterialType.EMF || e.materials[mx_index][my_index].type == MaterialType.SWITCH)
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

				if (e.materials[mx_index][my_index].type == MaterialType.EMF && turn_on_EMF) {
					floodFillSelectEMF(mx_index, my_index, true);
					EMF_selected = true;
					int setting = (int)(Math.round(50*e.materials[mx_index][my_index].emf/max_EMF));
					e.opts.gui_parameter3.setValue(setting);
					prev_EMF_setting = setting;
				}

				if (e.materials[mx_index][my_index].type == MaterialType.SWITCH) {
					this.floodFillToggleSwitch(mx_index, my_index, 1-e.materials[mx_index][my_index].activated);
					update = true;
				}
			}
			break;
		case FLOODSELECT:
		case SELECT:
			if (pressing) {
				if (brush == Brush.FLOODSELECT && !moving_selection) {
					floodFillSelect(mx_index, my_index, e.materials[mx_index][my_index].type, !selected[mx_index][my_index]);
				}
				else if (moving_selection && !dragging_selection) {
					for (int i = 0; i < e.nx; i++)
					{
						for (int j = 0; j < e.ny; j++)
						{
							int si = i-delta_mx_index;
							int sj = j-delta_my_index;
							if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && selection[si][sj].type != MaterialType.VACUUM) {
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
								if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && selection[si][sj].type != MaterialType.VACUUM) {
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
				texting = true;
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
			texting = false;
		}
		
		for (Plot p : e.plots) {
			p.updatePlot(e);
		}

		setEMFs();

		if ((releasing && Brush.isMaterialModifyingBrush(brush)) || (BoundaryCondition)e.opts.gui_bc.getSelectedItem() != prev_boundary || update) {

			//e.resetFields(false);
			//e.constructBoundary();
			e.updateAllMaterials(false);
			e.multigridSolve(true, false);
		}

		prev_boundary = (BoundaryCondition)e.opts.gui_bc.getSelectedItem();

		mxp_realspace = mx_realspace;
		myp_realspace = my_realspace;
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
						if (mat == MaterialType.EMF) e.materials[i][j].emf_direction = EMF_angle;
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
					if (e.materials[i][j].type == MaterialType.EMF && selected_EMF[i][j]) {
						e.materials[i][j].emf = new_EMF;
					}
				}
			}

			e.updateJustEMFs();
		}

		prev_EMF_setting = EMF_setting;
	}

	public void floodFillSet(int i, int j, MaterialType old_mat, MaterialType new_mat, double EMF_angle) {
		if (old_mat == new_mat)
			return;

		Queue<FloodFillCoordinate> queue = new LinkedList<>();
		queue.add(new FloodFillCoordinate(i, j));

		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < e.nx && coord.j >= 0 && coord.j < e.ny && e.materials[coord.i][coord.j].type == old_mat && e.materials[coord.i][coord.j].type != new_mat) {
				e.eraseMaterial(coord.i, coord.j);
				e.initializeMaterial(coord.i, coord.j, new_mat);
				if (new_mat == MaterialType.EMF) e.materials[coord.i][coord.j].emf_direction = EMF_angle;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}

	public void floodFillSelect(int i, int j, MaterialType mat, boolean select) {
		Queue<FloodFillCoordinate> queue = new LinkedList<>();
		queue.add(new FloodFillCoordinate(i, j));

		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < e.nx && coord.j >= 0 && coord.j < e.ny && e.materials[coord.i][coord.j].type == mat && selected[coord.i][coord.j] != select) {
				selected[coord.i][coord.j] = select;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}

	public void floodFillSelectEMF(int i, int j, boolean select) {
		Queue<FloodFillCoordinate> queue = new LinkedList<>();
		queue.add(new FloodFillCoordinate(i, j));

		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < e.nx && coord.j >= 0 && coord.j < e.ny && e.materials[coord.i][coord.j].type == MaterialType.EMF && selected_EMF[coord.i][coord.j] != select) {
				selected_EMF[coord.i][coord.j] = select;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}

	public void floodFillToggleSwitch(int i, int j, int active) {
		Queue<FloodFillCoordinate> queue = new LinkedList<>();
		queue.add(new FloodFillCoordinate(i, j));

		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < e.nx && coord.j >= 0 && coord.j < e.ny && e.materials[coord.i][coord.j].type == MaterialType.SWITCH && e.materials[coord.i][coord.j].activated != active) {
				e.materials[coord.i][coord.j].activated = active;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
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

	@SuppressWarnings("serial")
	private Action key_pause = new AbstractAction(null) {
		@Override
		public void actionPerformed(ActionEvent ev) {
			if (texting) return;
			e.opts.gui_paused.setSelected(!e.opts.gui_paused.isSelected());
		}
	};

    @SuppressWarnings("serial")
    private Action key_frame = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
			advanceframe = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_dbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
			debugging = !debugging;
			Timer.allEnabled = debugging;
        }
    };

    @SuppressWarnings("serial")
    private Action key_changebrush = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		e.opts.gui_brush_1.setSelectedIndex((e.opts.gui_brush_1.getSelectedIndex()+1)%2);
        }
    };

    @SuppressWarnings("serial")
    private Action key_shift = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		shift_down = true;
        }
    };
    @SuppressWarnings("serial")
    private Action key_shift_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		shift_down = false;
        }
    };
    @SuppressWarnings("serial")
    private Action key_ctrl = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		ctrl_down = true;
        }
    };
    @SuppressWarnings("serial")
    private Action key_ctrl_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		ctrl_down = false;
        }
    };

    @SuppressWarnings("serial")
    private Action key_cut = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		cut = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_copy = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		copy = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_paste = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		paste = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_delete = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		delete = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_color = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		e.opts.gui_elem_colors.setSelected(!e.opts.gui_elem_colors.isSelected());
        }
    };

    ScalarView prev_scalar_view = ScalarView.NONE;
    VectorView prev_vector_view = VectorView.NONE;

    @SuppressWarnings("serial")
    private Action key_scalar_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		if (e.opts.gui_view.getSelectedItem() == ScalarView.NONE)
    			e.opts.gui_view.setSelectedItem(prev_scalar_view);
    		else
    		{
    			prev_scalar_view = (ScalarView) e.opts.gui_view.getSelectedItem();
    			e.opts.gui_view.setSelectedItem(ScalarView.NONE);
    		}
        }
    };

    @SuppressWarnings("serial")
    private Action key_vector_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		if (e.opts.gui_view_vec.getSelectedItem() == VectorView.NONE)
    			e.opts.gui_view_vec.setSelectedItem(prev_vector_view);
    		else
    		{
    			prev_vector_view = (VectorView) e.opts.gui_view_vec.getSelectedItem();
    			e.opts.gui_view_vec.setSelectedItem(VectorView.NONE);
    		}
        }
    };

    @SuppressWarnings("serial")
    private Action key_tooltip = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		e.opts.gui_tooltip.setSelected(!e.opts.gui_tooltip.isSelected());
        }
    };

    @SuppressWarnings("serial")
    private Action key_textbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			if (texting) return;
    		e.opts.gui_text_bg.setSelected(!e.opts.gui_text_bg.isSelected());
        }
    };

    @SuppressWarnings("serial")
    private Action key_alt = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		alt_down = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_alt_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			alt_down = false;
        }
    };

    @SuppressWarnings("serial")
    private Action key_logdata = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			logdata = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_rendertext = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
			e.opts.gui_interface.setSelected(!e.opts.gui_interface.isSelected());
        }
    };

    public void addKeyBinds(JPanel contentPane) {
    	InputMap map = contentPane.getInputMap(JComponent.WHEN_FOCUSED);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_P, 0), key_pause);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0), key_pause);
    	contentPane.getActionMap().put(key_pause, key_pause);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_F, 0), key_frame);
    	contentPane.getActionMap().put(key_frame, key_frame);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, 0), key_dbg);
    	contentPane.getActionMap().put(key_dbg, key_dbg);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Q, 0), key_changebrush);
    	contentPane.getActionMap().put(key_changebrush, key_changebrush);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, InputEvent.SHIFT_DOWN_MASK), key_shift);
    	contentPane.getActionMap().put(key_shift, key_shift);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, 0, true), key_shift_up);
    	contentPane.getActionMap().put(key_shift_up, key_shift_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, InputEvent.CTRL_DOWN_MASK), key_ctrl);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, InputEvent.META_DOWN_MASK), key_ctrl);
    	contentPane.getActionMap().put(key_ctrl, key_ctrl);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, 0, true), key_ctrl_up);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, 0, true), key_ctrl_up);
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

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, InputEvent.ALT_DOWN_MASK), key_alt);
    	contentPane.getActionMap().put(key_alt, key_alt);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, 0, true), key_alt_up);
    	contentPane.getActionMap().put(key_alt_up, key_alt_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_R, 0), key_logdata);
    	contentPane.getActionMap().put(key_logdata, key_logdata);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_H, 0), key_rendertext);
    	contentPane.getActionMap().put(key_rendertext, key_rendertext);

    }

	@Override
	public void mouseWheelMoved(MouseWheelEvent ev) {
		e.opts.gui_brushsize.setValue(e.opts.gui_brushsize.getValue() - (int)(10*ev.getPreciseWheelRotation()));
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
		DRAW("Draw"),
		VOLTAGE("Voltage probe"),
		CURRENT("Current probe"),
		GROUND("Ground"),
		DELETEPROBE("Delete probe"),
		BANDS("Plot bands"),
		SCALARPLOT("Plot scalar field"),
		CARRIERPLOT("Plot carriers"),
		LIGHT("Flashlight"),
		REPLACE("Replace"),
		LINE("Line"),
		FILL("Fill"),
		ERASE("Eraser"),
		SELECT("Select and move"),
		FLOODSELECT("Flood select"),
		TEXT("Text");
	
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
