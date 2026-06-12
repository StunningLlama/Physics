// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.FlowLayout;
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
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.Queue;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;

import electrodynamics.Renderer.ScalarMode;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorMode;
import electrodynamics.Renderer.VectorView;
import electrodynamics.Simulation.BoundaryCondition;
import electrodynamics.plot.LinePath;
import electrodynamics.plot.Path;
import electrodynamics.plot.Plot;
import electrodynamics.plot.ProbePlot;
import electrodynamics.plot.SegmentedPath;
import electrodynamics.probe.AreaProbe;
import electrodynamics.probe.AreaProbe.QuantityType;
import electrodynamics.probe.ChargeProbe;
import electrodynamics.probe.CurrentProbe;
import electrodynamics.probe.FluxProbe;
import electrodynamics.probe.Ground;
import electrodynamics.probe.LineProbe;
import electrodynamics.probe.PointProbe;
import electrodynamics.probe.Probe;
import electrodynamics.probe.Ruler;
import electrodynamics.probe.VoltageProbe;
import electrodynamics.units.Quantity;
import electrodynamics.util.CustJRadioButtonMenuItem;
import electrodynamics.util.Font7x5;
import electrodynamics.util.MenuCheckList;
import electrodynamics.util.Utils;
import electrodynamics.util.Vector;

public class Controls implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener, ItemListener, WindowListener {
	Simulation e;
	
	/* Keyboard controls */

	public boolean advanceframe = false;
	public boolean clear = false;
	public boolean reset = false;
	public boolean save = false;
	public boolean saveas = false;
	public boolean load = false;
	public boolean debugging = false;
	public boolean cut = false;
	public boolean copy = false;
	public boolean paste = false;
	public boolean undo = false;
	public boolean redo = false;
	public boolean selectall = false;
	public boolean deselectall = false;
	public boolean delete = false;
    public boolean shift_down = false;
    public boolean ctrl_down = false;
    public boolean alt_down = false;
    public boolean logdata = false;
    public boolean rotate_selection = false;
    public boolean flip_h_selection = false;
    public boolean flip_v_selection = false;
    public boolean exit = false;
    public boolean updateimagesize = false;
    public boolean test = false;


	/* Mouse controls */

	PointerInfo pointerinfo = MouseInfo.getPointerInfo();
	public boolean mouse_pressed = false;
	public boolean mouse_pressed_prev = false;
	public boolean pressing = false;
	public boolean releasing = false;
	public boolean moving_selection = false;
	public boolean dragging_selection = false;
	public boolean brush_changed = false;
	public boolean changesmade = false;

	public int mousebutton = 0;
	public int mx_screen = 0;
	public int my_screen = 0;
	public int mx_start_screen = 0;
	public int my_start_screen = 0;

	public int mx = 0;
	public int my = 0;
	public int mx_start = 0;
	public int my_start = 0;
	public int mxp = 0;
	public int myp = 0;

	public int delta_mx = 0;
	public int delta_my = 0;
	
	public int zoom_i1 = 0;
	public int zoom_j1 = 0;
	public int zoom_i2 = 0;
	public int zoom_j2 = 0;
	
	public int zoom_i1_pan = 0;
	public int zoom_j1_pan = 0;
	public int zoom_i2_pan = 0;
	public int zoom_j2_pan = 0;
	public boolean zoomed = false;

	public boolean EMF_selected = false;

	Brush prev_brush;
	public double brushsize = 0;
	public int prev_EMF_setting = 0;
	BoundaryCondition prev_boundary = BoundaryCondition.DISSIPATIVE;

	public boolean[][] under_brush;
	public boolean[][] selected;
	public boolean[][] selected_EMF;
	
	public boolean iscurrentselected = false;

	public Clipboard selection;
	public Clipboard clipboard;
	public boolean selectionempty = true;
	public boolean clipboardempty = true;

	public int text_x = 0;
	public int text_y = 0;
	public boolean texting = false;

	public double flashlight_strength;
	
	public int plotinterval = 10;

	Cursor HAND_CURSOR = new Cursor(Cursor.HAND_CURSOR);
	Cursor DEFAULT_CURSOR = new Cursor(Cursor.DEFAULT_CURSOR);
	
	public MenuCheckList<Brush, JRadioButtonMenuItem> brushes = new MenuCheckList<Brush, JRadioButtonMenuItem>();
	public MenuCheckList<ScalarView, CustJRadioButtonMenuItem> scalarview = new MenuCheckList<ScalarView, CustJRadioButtonMenuItem>();
	public MenuCheckList<VectorView, CustJRadioButtonMenuItem> vectorview = new MenuCheckList<VectorView, CustJRadioButtonMenuItem>();
	public MenuCheckList<ScalarMode, CustJRadioButtonMenuItem> scalarmode = new MenuCheckList<ScalarMode, CustJRadioButtonMenuItem>();
	public MenuCheckList<VectorMode, CustJRadioButtonMenuItem> vectormode = new MenuCheckList<VectorMode, CustJRadioButtonMenuItem>();
	//public JCheckBoxMenuItem carriers = new JCheckBoxMenuItem("Show charge carriers");
	
	public UndoRedo undoredo = new UndoRedo(4);
	
	public Probe.LabelCoord labelcoord = null;
	
	public Path plotpath = null;
	
	public Controls(Simulation e) {
		this.e = e;
	}
	
	void setResolution(int resolution) {
		selection = new Clipboard(e);
		clipboard = new Clipboard(e);

		under_brush = new boolean[e.nx][e.ny];
		selected = new boolean[e.nx][e.ny];
		selected_EMF = new boolean[e.nx][e.ny];
		
		text_x = 0;
		text_y = 0;
		endTextInput();
	}

	public void handleMouseInput() {
		transformMouseCoords();
		processKeyboardCommands();
		makeUIchanges();
		applyTool();
	}
	
	public void transformMouseCoords() {
		pressing = false;
		releasing = false;

		if (mouse_pressed) {
			if (!mouse_pressed_prev) {
				pressing = true;
			}
		} else {
			if (mouse_pressed_prev) {
				releasing = true;
			}
		}
		mouse_pressed_prev = mouse_pressed;


		double sf_x = (zoom_i2-zoom_i1+1)/(double)e.canvas.zoom_bound_x;
		double sf_y = (zoom_j2-zoom_j1+1)/(double)e.canvas.zoom_bound_y;

		mx = (int)Math.round(zoom_i1 + (mx_screen-e.canvas.offset_x - 1)*sf_x - 0.5);
		my = (int)Math.round(zoom_j1 + (my_screen-e.canvas.offset_y - 2)*sf_y - 0.5);

		mx_start = (int)Math.round(zoom_i1 + (mx_start_screen-e.canvas.offset_x - 1)*sf_x - 0.5);
		my_start = (int)Math.round(zoom_j1 + (my_start_screen-e.canvas.offset_y - 2)*sf_y - 0.5);
		
		if (alt_down && (mouse_pressed || releasing))
			snapToCardinals(mx, my);

		if (mx < 0) mx = 0;
		if (my < 0) my = 0;
		if (mx >= e.nx) mx = e.nx-1;
		if (my >= e.ny) my = e.ny-1;

		if (mx_start < 0) mx_start = 0;
		if (my_start < 0) my_start = 0;
		if (mx_start >= e.nx) mx_start = e.nx-1;
		if (my_start >= e.ny) my_start = e.ny-1;
	}

	
	public void processKeyboardCommands() {
		Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();
		
		if (brush_changed) {
			if (!(brush == Brush.SELECT || brush == Brush.FLOODSELECT)) {
				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						selected[i][j] = false;
					}
				}
				
				e.probes.forEach((p) -> p.selected = false);
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

		if (cut || copy) {
			clipboardempty = clipboard.copy(cut);

			if (cut) flagChanges(true);

			cut = false;
			copy = false;
		}

		if (paste) {
			clipboard.transfer(selection);
			
			moving_selection = true;
			dragging_selection = false;
			paste = false;
			flagChanges(true);
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
			e.probes.removeIf((p) -> p.selected);
			e.relabelProbes();
			
			delete = false;
			flagChanges(true);
		}
		
		if (rotate_selection || flip_h_selection || flip_v_selection) {

			if (rotate_selection) {
				selection.rotate90();
				rotate_selection = false;
			}
			
			if (flip_h_selection) {
				selection.flip_h();
				flip_h_selection = false;
			}
			
			if (flip_v_selection) {
				selection.flip_v();
				flip_v_selection = false;
			}
		}
		
		if (selectall) {
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					if (!(e.materials[i][j].auto_placed && e.materials[i][j].type == MaterialType.ABSORBER)) {
						selected[i][j] = true;
					}
				}
			}
			e.probes.forEach((p) -> p.selected = true);
			e.opts.gui_brush.setSelectedItem(Brush.SELECT);
			selectall = false;
		}
		
		if (deselectall) {
			for (int i = 0; i < e.nx; i++)
			{
				for (int j = 0; j < e.ny; j++)
				{
					selected[i][j] = false;
				}
			}
			e.probes.forEach((p) -> p.selected = false);
			deselectall = false;
		}
		
		selectionempty = true;

		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (selected[i][j]) {
					selectionempty = false;
				}
			}
		}
		for (Probe p : e.probes)
			if (p.selected)
				selectionempty = false;
		
		e.opts.menu_cut.setEnabled(!selectionempty);
		e.opts.menu_copy.setEnabled(!selectionempty);
		e.opts.menu_paste.setEnabled(!clipboardempty);
		e.opts.menu_rotate.setEnabled(moving_selection);
		e.opts.menu_flip_h.setEnabled(moving_selection);
		e.opts.menu_flip_v.setEnabled(moving_selection);
		e.opts.menu_undo.setEnabled(undoredo.canUndo());
		e.opts.menu_redo.setEnabled(undoredo.canRedo());
		e.opts.menu_deselectall.setEnabled(!selectionempty);
	}
	
	public void makeUIchanges() {
		Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();
		
		if (pressing) {
			e.canvas.requestFocus();
			if (Brush.isMaterialModifyingBrush(brush) || brush == Brush.SELECT)
				e.opts.gui_paused.setSelected(true);
		}
		
		e.opts.gui_stepsizelbl.setText("Timestep: " + e.units.toString(e.dt, Quantity.TIME));
		e.opts.gui_stepslbl.setText("Sim steps/frame: " + e.opts.gui_simspeed_2.getValue());

		brushsize = Math.pow(10.0, 2*e.opts.gui_brushsize.getValue()/(50.0*10.0) - 0.75) + e.opts.gui_brushsize.getValue()/10.0 + 0.5;
		e.opts.lblBrushSize.setText("Brush size: " + (int)Math.ceil(brushsize));

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
		
		if (e.opts.gui_carriers.isSelected()) {
			e.opts.gui_carrierlbl.setVisible(true);
			e.opts.gui_carrier_density.setVisible(true);
		} else {
			e.opts.gui_carrierlbl.setVisible(false);
			e.opts.gui_carrier_density.setVisible(false);
		}
		
		if (brush == Brush.PROBEPLOT) {
			e.opts.gui_plotinterval.setVisible(true);
			e.opts.gui_plotinterval_text.setVisible(true);
		} else {
			e.opts.gui_plotinterval.setVisible(false);
			e.opts.gui_plotinterval_text.setVisible(false);
		}

		if (brush == Brush.LIGHT) {
			e.opts.gui_light.setVisible(true);
			e.opts.gui_light_text.setVisible(true);
		} else {
			e.opts.gui_light.setVisible(false);
			e.opts.gui_light_text.setVisible(false);
		}

		if (brush == Brush.PROBE) {
			e.opts.gui_probetype.setVisible(true);
		} else {
			e.opts.gui_probetype.setVisible(false);
		}
	}
	
	public void applyTool() {
		
		Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();
		BrushShape brushshape = (BrushShape) e.opts.gui_brush_1.getSelectedItem();
		
		switch(brush) {
		case DRAW:
		case LINE:
		case REPLACE:
		case ERASE:
		case FILL:
		case LIGHT:

			boolean pick_material = alt_down && mx-mx_start == 0 && my-my_start == 0;
			
			if ((mousebutton == MouseEvent.BUTTON2 || pick_material) && releasing) {
				e.opts.gui_material.setSelectedItem(new GeneralMaterialType(e.materials[mx][my]));
			}

			double angle = 0;
			GeneralMaterialType mat = (GeneralMaterialType) e.opts.gui_material.getSelectedItem();

			if (mat.type.hasEMF() && brush != Brush.LIGHT) {
				e.opts.gui_parameter2.setVisible(true);
				e.opts.gui_parameter2_text.setVisible(true);

				int directionval = e.opts.gui_parameter2.getValue()/6;
				if (directionval == 0) {
					e.opts.gui_parameter2_text.setText("Direction: Up");
					angle = -Math.PI/2;
				}
				if (directionval == 1) {
					e.opts.gui_parameter2_text.setText("Direction: Right");
					angle = 0;
				}
				if (directionval == 2) {
					e.opts.gui_parameter2_text.setText("Direction: Down");
					angle = Math.PI/2;
				}
				if (directionval == 3) {
					e.opts.gui_parameter2_text.setText("Direction: Left");
					angle = Math.PI;
				}
				//e.opts.gui_parameter2_text.setText("Brush orientation: " + directionval*(360/24) + " deg");
				//angle = Math.PI * directionval/12.0;
			} else {
				e.opts.gui_parameter2.setVisible(false);
				e.opts.gui_parameter2_text.setVisible(false);
				e.opts.gui_parameter2_text.setText("");
			}
			
			if (brush == Brush.LIGHT) {
				flashlight_strength = e.default_flashlight_strength*Math.pow(10, e.opts.gui_light.getValue()/10.0);
				e.opts.gui_light_text.setText("Light: " + e.units.toString(flashlight_strength*e.Eg_semi*e.depth, Quantity.INTENSITY));
			}


			if (mousebutton == MouseEvent.BUTTON3 || brush == Brush.ERASE)
				mat = GeneralMaterialType.EMPTY;

			if (!(mousebutton == MouseEvent.BUTTON2 || pick_material)) {
				if (brush == Brush.LINE) {
					if (releasing) {
						drawMaterialLine(mx_start, my_start, mx, my, brush, brushshape, mat, brushsize, angle);
					}
				} else if (brush == Brush.FILL) {
					if (pressing) {
						GeneralMaterialType old_mat = new GeneralMaterialType(e.materials[mx][my]);
						GeneralMaterialType new_mat = mat;
						double new_angle = angle;
						if (!new_mat.equals(old_mat)) {
							this.floodFill(mx, my, new FloodFillFunc() {
								@Override
								public boolean isValid(int i, int j) {
									return e.materials[i][j].type == old_mat.type && e.materials[i][j].cust_id == old_mat.cust_id;
								}

								@Override
								public void fill(int i, int j) {
									e.eraseMaterial(i, j);
									e.initializeMaterial(e.materials[i][j], new_mat);
									if (new_mat.type.hasEMF()) e.materials[i][j].emf_direction = new_angle;
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
								cx = i;
								cy = j;

								double px = (cx-mx);
								double py = (cy-my);
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
					drawMaterialLine(mxp, myp, mx, my, brush, brushshape, mat, brushsize, angle);
				}
			}

			if (Brush.isBrushShapeImportant(brush)) {
				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						double cx = 0;
						double cy = 0;
						cx = i;
						cy = j;

						double px = cx-mx;
						double py = cy-my;
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
			if (e.materials[mx][my].type.isInteractable())
				e.canvas.setCursor(HAND_CURSOR);
			else
				e.canvas.setCursor(DEFAULT_CURSOR);

			if (pressing) {
				boolean turn_on_EMF = !selected_EMF[mx][my];

				for (int i = 0; i < e.nx; i++)
				{
					for (int j = 0; j < e.ny; j++)
					{
						selected_EMF[i][j] = false;
					}
				}
				EMF_selected = false;

				if ((e.materials[mx][my].type.hasEMF()) && turn_on_EMF) {
					this.floodFill(mx, my, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == e.materials[mx][my].type && selected_EMF[i][j] != true;
						}

						@Override
						public void fill(int i, int j) {
							selected_EMF[i][j] = true;
						}
					});
					
					EMF_selected = true;
					iscurrentselected = e.materials[mx][my].type == MaterialType.CURRENT;
					
					int setting = 0;
					if (!iscurrentselected)
						setting = (int)(Math.round(50*e.materials[mx][my].emf/e.max_EMF));
					else {
						setting = (int)(Math.round(50*e.materials[mx][my].emf/(e.max_current/e.currentsource_sigma)));
					}
					
					e.opts.gui_parameter3.setValue(setting);
					prev_EMF_setting = setting;
				}

				if (e.materials[mx][my].type == MaterialType.SWITCH) {
					int active = 1-e.materials[mx][my].activated;
					
					this.floodFill(mx, my, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == MaterialType.SWITCH && e.materials[i][j].activated != active;
						}

						@Override
						public void fill(int i, int j) {
							e.materials[i][j].activated = active;
						}
					});
				
					updatematerials = true;
				}
			}
			break;
		case ZOOM:
			if (pressing && shift_down) {
				zoom_i1_pan = zoom_i1;
				zoom_j1_pan = zoom_j1;
				zoom_i2_pan = zoom_i2;
				zoom_j2_pan = zoom_j2;
			} else if (mouse_pressed && shift_down) {
				double sf_x = (zoom_i2_pan-zoom_i1_pan+1)/(double)e.canvas.zoom_bound_x;
				double sf_y = (zoom_j2_pan-zoom_j1_pan+1)/(double)e.canvas.zoom_bound_y;

				int mx_tmp = (int)Math.round(zoom_i1_pan + (mx_screen-e.canvas.offset_x - 1)*sf_x - 0.5);
				int my_tmp = (int)Math.round(zoom_j1_pan + (my_screen-e.canvas.offset_y - 2)*sf_y - 0.5);

				int mx_start_tmp = (int)Math.round(zoom_i1_pan + (mx_start_screen-e.canvas.offset_x - 1)*sf_x - 0.5);
				int my_start_tmp = (int)Math.round(zoom_j1_pan + (my_start_screen-e.canvas.offset_y - 2)*sf_y - 0.5);
				
				zoom_i1 = zoom_i1_pan - (mx_tmp - mx_start_tmp);
				zoom_j1 = zoom_j1_pan - (my_tmp - my_start_tmp);
				zoom_i2 = zoom_i2_pan - (mx_tmp - mx_start_tmp);
				zoom_j2 = zoom_j2_pan - (my_tmp - my_start_tmp);
			} else if (releasing && !shift_down) {
				if (mx == mx_start && my == my_start) {
					resetZoom();
				} else {
					zoom_i1 = Math.min(mx_start, mx);
					zoom_j1 = Math.min(my_start, my);
					zoom_i2 = Math.max(mx_start, mx);
					zoom_j2 = Math.max(my_start, my);
					zoomed = true;
				}
			}
			break;
		case FLOODSELECT:
		case SELECT:
			if (pressing) {
				if (brush == Brush.FLOODSELECT && !moving_selection) {

					boolean isselected = selected[mx][my];
					this.floodFill(mx, my, new FloodFillFunc() {
						@Override
						public boolean isValid(int i, int j) {
							return e.materials[i][j].type == e.materials[mx][my].type && selected[i][j] == isselected;
						}

						@Override
						public void fill(int i, int j) {
							selected[i][j] = !isselected;
						}
					});
				}
				else if (moving_selection && !dragging_selection) {
					selection.paste(delta_mx, delta_my); 
					flagChanges(true);
					moving_selection = false;
					dragging_selection = false;
				} else if (!moving_selection && selected[mx][my]) {
					selection.cut();
					flagChanges(true);
					moving_selection = true;
					dragging_selection = true;
					delta_mx = 0;
					delta_my = 0;
				}
			} else if (mouse_pressed) {
				if (brush != Brush.FLOODSELECT) {
					if (dragging_selection) {
						delta_mx = mx - mx_start;
						delta_my = my - my_start;
					} else {
						int mx0 = Math.min(mx_start, mx);
						int my0 = Math.min(my_start, my);
						int mx1 = Math.max(mx_start, mx);
						int my1 = Math.max(my_start, my);

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
						
						for (Probe p : e.probes) {
							p.selected = p.intersects(mx0, my0, mx1, my1);
						}
					}
				}
			} else if (releasing) {
				if (brush != Brush.FLOODSELECT) {
					delta_mx = mx - mx_start;
					delta_my = my - my_start;
					if (dragging_selection) {
						selection.paste(delta_mx, delta_my); 
						moving_selection = false;
						dragging_selection = false;
						flagChanges(true);
					} else {
						if (delta_mx == 0 && delta_my == 0) {
							for (int i = 0; i < e.nx; i++)
							{
								for (int j = 0; j < e.ny; j++)
								{
									selected[i][j] = false;
								}
							}
							e.probes.forEach((p) -> p.selected = false);
						}
					}
				}
			} else {
				delta_mx = mx;
				delta_my = my;
			}
			break;
		case CURRENT:
			if (pressing) {
				CurrentProbe p = new CurrentProbe(mx_start, my_start);
				p.name = e.getProbeName(); e.probe_index++;
				e.addProbe(p);
			} else if (mouse_pressed) {
				e.probes.get(e.probes.size()-1).drag(mx, my);
			} else if (releasing) {
				flagChanges(false);
			}
			break;
		case VOLTAGE:
			if (pressing) {
				VoltageProbe p = new VoltageProbe(mx_start, my_start);
				p.name = e.getProbeName(); e.probe_index++;
				e.addProbe(p);
			} else if (mouse_pressed) {
				e.probes.get(e.probes.size()-1).drag(mx, my);
			} else if (releasing) {
				flagChanges(false);
			}
			break;
		case CHARGE:
			if (pressing) {
				ChargeProbe p = new ChargeProbe(mx_start, my_start);
				p.name = e.getProbeName(); e.probe_index++;
				e.addProbe(p);
			} else if (mouse_pressed) {
				e.probes.get(e.probes.size()-1).drag(mx, my);
			} else if (releasing) {
				flagChanges(false);
			}
			break;
		case FLUX:
			if (pressing) {
				FluxProbe p = new FluxProbe(mx_start, my_start);
				p.name = e.getProbeName(); e.probe_index++;
				e.addProbe(p);
			} else if (mouse_pressed) {
				e.probes.get(e.probes.size()-1).drag(mx, my);
			} else if (releasing) {
				flagChanges(false);
			}
			break;
		case PROBE:
			if (pressing) {
				Probe p = null;
				switch ((CustProbeType) e.opts.gui_probetype.getSelectedItem()) {
				case AREA:
					p = new AreaProbe(mx_start, my_start);
					break;
				case LINE:
					p = new LineProbe(mx_start, my_start);
					break;
				case POINT:
					p = new PointProbe(mx_start, my_start);
					break;
				default:
					break;
				}
				if (p != null) {
					p.name = e.getProbeName(); e.probe_index++;
					e.addProbe(p);
				}
			} else if (mouse_pressed) {
				Probe p = e.probes.get(e.probes.size()-1);
				if (p != null) {
					p.drag(mx, my);
				}
			} else if (releasing) {
				Probe p = e.probes.get(e.probes.size()-1);
				if (p != null) {
					if (p instanceof AreaProbe) {

						JList<ScalarView> tmplist = new JList<>(ScalarView.values());

						tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
						tmplist.setVisibleRowCount(5);
						tmplist.setSelectedValue(Preset.DEFAULT, true);
						JScrollPane scrollPane = new JScrollPane(tmplist);
						int result = JOptionPane.showConfirmDialog(null, scrollPane, "Select quantity", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

						if (result == JOptionPane.OK_OPTION) {
							ScalarView selected = tmplist.getSelectedValue();

							if (selected != null) {
								Quantity quantity = selected.unit;
								QuantityType quantitytype = QuantityType.SCALAR;
								
								Quantity new_quantity = quantity.multiplyArea();
								if (new_quantity != null) {
									quantity = new_quantity;
									quantitytype = QuantityType.FLUX_DENSITY;
								} else {
									new_quantity = quantity.multiplyVolume();
									if (new_quantity != null) {
										quantity = new_quantity;
										quantitytype = QuantityType.DENSITY;
									} 
								}
								
								if (quantity != null)
								{
									((AreaProbe)p).scalarname = selected;
									((AreaProbe)p).shorthand = quantity.shorthand;
									((AreaProbe)p).quantity = quantity;
									((AreaProbe)p).quantitytype = quantitytype;
									((AreaProbe)p).user_placed = true;
								}
							}
						} else {
							e.removeProbe(p);
						}
					} else if (p instanceof PointProbe) {
						JList<ScalarView> tmplist = new JList<>(ScalarView.values());

						tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
						tmplist.setVisibleRowCount(5);
						tmplist.setSelectedValue(Preset.DEFAULT, true);
						JScrollPane scrollPane = new JScrollPane(tmplist);
						int result = JOptionPane.showConfirmDialog(null, scrollPane, "Select quantity", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

						if (result == JOptionPane.OK_OPTION) {
							ScalarView selected = tmplist.getSelectedValue();

							if (selected != null) {
								((PointProbe)p).scalarname = selected;
								((PointProbe)p).shorthand = selected.shorthand;
								((PointProbe)p).quantity = selected.unit;
							}
						} else {
							e.removeProbe(p);
						}
					} else if (p instanceof LineProbe) {
						JList<VectorView> tmplist = new JList<>(VectorView.values());

						tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
						tmplist.setVisibleRowCount(5);
						tmplist.setSelectedValue(Preset.DEFAULT, true);
						JScrollPane scrollPane = new JScrollPane(tmplist);
						int result = JOptionPane.showConfirmDialog(null, scrollPane, "Select quantity", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

						if (result == JOptionPane.OK_OPTION) {
							VectorView selected = tmplist.getSelectedValue();

							if (selected != null) {
								Quantity quantity = selected.unit;
								quantity = quantity.multiplyArea();
								if (quantity != null)
								{
									((LineProbe)p).vectorname = selected;
									((LineProbe)p).shorthand = quantity.shorthand;
									((LineProbe)p).quantity = quantity;
								} else {
									e.removeProbe(p);
								}
							}
						} else {
							e.removeProbe(p);
						}
					}
				}
				flagChanges(false);
			}
			break;
		case DELETEPROBE:
			e.canvas.setCursor(DEFAULT_CURSOR);
			int i = 0;
			while(i < e.probes.size()) {
				Probe p = e.probes.get(i);
				if (p.isMouseHovering(mx, my)) {
					e.canvas.setCursor(HAND_CURSOR);
					if (pressing) {
						e.removeProbe(p);
						flagChanges(false);
						i--;
					}
				}
				i++;
			}
			break;
		case MOVELABEL:
			if (pressing) {
				Probe.LabelCoord coord = null;
				for (Probe p : e.probes)
					coord = selectLabel(p.labelcoord, coord);
				
				if (coord != null)
					labelcoord = coord;
			} else if (mouse_pressed) {
				if (labelcoord != null) {
					labelcoord.x = mx;
					labelcoord.y = my;
				}
			} else if (releasing) {
				labelcoord = null;
				flagChanges(false);
			} else {
				e.canvas.setCursor(DEFAULT_CURSOR);
				Probe.LabelCoord coord = null;
				for (Probe p : e.probes)
					coord = selectLabel(p.labelcoord, coord);
				
				if (coord != null)
					e.canvas.setCursor(HAND_CURSOR);
			}
			break;
		case TEXT:
			e.canvas.setCursor(HAND_CURSOR);
			if (mouse_pressed) {
				startTextInput();
				text_x = mx;
				text_y = my;
			}
			break;
		case GROUND:
			if (pressing) {
				if (!e.hasGround()) {
					Ground ground = new Ground(mx_start, my_start);
					e.addProbe(ground);
				}
			} else if (mouse_pressed) {
				e.getGround().drag(mx, my);
			} else if (releasing) {
				flagChanges(false);
			}
			break;
		case RULER:
			if (pressing) {
				if (!e.hasRuler()) {
					Ruler p = new Ruler(mx_start, my_start);
					e.addProbe(p);
				} else {
					Ruler p = e.getRuler();
					p.x1 = mx_start;
					p.y1 = my_start;
				}
			} else if (mouse_pressed) {
				e.getRuler().drag(mx, my);
			} else if (releasing) {
				Ruler p = e.getRuler();
				if (p.x1 == p.x2 && p.y1 == p.y2) {
					e.removeProbe(e.getRuler());
				}
				flagChanges(false);
			}
			break;
		case BANDS:
			createPath();
			if (releasing && !shift_down) {
				e.bandplot.createPlot(e, plotpath);
				plotpath = null;
			}
			break;
		case SCALARPLOT:
			createPath();
			if (releasing && !shift_down) {
				e.scalarplot.createPlot(e, plotpath);
				plotpath = null;
			}
			break;
		case CARRIERPLOT:
			createPath();
			if (releasing && !shift_down) {
				e.carrierplot.createPlot(e, plotpath);
				plotpath = null;
			}
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

		if ((releasing && Brush.isMaterialModifyingBrush(brush)) || (BoundaryCondition)e.opts.gui_bc.getSelectedItem() != prev_boundary) {
			flagChanges(true);
		}
		prev_boundary = (BoundaryCondition)e.opts.gui_bc.getSelectedItem();
		
		if (updatematerials) {
			e.updateAllMaterials(false);
			e.multigridSolve(true, false);
			
			updatematerials = false;
		}

		mxp = mx;
		myp = my;
	}

	public void createPath() {
		if (plotpath == null) {
			if (pressing) {
				if (shift_down) {
					plotpath = new SegmentedPath();
					((SegmentedPath) plotpath).addNewJoint(mx_start, my_start);
				} else {
					plotpath = new LinePath();
					((LinePath) plotpath).x1 = mx_start;
					((LinePath) plotpath).y1 = my_start;
					((LinePath) plotpath).x2 = mx_start;
					((LinePath) plotpath).y2 = my_start;
				}
			}
		} else if (plotpath instanceof LinePath) {
			if (mouse_pressed) {
				((LinePath) plotpath).x2 = mx;
				((LinePath) plotpath).y2 = my;
			}
		} else if (plotpath instanceof SegmentedPath) {
			if (pressing) {
				((SegmentedPath) plotpath).addNewJoint(mx, my);
			}
		}
	}
	
	public boolean updatematerials = false;
	
	public void flagChanges(boolean updateMaterials) {
		changesmade = true;
		if (!e.opts.getTitle().endsWith("*"))
			e.opts.setTitle(e.opts.getTitle() + " *");

		undoredo.captureState(e);
		
		updatematerials = updatematerials || updateMaterials;
	}
	
	public Probe.LabelCoord selectLabel(Probe.LabelCoord c_in, Probe.LabelCoord c_opt) {
		if (c_opt == null) {
			if (Math.abs(c_in.x + 10 - mx) < 10 && Math.abs(c_in.y - 2 - my) < 4) {
				return c_in;
			} else {
				return null;
			}
		}
		
		return c_opt;
		//if (Utils.length(c_in.x - mx, c_in.y - my) < Utils.length(c_opt.x - mx, c_opt.y - my)) {
		//	return c_in;
		//} else {
		//	return c_opt;
		//}
	}
	
	public void resetZoom() {
		zoom_i1 = 0;
		zoom_j1 = 0;
		zoom_i2 = e.nx-1;
		zoom_j2 = e.ny-1;
		zoomed = false;
	}
	
	public void startTextInput() {
		if (!texting) {
			removeKeyBinds(e.canvas);
			removeKeyBinds(e.opts.panel);
			e.opts.menuBar.setEnabled(false);
			texting = true;
		}
	}

	public void endTextInput() {
		if (texting) {
			addKeyBinds(e.canvas);
			addKeyBinds(e.opts.panel);
			e.opts.menuBar.setEnabled(true);
			texting = false;
		}
	}
	
	public void drawMaterialLine(double x1, double y1, double x2, double y2, Brush brush, BrushShape brushshape, GeneralMaterialType mat, double brushsize, double EMF_angle) {
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
				cx = i;
				cy = j;

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
					if (mat.type == MaterialType.VACUUM) {
						e.eraseMaterial(i, j);
					} else if (e.materials[i][j].type == MaterialType.VACUUM ^ brush == Brush.REPLACE) {
						e.eraseMaterial(i, j);
						e.initializeMaterial(i, j, mat);
						if (mat.type.hasEMF()) e.materials[i][j].emf_direction = EMF_angle;
					}
				}
			}
		}
	}
	
	public void setEMFs() {
		int EMF_setting = e.opts.gui_parameter3.getValue();
		
		double new_EMF = 0;
		if (!iscurrentselected) {
			new_EMF = e.max_EMF*EMF_setting/50.0;
		} else {
			new_EMF = (e.max_current/e.currentsource_sigma)*EMF_setting/50.0;
		}

		if (e.opts.gui_brush.getSelectedItem() == Brush.INTERACT && EMF_selected) {
			e.opts.gui_parameter3.setVisible(true);
			e.opts.gui_parameter3_text.setVisible(true);
			if (!iscurrentselected)
				e.opts.gui_parameter3_text.setText("EMF: " + e.units.toString(new_EMF, Quantity.ELECTRIC_FIELD));
			else
				e.opts.gui_parameter3_text.setText("J: " + e.units.toString(new_EMF*e.currentsource_sigma, Quantity.CURRENT_DENSITY));
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
					if (e.materials[i][j].type.hasEMF() && selected_EMF[i][j]) {
						e.materials[i][j].emf = new_EMF;
					}
				}
			}

			e.updateJustEMFs();
		}

		prev_EMF_setting = EMF_setting;
		

		e.AC_freq = e.default_AC_freq*Math.pow(10, e.opts.gui_parameter1.getValue()/10.0);
		if (e.AC_source_exists) {
			e.opts.gui_parameter1.setEnabled(true);
			e.opts.gui_parameter1.setVisible(true);
			e.opts.gui_parameter1_text.setVisible(true);
		} else {
			e.opts.gui_parameter1.setEnabled(false);
			e.opts.gui_parameter1.setVisible(false);
			e.opts.gui_parameter1_text.setVisible(false);
		}
		
		e.opts.gui_parameter1_text.setText("AC Freq: " + e.units.toString(e.AC_freq, Quantity.FREQUENCY));

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
		else if (ev.getSource() == e.opts.menu_saveas)
			saveas = true;
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
		else if (ev.getSource() == e.opts.menu_github)
			try {
				java.awt.Desktop.getDesktop().browse(new URI("https://github.com/StunningLlama/SemiSim/tree/SemiSim"));
			} catch (IOException | URISyntaxException ex) {
				ex.printStackTrace();
			}
		else if (ev.getSource() == e.opts.menu_editdesc) {
			SwingUtilities.invokeLater(() -> {
				new DescDialog();
			});
		} else if (ev.getSource() == e.opts.gui_brush) {

			Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();

			if (brush == Brush.PROBEPLOT) {
				plotinterval = e.opts.gui_plotinterval.getValue();
				e.opts.gui_plotinterval_text.setText("Time resolution: " + e.units.toString(plotinterval*e.dt*e.iteration_multiplier, Quantity.TIME));
				for (Plot p : e.plots) {
					if (p instanceof ProbePlot) {
						if (!p.frame.isVisible() && ((ProbePlot)p).checkProbesExist(e))
							p.createPlot(e, null);
					}
				}
			} else if (brush == Brush.XYPLOT) {
				e.xyplot.createPlot(e, null);
				
				JList<Probe> tmplist = new JList<>(e.probes.toArray(new Probe[0]));

				tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
				tmplist.setVisibleRowCount(5);
				tmplist.setSelectedValue(Preset.DEFAULT, true);
				JScrollPane scrollPane = new JScrollPane(tmplist);
				int result = JOptionPane.showConfirmDialog(null, scrollPane, "Select X variable", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

				if (result == JOptionPane.OK_OPTION) {
					Probe selected = tmplist.getSelectedValue();

					if (selected != null) {
						e.xyplot.x = selected;
					}
				}
				
				tmplist = new JList<>(e.probes.toArray(new Probe[0]));

				tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
				tmplist.setVisibleRowCount(5);
				tmplist.setSelectedValue(Preset.DEFAULT, true);
				scrollPane = new JScrollPane(tmplist);
				result = JOptionPane.showConfirmDialog(null, scrollPane, "Select Y variable", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

				if (result == JOptionPane.OK_OPTION) {
					Probe selected = tmplist.getSelectedValue();

					if (selected != null) {
						e.xyplot.y = selected;
					}
				}
			}
			e.controls.brush_changed = true;
		} else if (ev.getSource() == e.opts.menu_advancedsettings) {
			e.adv_opts.storeAdvancedSettings();
			e.adv_opts.setVisible(true);
		} else if (ev.getSource() == e.opts.gui_carriers) {
			e.opts.menu_carriers.setSelected(e.opts.gui_carriers.isSelected());
		} else if (ev.getSource() == e.opts.menu_carriers) {
			e.opts.gui_carriers.setSelected(e.opts.menu_carriers.isSelected());
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
		} else if (ev.getSource() == e.opts.menu_undo) {
			undo = true;
		} else if (ev.getSource() == e.opts.menu_redo) {
			redo = true;
		} else if (ev.getSource() == e.opts.menu_selectall) {
			selectall = true;
		} else if (ev.getSource() == e.opts.menu_deselectall) {
			deselectall = true;
		} else if (ev.getSource() == e.opts.menu_about) {
			JOptionPane.showMessageDialog(e.opts, SemiSim.about, "About", JOptionPane.INFORMATION_MESSAGE);
		} else if (ev.getSource() == e.opts.menu_report) {
			JOptionPane.showMessageDialog(e.opts, "<html><body><p style='width: 300px;'>Please contact Brandon at brandonli.lex@gmail.com or go to https://github.com/StunningLlama/SemiSim/issues.</p></body></html>", "Report a bug", JOptionPane.INFORMATION_MESSAGE);
		} else if (ev.getSource() == e.opts.menu_pref) {
			e.prefs.getPrefs();
			e.prefs.setVisible(true);
		} else if (ev.getSource() == e.opts.menu_debug) {
			key_dbg.actionPerformed(null);
		} else if (ev.getSource() == e.opts.menu_exit) {
			exit = true;
		}else if (ev.getSource() == e.opts.menu_cust_material) {
			e.materialmanager.valueChanged(null);
			e.materialmanager.setVisible(true);
		} else if (ev.getSource() == e.opts.menu_view_materials) {
			e.materialviewer.updateUI();
			e.materialviewer.setVisible(true);
		} else if (ev.getSource() instanceof JRadioButtonMenuItem) {
			if (scalarview.containsButton((JRadioButtonMenuItem) ev.getSource()) != null || vectorview.containsButton((JRadioButtonMenuItem) ev.getSource()) != null)
				e.updateMiscFields = true;
			if (brushes.containsButton((JRadioButtonMenuItem) ev.getSource()) != null)
				e.opts.gui_brush.setSelectedItem(brushes.getOption());
		}
	}
	
	@Override
	public void itemStateChanged(ItemEvent ev) {
		if (ev.getSource() == e.opts.gui_brush) {
			brushes.setOption((Brush) e.opts.gui_brush.getSelectedItem());
		}
	}
	
	@Override
	public void mouseClicked(MouseEvent arg0) {}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent ev) {
		mouse_pressed = true;
		mousebutton = ev.getButton();
		mx_screen = ev.getX();
		my_screen = ev.getY();
		mx_start_screen = ev.getX();
		my_start_screen = ev.getY();

	}
	@Override
	public void mouseReleased(MouseEvent ev) {
		mouse_pressed = false;
	}

	@Override
	public void mouseDragged(MouseEvent ev) {
		mx_screen = ev.getX();
		my_screen = ev.getY();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		mx_screen = arg0.getX();
		my_screen = arg0.getY();
	}
	
	public void snapToCardinals(int mx, int my) {
		int dx = mx - mx_start;
		int dy = my - my_start;
		int min = Math.abs(dx) > Math.abs(dy)? dx : dy;
		int[] xc = {min, 0, min, min, -min, -min};
		int[] yc = {0, min, min, -min, min, -min};
		int dmin = Integer.MAX_VALUE;
		int imin = -1;
		
		for (int i = 0; i < 6; i++) {
			int d = (xc[i]-dx)*(xc[i]-dx) + (yc[i]-dy)*(yc[i]-dy);
			if (d < dmin) {
				dmin = d;
				imin = i;
			}
		}
		
		if (imin != -1) {
			this.mx = mx_start + xc[imin];
			this.my = my_start + yc[imin];
		}
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
			
			if (debugging) {
				e.controls.scalarview.addOption(ScalarView.DEBUG);
			} else {
				e.opts.textPane.setText(e.description);
				e.controls.scalarview.removeOption(ScalarView.DEBUG);
			}
			if (ev != null)
				e.opts.menu_debug.setSelected(debugging);
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
    
    private Action key_selectall = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		selectall = true;
        }
    };
    
    private Action key_deselectall = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		deselectall = true;
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
    		e.opts.menu_elem_colors.setSelected(!e.opts.menu_elem_colors.isSelected());
        }
    };

    ScalarMode prev_scalar_mode = ScalarMode.NONE;
    VectorMode prev_vector_mode = VectorMode.NONE;

    private Action key_scalar_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		if (e.controls.scalarmode.getOption() == ScalarMode.NONE)
    			e.controls.scalarmode.setOption(prev_scalar_mode);
    		else
    		{
    			prev_scalar_mode = e.controls.scalarmode.getOption();
    			e.controls.scalarmode.setOption(ScalarMode.NONE);
    		}
        }
    };

    private Action key_vector_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		if (e.controls.vectormode.getOption()  == VectorMode.NONE)
    			e.controls.vectormode.setOption(prev_vector_mode);
    		else
    		{
    			prev_vector_mode = e.controls.vectormode.getOption();
    			e.controls.vectormode.setOption(VectorMode.NONE);
    		}
        }
    };

    private Action key_tooltip = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.menu_tooltip.setSelected(!e.opts.menu_tooltip.isSelected());
        }
    };

    private Action key_textbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		e.opts.menu_text_bg.setSelected(!e.opts.menu_text_bg.isSelected());
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
			e.opts.menu_interface.setSelected(!e.opts.menu_interface.isSelected());
        }
    };
    
    private Action key_save = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		save = true;
        }
    };
    
    private Action key_saveas = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent ev) {
    		saveas = true;
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
    		//e.opts.gui_brush.setSelectedItem(Brush.SELECT);
			test = !test;
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

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, InputEvent.SHIFT_DOWN_MASK | InputEvent.ALT_DOWN_MASK), key_shift);
    	contentPane.getActionMap().put(key_shift, key_shift);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, InputEvent.ALT_DOWN_MASK, true), key_shift_up);
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
    	
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK), key_saveas);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.META_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK), key_saveas);
    	contentPane.getActionMap().put(key_saveas, key_saveas);

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

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_DOWN_MASK), key_selectall);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.META_DOWN_MASK), key_selectall);
    	contentPane.getActionMap().put(key_selectall, key_selectall);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_DOWN_MASK), key_deselectall);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.META_DOWN_MASK), key_deselectall);
    	contentPane.getActionMap().put(key_deselectall, key_deselectall);

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

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, InputEvent.ALT_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK), key_alt);
    	contentPane.getActionMap().put(key_alt, key_alt);

    	map2.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT,  InputEvent.SHIFT_DOWN_MASK, true), key_alt_up);
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
		e.opts.gui_brushsize.setValue(e.opts.gui_brushsize.getValue() - (int)(5*ev.getPreciseWheelRotation()));
	}
	

	@Override
	public void keyTyped(KeyEvent ev) {
		if (texting) {
			if (Font7x5.getCharacter(ev.getKeyChar()) != null) {
				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 7; j++) {
						if (text_x+i+1 >= 0 && text_x+i+1 < e.nx && text_y+j >= 0 && text_y+j < e.ny
								&& Font7x5.getPixel(ev.getKeyChar(), 4-i, j) == 1 && e.materials[text_x+i+1][text_y+j].type == MaterialType.VACUUM) {
							e.initializeMaterial(e.materials[text_x+i+1][text_y+j], MaterialType.DECO);
						}
					}
				}
				text_x += 6;
				flagChanges(false);
				e.updateAllMaterials(false);
			}
		}
	}

	@Override
	public void keyPressed(KeyEvent ev) {
		if (texting) {
			if (ev.getKeyCode() == KeyEvent.VK_ESCAPE || ev.getKeyCode() == KeyEvent.VK_ENTER) {
				endTextInput();
				return;
			}
			
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
				flagChanges(false);
				e.updateAllMaterials(false);
			}
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {}
	
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
		ZOOM("Zoom and Pan"),
		DRAW("Draw"),
		REPLACE("Replace"),
		LINE("Line"),
		FILL("Fill"),
		ERASE("Eraser"),
		SELECT("Select and Move"),
		FLOODSELECT("Flood select"),
		TEXT("Text"),
		VOLTAGE("Voltage probe"),
		CURRENT("Current probe"),
		CHARGE("Charge probe"),
		FLUX("Magnetic flux probe"),
		GROUND("Ground"),
		PROBE("Custom probe"),
		DELETEPROBE("Delete probe"),
		MOVELABEL("Move label"),
		BANDS("Plot bands"),
		SCALARPLOT("Plot scalar field"),
		CARRIERPLOT("Plot carriers"),
		PROBEPLOT("Plot probe data"),
		XYPLOT("XY plot"),
		RULER("Ruler");
	
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
	
	public enum CustProbeType {
		POINT("Type: point"), LINE("Type: line"), AREA("Type: area");

		public String name;
		CustProbeType(String name)
		{
			this.name = name;
		}

		@Override
		public String toString() {
			return name;
		}
	}

	class DescDialog extends JDialog {
		private static final long serialVersionUID = -1338819833865147129L;
		public JButton okButton = new JButton("Apply");
		public JButton cancelButton = new JButton("Cancel");
		
		public DescDialog() {
			String text = e.description;
			
			getContentPane().setLayout(new BorderLayout());
			JTextArea area = new JTextArea();
			JScrollPane scroll = new JScrollPane(area);
			setSize(500, 500);
			area.setText(text);
			area.setEditable(true);
			area.setLineWrap(true);
			area.setWrapStyleWord(true);
			getContentPane().add(scroll);

			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			
			buttonPane.add(okButton);
			getRootPane().setDefaultButton(okButton);
			
			buttonPane.add(cancelButton);
			
			okButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent ev) {
					e.description = area.getText();
					e.opts.textPane.setText(e.description);
					e.opts.textPane.setEditable(false);
					flagChanges(false);
					dispose();
				}
			});
			
			cancelButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent ev) {
					dispose();
				}
			});
			setVisible(true);
		}
	}

	interface FloodFillFunc {
		boolean isValid(int i, int j);
		void fill(int i, int j);
	}

	@Override
	public void windowOpened(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent ev) {
		if (e.controls.changesmade) {
			String[] options = {"Yes", "No"};
			int result = JOptionPane.showOptionDialog(e.opts, "There are unsaved changes. Do you still wish to quit?", "Message", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
			if (result == JOptionPane.OK_OPTION)
			{
				e.opts.dispose();
				System.exit(0);
			}
		} else {
			e.opts.dispose();
			System.exit(0);
		}
	}

	@Override
	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowIconified(WindowEvent e) {}

	@Override
	public void windowDeiconified(WindowEvent e) {}

	@Override
	public void windowActivated(WindowEvent e) {}

	@Override
	public void windowDeactivated(WindowEvent e) {}
}
