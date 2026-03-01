// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
import java.awt.Adjustable;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ListCellRenderer;
import javax.swing.ScrollPaneConstants;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import java.awt.event.KeyEvent;
import java.io.File;
import java.awt.event.InputEvent;

import electrodynamics.Controls.Brush;
import electrodynamics.Renderer.ScalarMode;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorMode;
import electrodynamics.Renderer.VectorView;
import electrodynamics.Simulation.BoundaryCondition;
import electrodynamics.util.MenuBuilder;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.SwingConstants;

public class MainWindow extends JFrame {

	/**
	 *
	 */
	private static final long serialVersionUID = -5756219569007074449L;

	public JPanel contentPane;
	public JButton gui_reset;
	public JScrollBar gui_simspeed;
	public JScrollBar gui_brightness;
	public JScrollBar gui_brushsize;
	public JScrollBar gui_parameter1;
	public JCheckBox gui_paused;
	public JPanel panel;
	public JLabel gui_parameter2_text;
	public JScrollBar gui_parameter2;
	public JLabel lblVectorBrightness;
	public JScrollBar gui_brightness_vec;
	public JScrollBar gui_simspeed_2;
	public JComboBox gui_brush_1;
	public JScrollBar gui_parameter3;
	public JTextArea textPane;
	public JScrollPane scrollPane;
	public JLabel gui_parameter1_text;
	public JLabel gui_parameter3_text;
	public JLabel gui_stepslbl;
	public JLabel gui_stepsizelbl;
	public JCheckBox gui_brush_highlight;
	public JComboBox gui_material;
	public JComboBox gui_brush;
	public JComboBox gui_bc;
	public JLabel lblBrushSize;
	public JMenuItem menu_open;
	public JMenuItem menu_save;
	public JMenuItem menu_about;
	public JMenuItem menu_help;
	public JMenuItem menu_cut;
	public JMenuItem menu_copy;
	public JMenuItem menu_paste;
	public JMenuItem menu_undo;
	public JMenuItem menu_redo;
	public JMenu menu_examples;
	public JMenu menu_tools;
	public JMenuItem menu_editdesc;
	public JMenuItem menu_new;
	public JMenuItem menu_rotate;
	public JMenuItem menu_flip_h;
	public JMenuItem menu_flip_v;
	public JSeparator separator_2;
	public JLabel gui_carrierlbl;
	public JScrollBar gui_carrier_density;
	public JMenu menu_view;
	public JCheckBox gui_carriers;
	public JMenuItem menu_img;
	public JMenu mnNewMenu_1;
	public JCheckBoxMenuItem menu_materialname;
	public JCheckBoxMenuItem menu_interface;
	public JCheckBoxMenuItem menu_tooltip;
	public JCheckBoxMenuItem menu_text_bg;
	public JCheckBoxMenuItem menu_elem_colors;
	public JCheckBoxMenuItem menu_borders;
	public JCheckBoxMenuItem menu_carriers;
	public JCheckBoxMenuItem menu_probes;
	public JCheckBoxMenuItem menu_time;
	public JMenuItem menu_advancedsettings;
	public JCheckBoxMenuItem menu_debug;
	public JScrollBar gui_plotinterval;
	public JLabel gui_plotinterval_text;
	public JMenuBar menuBar;
	private JLabel lblNewLabel;

	/**
	 * Create the frame.
	 */
	public MainWindow() {
		setTitle(SemiSim.name);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 590, 832);
		
		menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		JMenu mnNewMenu = new JMenu("File");
		menuBar.add(mnNewMenu);
		
		menu_new = new JMenuItem("New simulation...");
		menu_new.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, InputEvent.CTRL_DOWN_MASK));
		mnNewMenu.add(menu_new);
		
		menu_open = new JMenuItem("Open file...");
		menu_open.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_DOWN_MASK));
		mnNewMenu.add(menu_open);
		
		menu_save = new JMenuItem("Save as...");
		menu_save.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK));
		mnNewMenu.add(menu_save);
		
		menu_editdesc = new JMenuItem("Edit description...");
		mnNewMenu.add(menu_editdesc);
		
		JSeparator separator_1 = new JSeparator();
		mnNewMenu.add(separator_1);
		
		menu_about = new JMenuItem("About");
		mnNewMenu.add(menu_about);
		
		menu_help = new JMenuItem("Open manual");
		mnNewMenu.add(menu_help);
		
		JMenu menu_asdf = new JMenu("Edit");
		menuBar.add(menu_asdf);
		
		menu_undo = new JMenuItem("Undo");
		menu_undo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_undo);
		
		menu_redo = new JMenuItem("Redo");
		menu_redo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK));
		menu_asdf.add(menu_redo);
		
		JSeparator separator = new JSeparator();
		menu_asdf.add(separator);
		
		menu_cut = new JMenuItem("Cut");
		menu_cut.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_cut);
		
		menu_copy = new JMenuItem("Copy");
		menu_copy.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_copy);
		
		menu_paste = new JMenuItem("Paste");
		menu_paste.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_paste);
		
		menu_rotate = new JMenuItem("Rotate");
		menu_rotate.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_rotate);
		
		menu_flip_h = new JMenuItem("Flip horizontally");
		menu_flip_h.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_flip_h);
		
		menu_flip_v = new JMenuItem("Flip vertically");
		menu_flip_v.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_flip_v);
		
		separator_2 = new JSeparator();
		menu_asdf.add(separator_2);
		
		menu_img = new JMenuItem("Set image size");
		menu_asdf.add(menu_img);
		
		menu_advancedsettings = new JMenuItem("Advanced settings");
		menu_asdf.add(menu_advancedsettings);
		
		menu_tools = new JMenu("Tools");
		menuBar.add(menu_tools);
		
		menu_view = new JMenu("View");
		menuBar.add(menu_view);
		
		mnNewMenu_1 = new JMenu("Graphics");
		menuBar.add(mnNewMenu_1);
		
		menu_interface = new JCheckBoxMenuItem("Display interface");
		menu_interface.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_H, 0));
		menu_interface.setSelected(true);
		mnNewMenu_1.add(menu_interface);
		
		menu_time = new JCheckBoxMenuItem("Show time");
		menu_time.setSelected(true);
		mnNewMenu_1.add(menu_time);
		
		menu_materialname = new JCheckBoxMenuItem("Show material name");
		menu_materialname.setSelected(true);
		mnNewMenu_1.add(menu_materialname);
		
		menu_tooltip = new JCheckBoxMenuItem("Show simulation variables");
		menu_tooltip.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, 0));
		mnNewMenu_1.add(menu_tooltip);
		
		menu_probes = new JCheckBoxMenuItem("Show probes");
		menu_probes.setSelected(true);
		mnNewMenu_1.add(menu_probes);
		
		menu_elem_colors = new JCheckBoxMenuItem("Show material colors");
		menu_elem_colors.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, 0));
		menu_elem_colors.setSelected(true);
		mnNewMenu_1.add(menu_elem_colors);
		
		menu_borders = new JCheckBoxMenuItem("Show material borders");
		menu_borders.setSelected(true);
		mnNewMenu_1.add(menu_borders);
		
		menu_carriers = new JCheckBoxMenuItem("Show charge carriers");
		mnNewMenu_1.add(menu_carriers);
		
		menu_text_bg = new JCheckBoxMenuItem("Show text background");
		menu_text_bg.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G, 0));
		menu_text_bg.setSelected(true);
		mnNewMenu_1.add(menu_text_bg);
		
		menu_debug = new JCheckBoxMenuItem("Debug mode");
		menu_debug.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, 0));
		mnNewMenu_1.add(menu_debug);
		
		menu_examples = new JMenu("Examples");
		menuBar.add(menu_examples);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);

		panel = new JPanel();
		panel.setPreferredSize(new Dimension(380, 200));
		panel.setMinimumSize(new Dimension(200, 200));
		contentPane.add(panel, BorderLayout.EAST);
		panel.setLayout(null);

		gui_reset = new JButton("Reset fields");
		gui_reset.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
			}
		});
		gui_reset.setBounds(11, 42, 171, 23);
		panel.add(gui_reset);

		gui_paused = new JCheckBox("Pause (P)");
		gui_paused.setSelected(false);
		gui_paused.setBounds(11, 78, 101, 23);
		panel.add(gui_paused);

		gui_brushsize = new JScrollBar();
		gui_brushsize.setMaximum(750);
		gui_brushsize.setValue(250);
		gui_brushsize.setOrientation(Adjustable.HORIZONTAL);
		gui_brushsize.setBounds(202, 238, 171, 17);
		panel.add(gui_brushsize);

		gui_simspeed = new JScrollBar();
		gui_simspeed.setValue(20);
		gui_simspeed.setBlockIncrement(1);
		gui_simspeed.setMaximum(30);
		gui_simspeed.setOrientation(Adjustable.HORIZONTAL);
		gui_simspeed.setBounds(11, 137, 171, 17);
		panel.add(gui_simspeed);

		gui_brightness = new JScrollBar();
		gui_brightness.setValue(-20);
		gui_brightness.setBlockIncrement(1);
		gui_brightness.setMinimum(-45);
		gui_brightness.setMaximum(45);
		gui_brightness.setOrientation(Adjustable.HORIZONTAL);
		gui_brightness.setBounds(10, 238, 171, 17);
		panel.add(gui_brightness);

		gui_parameter1 = new JScrollBar();
		gui_parameter1.setVisible(false);
		gui_parameter1.setEnabled(false);
		gui_parameter1.setMaximum(25);
		gui_parameter1.setMinimum(-15);
		gui_parameter1.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter1.setBounds(202, 289, 171, 17);
		panel.add(gui_parameter1);

		gui_brush = new JComboBox();
		gui_brush.setMaximumRowCount(16);
		gui_brush.setModel(new DefaultComboBoxModel(Controls.Brush.values()));
		gui_brush.setSelectedIndex(0);
		gui_brush.setBounds(202, 74, 171, 22);
		panel.add(gui_brush);
		addTooltips(gui_brush);

		gui_stepsizelbl = new JLabel("Timestep");
		gui_stepsizelbl.setBounds(21, 115, 161, 14);
		panel.add(gui_stepsizelbl);

		JLabel label5 = new JLabel("Scalar brightness");
		label5.setBounds(20, 218, 150, 14);
		panel.add(label5);

		lblBrushSize = new JLabel("Brush size");
		lblBrushSize.setBounds(212, 218, 138, 14);
		panel.add(lblBrushSize);

		gui_parameter1_text = new JLabel("");
		gui_parameter1_text.setEnabled(false);
		gui_parameter1_text.setBounds(212, 267, 154, 14);
		panel.add(gui_parameter1_text);

		gui_parameter2_text = new JLabel("Direction");
		gui_parameter2_text.setBounds(212, 322, 161, 14);
		panel.add(gui_parameter2_text);

		gui_parameter2 = new JScrollBar();
		gui_parameter2.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter2.setMaximum(34);
		gui_parameter2.setBounds(202, 341, 171, 17);
		panel.add(gui_parameter2);

		lblVectorBrightness = new JLabel("Vector field brightness");
		lblVectorBrightness.setBounds(20, 269, 150, 14);
		panel.add(lblVectorBrightness);

		gui_brightness_vec = new JScrollBar();
		gui_brightness_vec.setValue(-10);
		gui_brightness_vec.setOrientation(Adjustable.HORIZONTAL);
		gui_brightness_vec.setMinimum(-45);
		gui_brightness_vec.setMaximum(45);
		gui_brightness_vec.setBlockIncrement(1);
		gui_brightness_vec.setBounds(11, 290, 171, 17);
		panel.add(gui_brightness_vec);

		gui_brush_1 = new JComboBox();
		gui_brush_1.setMaximumRowCount(16);
		gui_brush_1.setModel(new DefaultComboBoxModel(Controls.BrushShape.values()));
		gui_brush_1.setSelectedIndex(1);
		gui_brush_1.setBounds(202, 139, 171, 22);
		panel.add(gui_brush_1);
		addTooltips(gui_brush_1);

		gui_stepslbl = new JLabel("Sim steps/frame");
		gui_stepslbl.setBounds(21, 166, 144, 14);
		panel.add(gui_stepslbl);

		gui_simspeed_2 = new JScrollBar();
		gui_simspeed_2.setMinimum(1);
		gui_simspeed_2.setValue(25);
		gui_simspeed_2.setOrientation(Adjustable.HORIZONTAL);
		gui_simspeed_2.setMaximum(110);
		gui_simspeed_2.setBlockIncrement(1);
		gui_simspeed_2.setBounds(11, 186, 171, 17);
		panel.add(gui_simspeed_2);

		gui_parameter3 = new JScrollBar();
		gui_parameter3.setMinimum(-50);
		gui_parameter3.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter3.setMaximum(60);
		gui_parameter3.setBlockIncrement(1);
		gui_parameter3.setBounds(202, 390, 171, 17);
		panel.add(gui_parameter3);

		gui_parameter3_text = new JLabel("EMF");
		gui_parameter3_text.setBounds(212, 369, 150, 14);
		panel.add(gui_parameter3_text);

		gui_brush_highlight = new JCheckBox("Brush highlight");
		gui_brush_highlight.setSelected(false);
		gui_brush_highlight.setBounds(201, 179, 160, 23);
		panel.add(gui_brush_highlight);

		gui_material = new JComboBox();
		gui_material.setMaximumRowCount(16);
		gui_material.setModel(new DefaultComboBoxModel(electrodynamics.MaterialType.values()));
		gui_material.setSelectedIndex(0);
		gui_material.setBounds(202, 106, 171, 22);
		panel.add(gui_material);
		addTooltips(gui_material);

		gui_bc = new JComboBox();
		gui_bc.setModel(new DefaultComboBoxModel(BoundaryCondition.values()));
		gui_bc.setSelectedIndex(0);
		gui_bc.setBounds(202, 42, 171, 22);
		panel.add(gui_bc);
		addTooltips(gui_bc);
		
		gui_carrierlbl = new JLabel("Charge carrier density");
		gui_carrierlbl.setBounds(21, 369, 150, 14);
		panel.add(gui_carrierlbl);
		
		gui_carrier_density = new JScrollBar();
		gui_carrier_density.setValue(-45);
		gui_carrier_density.setOrientation(JScrollBar.HORIZONTAL);
		gui_carrier_density.setMinimum(-45);
		gui_carrier_density.setMaximum(45);
		gui_carrier_density.setBlockIncrement(1);
		gui_carrier_density.setBounds(11, 390, 171, 17);
		panel.add(gui_carrier_density);
		
		gui_carriers = new JCheckBox("Show charge carriers");
		gui_carriers.setSelected(false);
		gui_carriers.setBounds(11, 326, 171, 23);
		panel.add(gui_carriers);
		
		JPanel panel_1 = new JPanel();
		panel_1.setBounds(27, 432, 327, 313);
		panel.add(panel_1);
		panel_1.setLayout(new BorderLayout(0, 0));

				textPane = new JTextArea();
				scrollPane = new JScrollPane(textPane);
				scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
				panel_1.add(scrollPane, BorderLayout.CENTER);
				
						textPane.setColumns(1);
						textPane.setWrapStyleWord(true);
						textPane.setRows(16);
						textPane.setText("Description of simulation");
						textPane.setFont(new Font("SansSerif", Font.PLAIN, 13));
						textPane.setMargin(new Insets(4, 4, 4, 4));
						textPane.setLineWrap(true);
						textPane.setEditable(false);
						
						gui_plotinterval_text = new JLabel("Probe plot interval");
						gui_plotinterval_text.setBounds(206, 322, 167, 14);
						panel.add(gui_plotinterval_text);
						
						gui_plotinterval = new JScrollBar();
						gui_plotinterval.setValue(10);
						gui_plotinterval.setMinimum(1);
						gui_plotinterval.setOrientation(JScrollBar.HORIZONTAL);
						gui_plotinterval.setMaximum(60);
						gui_plotinterval.setBounds(202, 341, 171, 17);
						panel.add(gui_plotinterval);
						
						lblNewLabel = new JLabel("Simulation controls");
						lblNewLabel.setFont(new Font("Lucida Grande", Font.BOLD, 13));
						lblNewLabel.setHorizontalAlignment(SwingConstants.CENTER);
						lblNewLabel.setBounds(115, 11, 150, 16);
						panel.add(lblNewLabel);
	}

	public void addTooltips(JComboBox box) {
		ListCellRenderer<? super Object> originalRenderer = box.getRenderer();

		box.setRenderer(new ListCellRenderer<Object>() {
		    @Override
		    public Component getListCellRendererComponent(JList<?> list, Object value, int index,
		                                                  boolean isSelected, boolean cellHasFocus) {
		        Component c = originalRenderer.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);

		        if (c instanceof JComponent && value != null) {
		            ((JComponent) c).setToolTipText(value.toString());
		        }

		        return c;
		    }
		});
	}
	
	public void setDefaults() {
		menu_interface.setSelected(true);
		menu_materialname.setSelected(true);
		menu_tooltip.setSelected(false);
		menu_text_bg.setSelected(true);
		menu_elem_colors.setSelected(true);
		menu_borders.setSelected(true);
		menu_carriers.setSelected(false);
		menu_probes.setSelected(true);
		menu_time.setSelected(true);
		gui_carriers.setSelected(false);
		gui_brush.setSelectedItem(Brush.INTERACT);
	}
	
	public void initialize(Simulation e) {
		getContentPane().add(e.canvas, BorderLayout.CENTER);
		e.canvas.addMouseListener(e.controls);
		e.canvas.addMouseMotionListener(e.controls);
		e.canvas.addMouseWheelListener(e.controls);
		e.canvas.addKeyListener(e.controls);
		gui_reset.addActionListener(e.controls);
		gui_brush.addActionListener(e.controls);
		gui_material.addActionListener(e.controls);
		gui_carriers.addActionListener(e.controls);
		menu_advancedsettings.addActionListener(e.controls);
		
		menu_open.addActionListener(e.controls);
		menu_save.addActionListener(e.controls);
		menu_about.addActionListener(e.controls);
		menu_help.addActionListener(e.controls);
		menu_undo.addActionListener(e.controls);
		menu_redo.addActionListener(e.controls);
		menu_save.addActionListener(e.controls);
		menu_cut.addActionListener(e.controls);
		menu_copy.addActionListener(e.controls);
		menu_paste.addActionListener(e.controls);
		menu_editdesc.addActionListener(e.controls);
		menu_new.addActionListener(e.controls);
		menu_rotate.addActionListener(e.controls);
		menu_flip_v.addActionListener(e.controls);
		menu_flip_h.addActionListener(e.controls);
		menu_img.addActionListener(e.controls);
		menu_carriers.addActionListener(e.controls);
		menu_debug.addActionListener(e.controls);
		
		gui_brush.addItemListener(e.controls);
		
		
		e.controls.brushes.initialize(Controls.Brush.values(), menu_tools, e.controls, Controls.Brush.INTERACT, new Controls.Brush[] {Controls.Brush.DRAW, Controls.Brush.VOLTAGE, Controls.Brush.BANDS});

		e.controls.scalarview.initialize(ScalarView.values(), menu_view, e.controls, ScalarView.CHARGE, null);
		menu_view.add(new JSeparator());
		e.controls.scalarmode.initialize(ScalarMode.values(), menu_view, e.controls, ScalarMode.COLORS, null);
		menu_view.add(new JSeparator());
		e.controls.vectorview.initialize(VectorView.values(), menu_view, e.controls, VectorView.E_FIELD, null);
		menu_view.add(new JSeparator());
		e.controls.vectormode.initialize(VectorMode.values(), menu_view, e.controls, VectorMode.ARROWS, null);

		e.controls.brushes.buttonmap.get(Controls.Brush.INTERACT).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_1, 0));
		e.controls.brushes.buttonmap.get(Controls.Brush.DRAW).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_2, 0));
		e.controls.brushes.buttonmap.get(Controls.Brush.LINE).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_3, 0));
		e.controls.brushes.buttonmap.get(Controls.Brush.FILL).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_4, 0));
		e.controls.brushes.buttonmap.get(Controls.Brush.SELECT).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_5, 0));

		e.controls.scalarmode.buttonmap.get(ScalarMode.NONE).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, 0));
		e.controls.vectormode.buttonmap.get(VectorMode.NONE).setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, 0));
		
		MenuBuilder.addDirectoryToMenu(menu_examples, new File("examples"), e.savemanager.fileextension, (File f) -> e.savemanager.readfile(f));

		//gui_material.removeItem(MaterialType.ABSORBER);
		e.controls.scalarview.removeOption(ScalarView.DEBUG);
		e.controls.scalarview.removeOption(ScalarView.NONE);
		e.controls.vectorview.removeOption(VectorView.NONE);

		gui_parameter1.setEnabled(true);
		gui_parameter1.setVisible(true);
		gui_parameter1_text.setEnabled(true);
		gui_parameter1_text.setVisible(true);

		pack();

		e.controls.addKeyBinds(e.canvas);
		e.controls.addKeyBinds(panel);
		
		InputMap im = (InputMap)UIManager.get("Button.focusInputMap");
		im.put(KeyStroke.getKeyStroke("pressed SPACE"), "none");
		im.put(KeyStroke.getKeyStroke("released SPACE"), "none");
	}
}