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
import javax.swing.border.EmptyBorder;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JMenu;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import java.awt.event.KeyEvent;
import java.awt.event.InputEvent;

public class MainWindow extends JFrame {

	/**
	 *
	 */
	private static final long serialVersionUID = -5756219569007074449L;

	public JPanel contentPane;
	public JButton gui_reset;
	public JComboBox gui_view;
	public JComboBox gui_view_vec;
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
	public JCheckBox gui_tooltip;
	public JTextArea textPane;
	public JScrollPane scrollPane;
	public JLabel gui_parameter1_text;
	public JLabel gui_parameter3_text;
	public JLabel gui_stepslbl;
	public JLabel gui_stepsizelbl;
	public JComboBox gui_view_vec_mode;
	public JCheckBox gui_brush_highlight;
	public JComboBox gui_material;
	public JComboBox gui_brush;
	public JCheckBox gui_elem_colors;
	public JComboBox gui_bc;
	public JCheckBox gui_text_bg;
	public JLabel lblBrushSize;
	public JButton gui_adv_settings;
	public JCheckBox gui_interface;
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

	/**
	 * Create the frame.
	 */
	public MainWindow() {
		setTitle("Brandon's semiconductor simulator");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 590, 832);
		
		JMenuBar menuBar = new JMenuBar();
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
		
		menu_cut = new JMenuItem("Cut");
		menu_cut.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_cut);
		
		menu_copy = new JMenuItem("Copy");
		menu_copy.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_copy);
		
		menu_paste = new JMenuItem("Paste");
		menu_paste.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_paste);
		
		JSeparator separator = new JSeparator();
		menu_asdf.add(separator);
		
		menu_undo = new JMenuItem("Undo");
		menu_undo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK));
		menu_asdf.add(menu_undo);
		
		menu_redo = new JMenuItem("Redo");
		menu_redo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, InputEvent.CTRL_DOWN_MASK | InputEvent.SHIFT_DOWN_MASK));
		menu_asdf.add(menu_redo);
		
		menu_tools = new JMenu("Tools");
		menuBar.add(menu_tools);
		
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
		gui_reset.setBounds(201, 40, 171, 23);
		panel.add(gui_reset);

		gui_paused = new JCheckBox("Paused");
		gui_paused.setSelected(false);
		gui_paused.setBounds(10, 9, 101, 23);
		panel.add(gui_paused);

		gui_brushsize = new JScrollBar();
		gui_brushsize.setMaximum(1010);
		gui_brushsize.setValue(500);
		gui_brushsize.setOrientation(Adjustable.HORIZONTAL);
		gui_brushsize.setBounds(201, 320, 171, 17);
		panel.add(gui_brushsize);

		gui_simspeed = new JScrollBar();
		gui_simspeed.setValue(20);
		gui_simspeed.setBlockIncrement(1);
		gui_simspeed.setMaximum(30);
		gui_simspeed.setOrientation(Adjustable.HORIZONTAL);
		gui_simspeed.setBounds(10, 320, 171, 17);
		panel.add(gui_simspeed);

		gui_brightness = new JScrollBar();
		gui_brightness.setValue(-20);
		gui_brightness.setBlockIncrement(1);
		gui_brightness.setMinimum(-45);
		gui_brightness.setMaximum(45);
		gui_brightness.setOrientation(Adjustable.HORIZONTAL);
		gui_brightness.setBounds(10, 426, 171, 17);
		panel.add(gui_brightness);

		gui_parameter1 = new JScrollBar();
		gui_parameter1.setVisible(false);
		gui_parameter1.setEnabled(false);
		gui_parameter1.setMaximum(25);
		gui_parameter1.setMinimum(-15);
		gui_parameter1.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter1.setBounds(201, 373, 171, 17);
		panel.add(gui_parameter1);

		gui_brush = new JComboBox();
		gui_brush.setMaximumRowCount(16);
		gui_brush.setModel(new DefaultComboBoxModel(Controls.Brush.values()));
		gui_brush.setSelectedIndex(0);
		gui_brush.setBounds(201, 150, 171, 22);
		panel.add(gui_brush);
		addTooltips(gui_brush);

		gui_view = new JComboBox();
		gui_view.setMaximumRowCount(16);
		gui_view.setToolTipText("");
		gui_view.setModel(new DefaultComboBoxModel(Renderer.ScalarView.values()));
		gui_view.setSelectedIndex(3);
		gui_view.setBounds(10, 150, 171, 22);
		panel.add(gui_view);
		addTooltips(gui_view);

		gui_stepsizelbl = new JLabel("Timestep");
		gui_stepsizelbl.setBounds(20, 295, 161, 14);
		panel.add(gui_stepsizelbl);

		JLabel label5 = new JLabel("Scalar brightness");
		label5.setBounds(20, 401, 150, 14);
		panel.add(label5);

		lblBrushSize = new JLabel("Brush size");
		lblBrushSize.setBounds(211, 295, 138, 14);
		panel.add(lblBrushSize);

		gui_parameter1_text = new JLabel("");
		gui_parameter1_text.setEnabled(false);
		gui_parameter1_text.setBounds(211, 348, 154, 14);
		panel.add(gui_parameter1_text);

		gui_view_vec = new JComboBox();
		gui_view_vec.setMaximumRowCount(16);
		gui_view_vec.setModel(new DefaultComboBoxModel(Renderer.VectorView.values()));
		gui_view_vec.setSelectedIndex(1);
		gui_view_vec.setBounds(10, 184, 171, 22);
		panel.add(gui_view_vec);
		addTooltips(gui_view_vec);

		gui_parameter2_text = new JLabel("Direction");
		gui_parameter2_text.setBounds(211, 405, 161, 14);
		panel.add(gui_parameter2_text);

		gui_parameter2 = new JScrollBar();
		gui_parameter2.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter2.setMaximum(34);
		gui_parameter2.setBounds(201, 426, 171, 17);
		panel.add(gui_parameter2);

		lblVectorBrightness = new JLabel("Vector field brightness");
		lblVectorBrightness.setBounds(20, 458, 150, 14);
		panel.add(lblVectorBrightness);

		gui_brightness_vec = new JScrollBar();
		gui_brightness_vec.setValue(-10);
		gui_brightness_vec.setOrientation(Adjustable.HORIZONTAL);
		gui_brightness_vec.setMinimum(-45);
		gui_brightness_vec.setMaximum(45);
		gui_brightness_vec.setBlockIncrement(1);
		gui_brightness_vec.setBounds(10, 481, 171, 17);
		panel.add(gui_brightness_vec);

		gui_brush_1 = new JComboBox();
		gui_brush_1.setMaximumRowCount(16);
		gui_brush_1.setModel(new DefaultComboBoxModel(Controls.BrushShape.values()));
		gui_brush_1.setSelectedIndex(1);
		gui_brush_1.setBounds(201, 218, 171, 22);
		panel.add(gui_brush_1);
		addTooltips(gui_brush_1);

		gui_stepslbl = new JLabel("Sim steps/frame");
		gui_stepslbl.setBounds(20, 348, 144, 14);
		panel.add(gui_stepslbl);

		gui_simspeed_2 = new JScrollBar();
		gui_simspeed_2.setMinimum(1);
		gui_simspeed_2.setValue(25);
		gui_simspeed_2.setOrientation(Adjustable.HORIZONTAL);
		gui_simspeed_2.setMaximum(110);
		gui_simspeed_2.setBlockIncrement(1);
		gui_simspeed_2.setBounds(10, 373, 171, 17);
		panel.add(gui_simspeed_2);

		gui_parameter3 = new JScrollBar();
		gui_parameter3.setMinimum(-50);
		gui_parameter3.setOrientation(Adjustable.HORIZONTAL);
		gui_parameter3.setMaximum(60);
		gui_parameter3.setBlockIncrement(1);
		gui_parameter3.setBounds(201, 481, 171, 17);
		panel.add(gui_parameter3);

		gui_parameter3_text = new JLabel("EMF");
		gui_parameter3_text.setBounds(211, 458, 150, 14);
		panel.add(gui_parameter3_text);

		gui_tooltip = new JCheckBox("Show detailed info");
		gui_tooltip.setSelected(true);
		gui_tooltip.setBounds(10, 38, 188, 23);
		panel.add(gui_tooltip);

		scrollPane = new JScrollPane();
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.setBounds(16, 523, 349, 232);
		panel.add(scrollPane);

		textPane = new JTextArea();
		textPane.setWrapStyleWord(true);
		textPane.setText("Description of simulation scenario");
		textPane.setFont(new Font("SansSerif", Font.PLAIN, 13));
		textPane.setMargin(new Insets(4, 4, 4, 4));
		textPane.setLineWrap(true);
		textPane.setEditable(false);
		scrollPane.setColumnHeaderView(textPane);

		gui_view_vec_mode = new JComboBox();
		gui_view_vec_mode.setMaximumRowCount(16);
		gui_view_vec_mode.setModel(new DefaultComboBoxModel(Renderer.VectorMode.values()));
		gui_view_vec_mode.setSelectedIndex(0);
		gui_view_vec_mode.setBounds(10, 218, 171, 22);
		panel.add(gui_view_vec_mode);
		addTooltips(gui_view_vec_mode);

		gui_brush_highlight = new JCheckBox("Brush highlight");
		gui_brush_highlight.setSelected(true);
		gui_brush_highlight.setBounds(205, 252, 160, 23);
		panel.add(gui_brush_highlight);

		gui_elem_colors = new JCheckBox("Show material colors");
		gui_elem_colors.setSelected(true);
		gui_elem_colors.setBounds(10, 252, 171, 23);
		panel.add(gui_elem_colors);

		gui_material = new JComboBox();
		gui_material.setMaximumRowCount(16);
		gui_material.setModel(new DefaultComboBoxModel(electrodynamics.MaterialType.values()));
		gui_material.setSelectedIndex(0);
		gui_material.setBounds(201, 184, 171, 22);
		panel.add(gui_material);
		addTooltips(gui_material);

		gui_bc = new JComboBox();
		gui_bc.setModel(new DefaultComboBoxModel(BoundaryCondition.values()));
		gui_bc.setSelectedIndex(0);
		gui_bc.setBounds(201, 70, 171, 22);
		panel.add(gui_bc);
		addTooltips(gui_bc);

		gui_text_bg = new JCheckBox("Show text background");
		gui_text_bg.setSelected(true);
		gui_text_bg.setBounds(10, 68, 179, 23);
		panel.add(gui_text_bg);
		
		gui_adv_settings = new JButton("Advanced settings");
		gui_adv_settings.setBounds(201, 11, 171, 23);
		panel.add(gui_adv_settings);
		
		gui_interface = new JCheckBox("Display interface");
		gui_interface.setSelected(true);
		gui_interface.setBounds(10, 97, 179, 23);
		panel.add(gui_interface);
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
}