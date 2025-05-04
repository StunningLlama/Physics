package electrodynamics;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JScrollBar;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.DefaultComboBoxModel;
import java.awt.Dimension;
import javax.swing.JToggleButton;
import javax.swing.ListCellRenderer;
import javax.swing.JTextPane;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;
import java.awt.Insets;
import java.awt.Font;
import electrodynamics.Electrodynamics.Brush;
import electrodynamics.Electrodynamics.ScalarView;
import electrodynamics.Electrodynamics.VectorView;
import electrodynamics.Electrodynamics.BrushShape;
import electrodynamics.Electrodynamics.VectorMode;
import electrodynamics.Electrodynamics.BoundaryCondition;

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
	public JButton gui_resetall;
	public JButton gui_save;
	public JButton gui_open;
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
	public JButton gui_editdesc;
	public JButton gui_help;
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

	/**
	 * Create the frame.
	 */
	public MainWindow() {
		setTitle("Brandon's semiconductor simulator");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 589, 746);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		
		panel = new JPanel();
		panel.setPreferredSize(new Dimension(380, 200));
		panel.setMinimumSize(new Dimension(200, 200));
		contentPane.add(panel, BorderLayout.EAST);
		panel.setLayout(null);
		
		gui_reset = new JButton("Set fields to zero");
		gui_reset.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
			}
		});
		gui_reset.setBounds(201, 99, 171, 23);
		panel.add(gui_reset);
		
		gui_paused = new JCheckBox("Paused");
		gui_paused.setSelected(true);
		gui_paused.setBounds(10, 11, 101, 23);
		panel.add(gui_paused);
		
		gui_brushsize = new JScrollBar();
		gui_brushsize.setMaximum(1010);
		gui_brushsize.setValue(500);
		gui_brushsize.setOrientation(JScrollBar.HORIZONTAL);
		gui_brushsize.setBounds(201, 321, 171, 17);
		panel.add(gui_brushsize);
		
		gui_simspeed = new JScrollBar();
		gui_simspeed.setValue(20);
		gui_simspeed.setBlockIncrement(1);
		gui_simspeed.setMaximum(30);
		gui_simspeed.setOrientation(JScrollBar.HORIZONTAL);
		gui_simspeed.setBounds(10, 321, 171, 17);
		panel.add(gui_simspeed);
		
		gui_brightness = new JScrollBar();
		gui_brightness.setValue(-20);
		gui_brightness.setBlockIncrement(1);
		gui_brightness.setMinimum(-45);
		gui_brightness.setMaximum(45);
		gui_brightness.setOrientation(JScrollBar.HORIZONTAL);
		gui_brightness.setBounds(10, 427, 171, 17);
		panel.add(gui_brightness);
		
		gui_parameter1 = new JScrollBar();
		gui_parameter1.setVisible(false);
		gui_parameter1.setEnabled(false);
		gui_parameter1.setMaximum(25);
		gui_parameter1.setMinimum(-15);
		gui_parameter1.setOrientation(JScrollBar.HORIZONTAL);
		gui_parameter1.setBounds(201, 374, 171, 17);
		panel.add(gui_parameter1);
		
		gui_brush = new JComboBox();
		gui_brush.setMaximumRowCount(16);
		gui_brush.setModel(new DefaultComboBoxModel(Brush.values()));
		gui_brush.setSelectedIndex(0);
		gui_brush.setBounds(201, 151, 171, 22);
		panel.add(gui_brush);
		addTooltips(gui_brush);
		
		gui_view = new JComboBox();
		gui_view.setMaximumRowCount(16);
		gui_view.setToolTipText("");
		gui_view.setModel(new DefaultComboBoxModel(ScalarView.values()));
		gui_view.setSelectedIndex(3);
		gui_view.setBounds(10, 151, 171, 22);
		panel.add(gui_view);
		addTooltips(gui_view);
		
		gui_stepsizelbl = new JLabel("Step size");
		gui_stepsizelbl.setBounds(20, 296, 161, 14);
		panel.add(gui_stepsizelbl);
		
		JLabel label5 = new JLabel("Scalar Brightness");
		label5.setBounds(20, 402, 150, 14);
		panel.add(label5);
		
		lblBrushSize = new JLabel("Brush size");
		lblBrushSize.setBounds(211, 296, 138, 14);
		panel.add(lblBrushSize);
		
		gui_parameter1_text = new JLabel("");
		gui_parameter1_text.setEnabled(false);
		gui_parameter1_text.setBounds(211, 349, 154, 14);
		panel.add(gui_parameter1_text);
		
		gui_view_vec = new JComboBox();
		gui_view_vec.setMaximumRowCount(16);
		gui_view_vec.setModel(new DefaultComboBoxModel(VectorView.values()));
		gui_view_vec.setSelectedIndex(1);
		gui_view_vec.setBounds(10, 185, 171, 22);
		panel.add(gui_view_vec);
		addTooltips(gui_view_vec);
		
		gui_save = new JButton("Save scenario");
		gui_save.setBounds(201, 11, 171, 23);
		panel.add(gui_save);
		
		gui_open = new JButton("Load scenario");
		gui_open.setBounds(201, 40, 171, 23);
		panel.add(gui_open);
		
		gui_resetall = new JButton("Clear all");
		gui_resetall.setBounds(201, 70, 171, 23);
		panel.add(gui_resetall);
		
		gui_parameter2_text = new JLabel("Direction");
		gui_parameter2_text.setBounds(211, 406, 161, 14);
		panel.add(gui_parameter2_text);
		
		gui_parameter2 = new JScrollBar();
		gui_parameter2.setOrientation(JScrollBar.HORIZONTAL);
		gui_parameter2.setMaximum(34);
		gui_parameter2.setBounds(201, 427, 171, 17);
		panel.add(gui_parameter2);
		
		lblVectorBrightness = new JLabel("Vector field brightness");
		lblVectorBrightness.setBounds(20, 459, 150, 14);
		panel.add(lblVectorBrightness);
		
		gui_brightness_vec = new JScrollBar();
		gui_brightness_vec.setValue(-10);
		gui_brightness_vec.setOrientation(JScrollBar.HORIZONTAL);
		gui_brightness_vec.setMinimum(-45);
		gui_brightness_vec.setMaximum(45);
		gui_brightness_vec.setBlockIncrement(1);
		gui_brightness_vec.setBounds(10, 482, 171, 17);
		panel.add(gui_brightness_vec);
		
		gui_brush_1 = new JComboBox();
		gui_brush_1.setMaximumRowCount(16);
		gui_brush_1.setModel(new DefaultComboBoxModel(BrushShape.values()));
		gui_brush_1.setSelectedIndex(1);
		gui_brush_1.setBounds(201, 219, 171, 22);
		panel.add(gui_brush_1);
		addTooltips(gui_brush_1);
		
		gui_stepslbl = new JLabel("Steps/Frame");
		gui_stepslbl.setBounds(20, 349, 144, 14);
		panel.add(gui_stepslbl);
		
		gui_simspeed_2 = new JScrollBar();
		gui_simspeed_2.setMinimum(1);
		gui_simspeed_2.setValue(25);
		gui_simspeed_2.setOrientation(JScrollBar.HORIZONTAL);
		gui_simspeed_2.setMaximum(110);
		gui_simspeed_2.setBlockIncrement(1);
		gui_simspeed_2.setBounds(10, 374, 171, 17);
		panel.add(gui_simspeed_2);
		
		gui_parameter3 = new JScrollBar();
		gui_parameter3.setMinimum(-50);
		gui_parameter3.setOrientation(JScrollBar.HORIZONTAL);
		gui_parameter3.setMaximum(60);
		gui_parameter3.setBlockIncrement(1);
		gui_parameter3.setBounds(201, 482, 171, 17);
		panel.add(gui_parameter3);
		
		gui_parameter3_text = new JLabel("EMF");
		gui_parameter3_text.setBounds(211, 459, 150, 14);
		panel.add(gui_parameter3_text);
		
		gui_tooltip = new JCheckBox("Show physical quantities");
		gui_tooltip.setSelected(true);
		gui_tooltip.setBounds(10, 40, 188, 23);
		panel.add(gui_tooltip);
		
		scrollPane = new JScrollPane();
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.setBounds(27, 521, 345, 146);
		panel.add(scrollPane);
		
		textPane = new JTextArea();
		textPane.setWrapStyleWord(true);
		textPane.setText("Description of simulation scenario");
		textPane.setFont(new Font("SansSerif", Font.PLAIN, 13));
		textPane.setMargin(new Insets(4, 4, 4, 4));
		textPane.setLineWrap(true);
		scrollPane.setColumnHeaderView(textPane);
		
		gui_editdesc = new JButton("Edit description");
		gui_editdesc.setBounds(196, 679, 165, 23);
		panel.add(gui_editdesc);
		
		gui_help = new JButton("Help/About");
		gui_help.setBounds(19, 679, 162, 23);
		panel.add(gui_help);
		
		gui_view_vec_mode = new JComboBox();
		gui_view_vec_mode.setMaximumRowCount(16);
		gui_view_vec_mode.setModel(new DefaultComboBoxModel(VectorMode.values()));
		gui_view_vec_mode.setSelectedIndex(0);
		gui_view_vec_mode.setBounds(10, 219, 171, 22);
		panel.add(gui_view_vec_mode);
		addTooltips(gui_view_vec_mode);
		
		gui_brush_highlight = new JCheckBox("Brush highlight");
		gui_brush_highlight.setSelected(true);
		gui_brush_highlight.setBounds(205, 253, 160, 23);
		panel.add(gui_brush_highlight);
		
		gui_elem_colors = new JCheckBox("Show material colors");
		gui_elem_colors.setSelected(true);
		gui_elem_colors.setBounds(10, 253, 171, 23);
		panel.add(gui_elem_colors);
		
		gui_material = new JComboBox();
		gui_material.setMaximumRowCount(16);
		gui_material.setModel(new DefaultComboBoxModel(electrodynamics.MaterialType.values()));
		gui_material.setSelectedIndex(0);
		gui_material.setBounds(201, 185, 171, 22);
		panel.add(gui_material);
		addTooltips(gui_material);
		
		gui_bc = new JComboBox();
		gui_bc.setModel(new DefaultComboBoxModel(BoundaryCondition.values()));
		gui_bc.setSelectedIndex(0);
		gui_bc.setBounds(10, 99, 171, 22);
		panel.add(gui_bc);
		addTooltips(gui_bc);
		
		gui_text_bg = new JCheckBox("Black text background");
		gui_text_bg.setSelected(true);
		gui_text_bg.setBounds(10, 68, 179, 23);
		panel.add(gui_text_bg);
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