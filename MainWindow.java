package electrodynamics;
import java.awt.BorderLayout;
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
import javax.swing.JLabel;
import javax.swing.DefaultComboBoxModel;
import java.awt.Dimension;
import javax.swing.JToggleButton;
import javax.swing.JTextPane;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;
import java.awt.Insets;
import java.awt.Font;

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
	public JComboBox gui_brush;
	public JScrollBar gui_brushsize;
	public JScrollBar gui_brushintensity2;
	public JCheckBox gui_paused;
	public JButton gui_resetall;
	public JButton gui_save;
	public JButton gui_open;
	public JPanel panel;
	public JLabel gui_label_direction;
	public JScrollBar gui_direction;
	public JLabel lblVectorBrightness;
	public JScrollBar gui_brightness_1;
	public JScrollBar gui_simspeed_2;
	public JComboBox gui_brush_1;
	public JScrollBar gui_emf_freq;
	public JCheckBox gui_tooltip;
	public JTextArea textPane;
	public JScrollPane scrollPane;
	public JButton gui_editdesc;
	public JButton gui_help;
	public JLabel gui_brushintensity;
	public JLabel gui_acfreqlbl;
	public JLabel gui_stepslbl;
	public JLabel gui_stepsizelbl;

	/**
	 * Create the frame.
	 */
	public MainWindow() {
		setTitle("Electrodynamics");
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
		
		gui_reset = new JButton("Clear fields");
		gui_reset.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
			}
		});
		gui_reset.setBounds(10, 70, 171, 23);
		panel.add(gui_reset);
		
		gui_paused = new JCheckBox("Paused");
		gui_paused.setSelected(true);
		gui_paused.setBounds(10, 11, 68, 23);
		panel.add(gui_paused);
		
		gui_brushsize = new JScrollBar();
		gui_brushsize.setMinimum(-15);
		gui_brushsize.setMaximum(25);
		gui_brushsize.setOrientation(JScrollBar.HORIZONTAL);
		gui_brushsize.setBounds(201, 238, 171, 17);
		panel.add(gui_brushsize);
		
		gui_simspeed = new JScrollBar();
		gui_simspeed.setValue(10);
		gui_simspeed.setBlockIncrement(1);
		gui_simspeed.setMaximum(30);
		gui_simspeed.setOrientation(JScrollBar.HORIZONTAL);
		gui_simspeed.setBounds(10, 238, 171, 17);
		panel.add(gui_simspeed);
		
		gui_brightness = new JScrollBar();
		gui_brightness.setBlockIncrement(1);
		gui_brightness.setMinimum(-15);
		gui_brightness.setMaximum(25);
		gui_brightness.setOrientation(JScrollBar.HORIZONTAL);
		gui_brightness.setBounds(10, 344, 171, 17);
		panel.add(gui_brightness);
		
		gui_brushintensity2 = new JScrollBar();
		gui_brushintensity2.setMaximum(25);
		gui_brushintensity2.setMinimum(-15);
		gui_brushintensity2.setOrientation(JScrollBar.HORIZONTAL);
		gui_brushintensity2.setBounds(201, 291, 171, 17);
		panel.add(gui_brushintensity2);
		
		gui_brush = new JComboBox();
		gui_brush.setModel(new DefaultComboBoxModel(new String[] {"Mouse: Interact", "Mouse: Set charge density", "Mouse: Add conductor", "Mouse: Add battery", "Mouse: Add dielectric", "Mouse: Add paramagnet/diamagnet", "Mouse: Add switch"}));
		gui_brush.setSelectedIndex(1);
		gui_brush.setBounds(201, 132, 171, 22);
		panel.add(gui_brush);
		
		gui_view = new JComboBox();
		gui_view.setModel(new DefaultComboBoxModel(new String[] {"View: None", "View: Ex/Ey (Electric field)", "View: B (Magetic B field)", "View: \u03C1 (Charge)", "View: Jx/Jy (Current)", "View: H (Magnetic H field)", "View: u (Energy density)"}));
		gui_view.setSelectedIndex(3);
		gui_view.setBounds(10, 132, 171, 22);
		panel.add(gui_view);
		
		gui_stepsizelbl = new JLabel("Step size");
		gui_stepsizelbl.setBounds(20, 213, 161, 14);
		panel.add(gui_stepsizelbl);
		
		JLabel label5 = new JLabel("Scalar Brightness");
		label5.setBounds(20, 319, 150, 14);
		panel.add(label5);
		
		JLabel lblBrushSize = new JLabel("Brush size");
		lblBrushSize.setBounds(211, 213, 87, 14);
		panel.add(lblBrushSize);
		
		gui_brushintensity = new JLabel("Brush intensity");
		gui_brushintensity.setBounds(211, 266, 154, 14);
		panel.add(gui_brushintensity);
		
		gui_view_vec = new JComboBox();
		gui_view_vec.setModel(new DefaultComboBoxModel(new String[] {"View: None", "View: E (Electric field)", "View: J (Current)", "View: D (Displacement field)", "View: S (Poynting vector)"}));
		gui_view_vec.setSelectedIndex(1);
		gui_view_vec.setBounds(10, 165, 171, 22);
		panel.add(gui_view_vec);
		
		gui_save = new JButton("Save to file");
		gui_save.setBounds(201, 11, 171, 23);
		panel.add(gui_save);
		
		gui_open = new JButton("Load scenario");
		gui_open.setBounds(201, 40, 171, 23);
		panel.add(gui_open);
		
		gui_resetall = new JButton("Reset all");
		gui_resetall.setBounds(201, 70, 171, 23);
		panel.add(gui_resetall);
		
		gui_label_direction = new JLabel("Direction");
		gui_label_direction.setBounds(211, 323, 161, 14);
		panel.add(gui_label_direction);
		
		gui_direction = new JScrollBar();
		gui_direction.setOrientation(JScrollBar.HORIZONTAL);
		gui_direction.setMaximum(34);
		gui_direction.setBounds(201, 344, 171, 17);
		panel.add(gui_direction);
		
		lblVectorBrightness = new JLabel("Vector field brightness");
		lblVectorBrightness.setBounds(20, 376, 150, 14);
		panel.add(lblVectorBrightness);
		
		gui_brightness_1 = new JScrollBar();
		gui_brightness_1.setOrientation(JScrollBar.HORIZONTAL);
		gui_brightness_1.setMinimum(-15);
		gui_brightness_1.setMaximum(25);
		gui_brightness_1.setBlockIncrement(1);
		gui_brightness_1.setBounds(10, 399, 171, 17);
		panel.add(gui_brightness_1);
		
		gui_brush_1 = new JComboBox();
		gui_brush_1.setModel(new DefaultComboBoxModel(new String[] {"Shape: Circle", "Shape: Square"}));
		gui_brush_1.setSelectedIndex(0);
		gui_brush_1.setBounds(201, 165, 171, 22);
		panel.add(gui_brush_1);
		
		gui_stepslbl = new JLabel("Steps/Frame");
		gui_stepslbl.setBounds(20, 266, 144, 14);
		panel.add(gui_stepslbl);
		
		gui_simspeed_2 = new JScrollBar();
		gui_simspeed_2.setMinimum(1);
		gui_simspeed_2.setValue(10);
		gui_simspeed_2.setOrientation(JScrollBar.HORIZONTAL);
		gui_simspeed_2.setMaximum(30);
		gui_simspeed_2.setBlockIncrement(1);
		gui_simspeed_2.setBounds(10, 291, 171, 17);
		panel.add(gui_simspeed_2);
		
		gui_emf_freq = new JScrollBar();
		gui_emf_freq.setOrientation(JScrollBar.HORIZONTAL);
		gui_emf_freq.setMaximum(30);
		gui_emf_freq.setBlockIncrement(1);
		gui_emf_freq.setBounds(201, 399, 171, 17);
		panel.add(gui_emf_freq);
		
		gui_acfreqlbl = new JLabel("AC Frequency");
		gui_acfreqlbl.setBounds(211, 376, 150, 14);
		panel.add(gui_acfreqlbl);
		
		gui_tooltip = new JCheckBox("Show tooltip");
		gui_tooltip.setSelected(true);
		gui_tooltip.setBounds(10, 40, 101, 23);
		panel.add(gui_tooltip);
		
		scrollPane = new JScrollPane();
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.setBounds(20, 434, 345, 218);
		panel.add(scrollPane);
		
		textPane = new JTextArea();
		textPane.setEnabled(false);
		textPane.setText("Description of simulation scenario");
		textPane.setFont(new Font("SansSerif", Font.PLAIN, 13));
		textPane.setMargin(new Insets(4, 4, 4, 4));
		textPane.setLineWrap(true);
		scrollPane.setViewportView(textPane);
		
		gui_editdesc = new JButton("Edit description");
		gui_editdesc.setBounds(197, 661, 165, 23);
		panel.add(gui_editdesc);
		
		gui_help = new JButton("Help/About");
		gui_help.setBounds(24, 661, 162, 23);
		panel.add(gui_help);
	}
}
