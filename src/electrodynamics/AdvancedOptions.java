// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.JButton;

public class AdvancedOptions extends JFrame {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
	public JTextField width;
	public JTextField resolution;
	public JTextField depth;
	public JTextField mu_electron;
	public JTextField mu_hole;
	public JTextField ni_semi;
	public JTextField W_semi;
	public JTextField E_b_semi;
	public JTextField ni_metal;
	public JTextField W_metal;
	public JTextField E_b_metal;
	public JTextField W_metal_high;
	public JTextField W_metal_low;
	public JTextField recomb_rate_semi;
	public JTextField recomb_rate_metal;
	public JTextField T;
	public JButton btn_apply;
	public JButton btn_cancel;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					AdvancedOptions frame = new AdvancedOptions();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public AdvancedOptions() {
		setTitle("Advanced settings");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 615, 401);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JLabel lblNewLabel = new JLabel("width [m]");
		lblNewLabel.setToolTipText("Width of simulation domain");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(6, 11, 168, 16);
		contentPane.add(lblNewLabel);
		
		width = new JTextField();
		width.setBounds(186, 6, 98, 26);
		contentPane.add(width);
		width.setColumns(10);
		
		btn_apply = new JButton("Apply changes");
		btn_apply.setBounds(310, 321, 144, 29);
		contentPane.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(452, 321, 144, 29);
		contentPane.add(btn_cancel);
		
		JLabel lblNewLabel_1 = new JLabel("resolution");
		lblNewLabel_1.setToolTipText("Number of grid points in x or y direction. Must be power of 2");
		lblNewLabel_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1.setBounds(6, 44, 168, 16);
		contentPane.add(lblNewLabel_1);
		
		resolution = new JTextField();
		resolution.setColumns(10);
		resolution.setBounds(186, 39, 98, 26);
		contentPane.add(resolution);
		
		JLabel lblNewLabel_2 = new JLabel("depth [m]");
		lblNewLabel_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2.setBounds(6, 77, 168, 16);
		contentPane.add(lblNewLabel_2);
		
		depth = new JTextField();
		depth.setColumns(10);
		depth.setBounds(186, 72, 98, 26);
		contentPane.add(depth);
		
		JLabel lblNewLabel_3 = new JLabel("mu_electron [m^2/(V s)]");
		lblNewLabel_3.setToolTipText("Electron mobility");
		lblNewLabel_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3.setBounds(6, 110, 168, 16);
		contentPane.add(lblNewLabel_3);
		
		mu_electron = new JTextField();
		mu_electron.setColumns(10);
		mu_electron.setBounds(186, 105, 98, 26);
		contentPane.add(mu_electron);
		
		JLabel lblNewLabel_4 = new JLabel("mu_hole [m^2/(V s)]");
		lblNewLabel_4.setToolTipText("Hole mobility");
		lblNewLabel_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4.setBounds(6, 143, 168, 16);
		contentPane.add(lblNewLabel_4);
		
		mu_hole = new JTextField();
		mu_hole.setColumns(10);
		mu_hole.setBounds(186, 138, 98, 26);
		contentPane.add(mu_hole);
		
		JLabel lblNewLabel_5 = new JLabel("ni_semi [1/m^3]");
		lblNewLabel_5.setToolTipText("Semiconductor equilibrium carrier concentration");
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5.setBounds(6, 176, 168, 16);
		contentPane.add(lblNewLabel_5);
		
		ni_semi = new JTextField();
		ni_semi.setColumns(10);
		ni_semi.setBounds(186, 171, 98, 26);
		contentPane.add(ni_semi);
		
		JLabel lblNewLabel_6 = new JLabel("W_semi [eV]");
		lblNewLabel_6.setToolTipText("Semiconductor work function");
		lblNewLabel_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6.setBounds(6, 209, 168, 16);
		contentPane.add(lblNewLabel_6);
		
		W_semi = new JTextField();
		W_semi.setColumns(10);
		W_semi.setBounds(186, 204, 98, 26);
		contentPane.add(W_semi);
		
		JLabel lblEbsemi = new JLabel("E_b_semi [eV]");
		lblEbsemi.setToolTipText("Semiconductor band gap");
		lblEbsemi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblEbsemi.setBounds(6, 242, 168, 16);
		contentPane.add(lblEbsemi);
		
		E_b_semi = new JTextField();
		E_b_semi.setColumns(10);
		E_b_semi.setBounds(186, 237, 98, 26);
		contentPane.add(E_b_semi);
		
		JLabel lblNimetal = new JLabel("ni_metal [1/m^3]");
		lblNimetal.setToolTipText("Metal equilibrium carrier concentration");
		lblNimetal.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNimetal.setBounds(310, 11, 176, 16);
		contentPane.add(lblNimetal);
		
		ni_metal = new JTextField();
		ni_metal.setColumns(10);
		ni_metal.setBounds(498, 6, 98, 26);
		contentPane.add(ni_metal);
		
		JLabel lblNewLabel_1_1 = new JLabel("W_metal [eV]");
		lblNewLabel_1_1.setToolTipText("Metal work function");
		lblNewLabel_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1_1.setBounds(310, 44, 176, 16);
		contentPane.add(lblNewLabel_1_1);
		
		W_metal = new JTextField();
		W_metal.setColumns(10);
		W_metal.setBounds(498, 39, 98, 26);
		contentPane.add(W_metal);
		
		JLabel lblNewLabel_2_1 = new JLabel("E_b_metal [eV]");
		lblNewLabel_2_1.setToolTipText("Metal \"band gap\"");
		lblNewLabel_2_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2_1.setBounds(310, 77, 176, 16);
		contentPane.add(lblNewLabel_2_1);
		
		E_b_metal = new JTextField();
		E_b_metal.setColumns(10);
		E_b_metal.setBounds(498, 72, 98, 26);
		contentPane.add(E_b_metal);
		
		JLabel lblNewLabel_3_1 = new JLabel("W_metal_high [eV]");
		lblNewLabel_3_1.setToolTipText("Workfunction of low-workfunction metal");
		lblNewLabel_3_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3_1.setBounds(310, 110, 176, 16);
		contentPane.add(lblNewLabel_3_1);
		
		W_metal_high = new JTextField();
		W_metal_high.setColumns(10);
		W_metal_high.setBounds(498, 105, 98, 26);
		contentPane.add(W_metal_high);
		
		JLabel lblNewLabel_4_1 = new JLabel("W_metal_low [eV]");
		lblNewLabel_4_1.setToolTipText("Workfunction of high-workfunction metal");
		lblNewLabel_4_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4_1.setBounds(310, 143, 176, 16);
		contentPane.add(lblNewLabel_4_1);
		
		W_metal_low = new JTextField();
		W_metal_low.setColumns(10);
		W_metal_low.setBounds(498, 138, 98, 26);
		contentPane.add(W_metal_low);
		
		JLabel lblNewLabel_5_1 = new JLabel("recomb_rate_semi [m^3/s]");
		lblNewLabel_5_1.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_1.setBounds(310, 176, 176, 16);
		contentPane.add(lblNewLabel_5_1);
		
		recomb_rate_semi = new JTextField();
		recomb_rate_semi.setColumns(10);
		recomb_rate_semi.setBounds(498, 171, 98, 26);
		contentPane.add(recomb_rate_semi);
		
		JLabel lblNewLabel_6_1 = new JLabel("recomb_rate_metal [m^3/s]");
		lblNewLabel_6_1.setToolTipText("Metal carrier recombination rate (in number density units)");
		lblNewLabel_6_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6_1.setBounds(310, 209, 176, 16);
		contentPane.add(lblNewLabel_6_1);
		
		recomb_rate_metal = new JTextField();
		recomb_rate_metal.setColumns(10);
		recomb_rate_metal.setBounds(498, 204, 98, 26);
		contentPane.add(recomb_rate_metal);
		
		JLabel lblT = new JLabel("T [K]");
		lblT.setToolTipText("Global temperature");
		lblT.setHorizontalAlignment(SwingConstants.TRAILING);
		lblT.setBounds(310, 242, 176, 16);
		contentPane.add(lblT);
		
		T = new JTextField();
		T.setColumns(10);
		T.setBounds(498, 237, 98, 26);
		contentPane.add(T);
	}
}
