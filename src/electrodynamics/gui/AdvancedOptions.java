// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.JsonObject;
import com.google.gson.JsonSyntaxException;
import com.google.gson.stream.JsonReader;

import electrodynamics.Preset;
import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.units.Units;
import electrodynamics.util.Utils;

import javax.swing.JTabbedPane;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.io.StringReader;

public class AdvancedOptions extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	private JPanel sim;
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
	public JTextField k_rad_semi;
	public JTextField recomb_rate_metal;
	public JTextField T;
	public JButton btn_apply;
	public JButton btn_cancel;
	private JTabbedPane tabbedPane;
	private JPanel other;
	private JLabel lblNewLabel_5_2;
	private JTextField n_default_doping;
	private JLabel lblNewLabel_5_3;
	private JTextField p_default_doping;
	private JLabel lblNewLabel_5_4;
	private JTextField n_light_doping;
	private JLabel lblNewLabel_5_5;
	private JTextField p_light_doping;
	private JLabel lblNewLabel_5_6;
	private JTextField n_heavy_doping;
	private JLabel lblNewLabel_5_7;
	private JTextField p_heavy_doping;
	private JLabel lblVacuumPermittivitysi;
	private JTextField eps0;
	private JLabel lblVacuumPermeabilitysi;
	private JTextField mu0;
	private JLabel lblBoltzmannConstantk;
	private JTextField e_charge;
	private JLabel lblNewLabel_7;
	private JTextField junction_size;
	private JLabel lblDielectricRelPermittivity;
	private JTextField dielectric_eps_r;
	private JLabel lblFerromagnetRelPermeability;
	private JTextField ferromagnet_mu_r;
	private JLabel lblStaticChargeDensity;
	private JTextField staticcharge_density;
	private JLabel lblCurrentSourceMobilit;
	private JTextField currentsource_mobility;
	private JPanel panel;
	public JButton btn_reset;
	private JLabel lblNewLabel_5_8;
	private JTextField v_sat_n;
	private JLabel lblNewLabel_5_9;
	private JTextField v_sat_p;
	private JLabel lblNewLabel_5_10;
	private JTextField k_SRH_n_semi;
	private JLabel lblNewLabel_5_11;
	private JTextField k_SRH_p_semi;
	private JLabel lblNewLabel_5_12;
	private JTextField k_aug_n_semi;
	private JLabel lblNewLabel_5_13;
	private JTextField k_aug_p_semi;
	private JButton btn_presets;
	private JLabel lblNewLabel_8;
	private JTextField dopant_smoothing_distance;

	public AdvancedOptions() {
		setTitle("Advanced settings");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 699, 447);

		panel = new JPanel();
		panel.setLayout(null);

		tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		tabbedPane.setBounds(0, 0, 699, 372);
		
		sim = new JPanel();
		tabbedPane.addTab("Simulation", null, sim, null);
		sim.setLayout(null);
		
		JPanel phys = new JPanel();
		tabbedPane.addTab("Physics", null, phys, null);
		phys.setLayout(null);
		
		JPanel metal = new JPanel();
		tabbedPane.addTab("Metal", null, metal, null);
		metal.setLayout(null);
		
		panel.add(tabbedPane);

		setContentPane(panel);
		
		JLabel lblNewLabel = new JLabel("Sim. region width [m]");
		lblNewLabel.setToolTipText("Width of simulation domain");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(45, 11, 168, 16);
		sim.add(lblNewLabel);
		
		width = new JTextField();
		width.setBounds(225, 6, 98, 26);
		sim.add(width);
		width.setColumns(10);
		
		JLabel lblNewLabel_1 = new JLabel("Resolution (Num. grid points)");
		lblNewLabel_1.setToolTipText("Number of grid points in x or y direction. Must be power of 2");
		lblNewLabel_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1.setBounds(27, 44, 186, 16);
		sim.add(lblNewLabel_1);
		
		resolution = new JTextField();
		resolution.setColumns(10);
		resolution.setBounds(225, 39, 98, 26);
		sim.add(resolution);
		
		JLabel lblNewLabel_2 = new JLabel("Depth [m]");
		lblNewLabel_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2.setBounds(45, 77, 168, 16);
		sim.add(lblNewLabel_2);
		
		depth = new JTextField();
		depth.setColumns(10);
		depth.setBounds(225, 72, 98, 26);
		sim.add(depth);
		
		lblNewLabel_7 = new JLabel("Junction smoothing [px]");
		lblNewLabel_7.setToolTipText("Smooths junction between materials with different chemical potential (ie. metal-semiconductor junctions)");
		lblNewLabel_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_7.setBounds(45, 110, 168, 16);
		sim.add(lblNewLabel_7);
		
		junction_size = new JTextField();
		junction_size.setColumns(10);
		junction_size.setBounds(225, 105, 98, 26);
		sim.add(junction_size);
		
		lblNewLabel_8 = new JLabel("Dopant smoothing [px]");
		lblNewLabel_8.setToolTipText("Smooths dopant density, use when modelling non abrupt PN junctions.");
		lblNewLabel_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_8.setBounds(45, 143, 168, 16);
		sim.add(lblNewLabel_8);
		
		dopant_smoothing_distance = new JTextField();
		dopant_smoothing_distance.setColumns(10);
		dopant_smoothing_distance.setBounds(225, 138, 98, 26);
		sim.add(dopant_smoothing_distance);
		
		JLabel lblNimetal = new JLabel("Intrinsic carrier conc. [1/m^3]");
		lblNimetal.setToolTipText("Metal equilibrium carrier concentration");
		lblNimetal.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNimetal.setBounds(6, 11, 203, 16);
		metal.add(lblNimetal);
		
		ni_metal = new JTextField();
		ni_metal.setColumns(10);
		ni_metal.setBounds(221, 6, 98, 26);
		metal.add(ni_metal);
		
		JLabel lblNewLabel_1_1 = new JLabel("Workfunction (default) [eV]");
		lblNewLabel_1_1.setToolTipText("");
		lblNewLabel_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1_1.setBounds(33, 77, 176, 16);
		metal.add(lblNewLabel_1_1);
		
		W_metal = new JTextField();
		W_metal.setColumns(10);
		W_metal.setBounds(221, 72, 98, 26);
		metal.add(W_metal);
		
		JLabel lblNewLabel_2_1 = new JLabel("Metal \"bandgap\" [eV]");
		lblNewLabel_2_1.setToolTipText("");
		lblNewLabel_2_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2_1.setBounds(33, 44, 176, 16);
		metal.add(lblNewLabel_2_1);
		
		E_b_metal = new JTextField();
		E_b_metal.setColumns(10);
		E_b_metal.setBounds(221, 39, 98, 26);
		metal.add(E_b_metal);
		
		JLabel lblNewLabel_3_1 = new JLabel("Workfunction (High WF metal) [eV]");
		lblNewLabel_3_1.setToolTipText("");
		lblNewLabel_3_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3_1.setBounds(330, 11, 232, 16);
		metal.add(lblNewLabel_3_1);
		
		W_metal_high = new JTextField();
		W_metal_high.setColumns(10);
		W_metal_high.setBounds(574, 6, 98, 26);
		metal.add(W_metal_high);
		
		JLabel lblNewLabel_4_1 = new JLabel("Workfunction (Low WF metal) [eV]");
		lblNewLabel_4_1.setToolTipText("");
		lblNewLabel_4_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4_1.setBounds(340, 44, 222, 16);
		metal.add(lblNewLabel_4_1);
		
		W_metal_low = new JTextField();
		W_metal_low.setColumns(10);
		W_metal_low.setBounds(574, 39, 98, 26);
		metal.add(W_metal_low);
		
		JLabel lblNewLabel_6_1 = new JLabel("Recomb. rate [m^3/s]");
		lblNewLabel_6_1.setToolTipText("Radiative recombination rate constant in metal");
		lblNewLabel_6_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6_1.setBounds(39, 110, 170, 16);
		metal.add(lblNewLabel_6_1);
		
		recomb_rate_metal = new JTextField();
		recomb_rate_metal.setColumns(10);
		recomb_rate_metal.setBounds(221, 105, 98, 26);
		metal.add(recomb_rate_metal);
		
		JLabel lblCarrierConchigh = new JLabel("Carrier conc. (High cond.) [1/m^3]");
		lblCarrierConchigh.setToolTipText("Metal equilibrium carrier concentration");
		lblCarrierConchigh.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCarrierConchigh.setBounds(331, 77, 231, 16);
		metal.add(lblCarrierConchigh);
		
		ni_metal_high = new JTextField();
		ni_metal_high.setColumns(10);
		ni_metal_high.setBounds(574, 72, 98, 26);
		metal.add(ni_metal_high);
		
		JLabel lblCarrierConclow = new JLabel("Carrier conc. (Low cond.) [1/m^3]");
		lblCarrierConclow.setToolTipText("Metal equilibrium carrier concentration");
		lblCarrierConclow.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCarrierConclow.setBounds(330, 110, 231, 16);
		metal.add(lblCarrierConclow);
		
		ni_metal_low = new JTextField();
		ni_metal_low.setColumns(10);
		ni_metal_low.setBounds(573, 105, 98, 26);
		metal.add(ni_metal_low);
		
		
		JLabel lblT = new JLabel("Temperature [K]");
		lblT.setToolTipText("Global temperature");
		lblT.setHorizontalAlignment(SwingConstants.TRAILING);
		lblT.setBounds(39, 110, 176, 16);
		phys.add(lblT);
		
		T = new JTextField();
		T.setColumns(10);
		T.setBounds(227, 105, 98, 26);
		phys.add(T);
		
		lblVacuumPermittivitysi = new JLabel("Vacuum permittivity [SI units]");
		lblVacuumPermittivitysi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVacuumPermittivitysi.setBounds(29, 11, 186, 16);
		phys.add(lblVacuumPermittivitysi);
		
		eps0 = new JTextField();
		eps0.setColumns(10);
		eps0.setBounds(227, 6, 98, 26);
		phys.add(eps0);
		
		lblVacuumPermeabilitysi = new JLabel("Vacuum permeability [SI units]");
		lblVacuumPermeabilitysi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVacuumPermeabilitysi.setBounds(17, 44, 198, 16);
		phys.add(lblVacuumPermeabilitysi);
		
		mu0 = new JTextField();
		mu0.setColumns(10);
		mu0.setBounds(227, 39, 98, 26);
		phys.add(mu0);
		
		lblBoltzmannConstantk = new JLabel("Elementary charge [C]");
		lblBoltzmannConstantk.setToolTipText("Charge of an electron or hole");
		lblBoltzmannConstantk.setHorizontalAlignment(SwingConstants.TRAILING);
		lblBoltzmannConstantk.setBounds(39, 77, 176, 16);
		phys.add(lblBoltzmannConstantk);
		
		e_charge = new JTextField();
		e_charge.setColumns(10);
		e_charge.setBounds(227, 72, 98, 26);
		phys.add(e_charge);
		
		JPanel semi = new JPanel();
		tabbedPane.addTab("Semiconductor", null, semi, null);
		semi.setLayout(null);
		
		mu_electron = new JTextField();
		mu_electron.setColumns(10);
		mu_electron.setBounds(216, 6, 98, 26);
		semi.add(mu_electron);
		
		JLabel lblNewLabel_3 = new JLabel("Electron mobility [m^2/(V s)]");
		lblNewLabel_3.setToolTipText("Electron mobility");
		lblNewLabel_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3.setBounds(12, 11, 192, 16);
		semi.add(lblNewLabel_3);
		
		JLabel lblNewLabel_4 = new JLabel("Hole mobility [m^2/(V s)]");
		lblNewLabel_4.setToolTipText("Hole mobility");
		lblNewLabel_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4.setBounds(36, 44, 168, 16);
		semi.add(lblNewLabel_4);
		
		mu_hole = new JTextField();
		mu_hole.setColumns(10);
		mu_hole.setBounds(216, 39, 98, 26);
		semi.add(mu_hole);
		
		JLabel lblNewLabel_5 = new JLabel("Intrinsic carrier conc. [1/m^3]");
		lblNewLabel_5.setToolTipText("Semiconductor equilibrium carrier concentration");
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5.setBounds(6, 77, 198, 16);
		semi.add(lblNewLabel_5);
		
		ni_semi = new JTextField();
		ni_semi.setColumns(10);
		ni_semi.setBounds(216, 72, 98, 26);
		semi.add(ni_semi);
		
		JLabel lblNewLabel_6 = new JLabel("Workfunction [eV]");
		lblNewLabel_6.setToolTipText("Semiconductor work function");
		lblNewLabel_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6.setBounds(36, 110, 168, 16);
		semi.add(lblNewLabel_6);
		
		W_semi = new JTextField();
		W_semi.setColumns(10);
		W_semi.setBounds(216, 105, 98, 26);
		semi.add(W_semi);
		
		JLabel lblEbsemi = new JLabel("Bandgap [eV]");
		lblEbsemi.setToolTipText("Semiconductor band gap");
		lblEbsemi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblEbsemi.setBounds(36, 143, 168, 16);
		semi.add(lblEbsemi);
		
		E_b_semi = new JTextField();
		E_b_semi.setColumns(10);
		E_b_semi.setBounds(216, 138, 98, 26);
		semi.add(E_b_semi);
		
		JLabel lblNewLabel_5_1 = new JLabel("Radiative recomb. rate [m^3/s]");
		lblNewLabel_5_1.setToolTipText("");
		lblNewLabel_5_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_1.setBounds(338, 11, 211, 16);
		semi.add(lblNewLabel_5_1);
		
		k_rad_semi = new JTextField();
		k_rad_semi.setColumns(10);
		k_rad_semi.setBounds(561, 6, 98, 26);
		semi.add(k_rad_semi);

		
		other = new JPanel();
		tabbedPane.addTab("Material parameters", null, other, null);
		other.setLayout(null);
		
		lblNewLabel_5_2 = new JLabel("n-type default doping conc. [1/m^3]");
		lblNewLabel_5_2.setToolTipText("");
		lblNewLabel_5_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_2.setBounds(326, 11, 236, 16);
		other.add(lblNewLabel_5_2);
		
		n_default_doping = new JTextField();
		n_default_doping.setColumns(10);
		n_default_doping.setBounds(574, 6, 98, 26);
		other.add(n_default_doping);
		
		lblNewLabel_5_3 = new JLabel("p-type default doping conc. [1/m^3]");
		lblNewLabel_5_3.setToolTipText("");
		lblNewLabel_5_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_3.setBounds(326, 44, 236, 16);
		other.add(lblNewLabel_5_3);
		
		p_default_doping = new JTextField();
		p_default_doping.setColumns(10);
		p_default_doping.setBounds(574, 39, 98, 26);
		other.add(p_default_doping);
		
		lblNewLabel_5_4 = new JLabel("n-type light doping conc. [1/m^3]");
		lblNewLabel_5_4.setToolTipText("");
		lblNewLabel_5_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_4.setBounds(326, 77, 236, 16);
		other.add(lblNewLabel_5_4);
		
		n_light_doping = new JTextField();
		n_light_doping.setColumns(10);
		n_light_doping.setBounds(574, 72, 98, 26);
		other.add(n_light_doping);
		
		lblNewLabel_5_5 = new JLabel("p-type light doping conc. [1/m^3]");
		lblNewLabel_5_5.setToolTipText("");
		lblNewLabel_5_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_5.setBounds(326, 110, 236, 16);
		other.add(lblNewLabel_5_5);
		
		p_light_doping = new JTextField();
		p_light_doping.setColumns(10);
		p_light_doping.setBounds(574, 105, 98, 26);
		other.add(p_light_doping);
		
		lblNewLabel_5_6 = new JLabel("n-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_6.setToolTipText("");
		lblNewLabel_5_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_6.setBounds(326, 143, 236, 16);
		other.add(lblNewLabel_5_6);
		
		n_heavy_doping = new JTextField();
		n_heavy_doping.setColumns(10);
		n_heavy_doping.setBounds(574, 138, 98, 26);
		other.add(n_heavy_doping);
		
		lblNewLabel_5_7 = new JLabel("p-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_7.setToolTipText("");
		lblNewLabel_5_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_7.setBounds(326, 176, 236, 16);
		other.add(lblNewLabel_5_7);
		
		p_heavy_doping = new JTextField();
		p_heavy_doping.setColumns(10);
		p_heavy_doping.setBounds(574, 171, 98, 26);
		other.add(p_heavy_doping);
		
		lblNewLabel_5_8 = new JLabel("Electron sat. velocity [m/s]");
		lblNewLabel_5_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_8.setBounds(12, 176, 192, 16);
		semi.add(lblNewLabel_5_8);
		
		v_sat_n = new JTextField();
		v_sat_n.setColumns(10);
		v_sat_n.setBounds(216, 171, 98, 26);
		semi.add(v_sat_n);
		
		lblNewLabel_5_9 = new JLabel("Hole sat. velocity [m/s]");
		lblNewLabel_5_9.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9.setBounds(12, 209, 192, 16);
		semi.add(lblNewLabel_5_9);
		
		v_sat_p = new JTextField();
		v_sat_p.setColumns(10);
		v_sat_p.setBounds(216, 204, 98, 26);
		semi.add(v_sat_p);
		
		lblNewLabel_5_10 = new JLabel("SRH recomb. rate n [1/s]");
		lblNewLabel_5_10.setToolTipText("");
		lblNewLabel_5_10.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_10.setBounds(357, 44, 192, 16);
		semi.add(lblNewLabel_5_10);
		
		k_SRH_n_semi = new JTextField();
		k_SRH_n_semi.setColumns(10);
		k_SRH_n_semi.setBounds(561, 39, 98, 26);
		semi.add(k_SRH_n_semi);
		
		lblNewLabel_5_11 = new JLabel("SRH recomb. rate p [1/s]");
		lblNewLabel_5_11.setToolTipText("");
		lblNewLabel_5_11.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_11.setBounds(367, 77, 182, 16);
		semi.add(lblNewLabel_5_11);
		
		k_SRH_p_semi = new JTextField();
		k_SRH_p_semi.setColumns(10);
		k_SRH_p_semi.setBounds(561, 72, 98, 26);
		semi.add(k_SRH_p_semi);
		
		lblNewLabel_5_12 = new JLabel("Auger recomb. rate n [m^6/s]");
		lblNewLabel_5_12.setToolTipText("");
		lblNewLabel_5_12.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12.setBounds(357, 110, 192, 16);
		semi.add(lblNewLabel_5_12);
		
		k_aug_n_semi = new JTextField();
		k_aug_n_semi.setColumns(10);
		k_aug_n_semi.setBounds(561, 105, 98, 26);
		semi.add(k_aug_n_semi);
		
		lblNewLabel_5_13 = new JLabel("Auger recomb. rate p [m^6/s]");
		lblNewLabel_5_13.setToolTipText("");
		lblNewLabel_5_13.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_13.setBounds(357, 143, 192, 16);
		semi.add(lblNewLabel_5_13);
		
		k_aug_p_semi = new JTextField();
		k_aug_p_semi.setColumns(10);
		k_aug_p_semi.setBounds(561, 138, 98, 26);
		semi.add(k_aug_p_semi);
		
		JLabel lblNewLabel_5_9_1 = new JLabel("Dielectric constant");
		lblNewLabel_5_9_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9_1.setBounds(12, 243, 192, 16);
		semi.add(lblNewLabel_5_9_1);
		
		eps_r_semi = new JTextField();
		eps_r_semi.setColumns(10);
		eps_r_semi.setBounds(216, 238, 98, 26);
		semi.add(eps_r_semi);
		
		JLabel lblNewLabel_5_12_1 = new JLabel("Doping dep. mobility factor n");
		lblNewLabel_5_12_1.setToolTipText("");
		lblNewLabel_5_12_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12_1.setBounds(357, 176, 192, 16);
		semi.add(lblNewLabel_5_12_1);
		
		a_factor_n = new JTextField();
		a_factor_n.setColumns(10);
		a_factor_n.setBounds(561, 171, 98, 26);
		semi.add(a_factor_n);
		
		JLabel lblNewLabel_5_12_2 = new JLabel("Doping dep. mobility factor p");
		lblNewLabel_5_12_2.setToolTipText("");
		lblNewLabel_5_12_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12_2.setBounds(357, 209, 192, 16);
		semi.add(lblNewLabel_5_12_2);
		
		a_factor_p = new JTextField();
		a_factor_p.setColumns(10);
		a_factor_p.setBounds(561, 204, 98, 26);
		semi.add(a_factor_p);
		
		lblDielectricRelPermittivity = new JLabel("Dielectric rel. permittivity");
		lblDielectricRelPermittivity.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDielectricRelPermittivity.setBounds(6, 11, 208, 16);
		other.add(lblDielectricRelPermittivity);
		
		dielectric_eps_r = new JTextField();
		dielectric_eps_r.setColumns(10);
		dielectric_eps_r.setBounds(226, 6, 98, 26);
		other.add(dielectric_eps_r);
		
		lblFerromagnetRelPermeability = new JLabel("Ferromagnet rel. permeability");
		lblFerromagnetRelPermeability.setHorizontalAlignment(SwingConstants.TRAILING);
		lblFerromagnetRelPermeability.setBounds(6, 44, 208, 16);
		other.add(lblFerromagnetRelPermeability);
		
		ferromagnet_mu_r = new JTextField();
		ferromagnet_mu_r.setColumns(10);
		ferromagnet_mu_r.setBounds(226, 39, 98, 26);
		other.add(ferromagnet_mu_r);
		
		lblStaticChargeDensity = new JLabel("Static charge density [C/m^3]");
		lblStaticChargeDensity.setHorizontalAlignment(SwingConstants.TRAILING);
		lblStaticChargeDensity.setBounds(6, 77, 208, 16);
		other.add(lblStaticChargeDensity);
		
		staticcharge_density = new JTextField();
		staticcharge_density.setColumns(10);
		staticcharge_density.setBounds(226, 72, 98, 26);
		other.add(staticcharge_density);
		
		lblCurrentSourceMobilit = new JLabel("Current source rel. mobility");
		lblCurrentSourceMobilit.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCurrentSourceMobilit.setBounds(6, 110, 208, 16);
		other.add(lblCurrentSourceMobilit);
		
		currentsource_mobility = new JTextField();
		currentsource_mobility.setColumns(10);
		currentsource_mobility.setBounds(226, 105, 98, 26);
		other.add(currentsource_mobility);
		
		JLabel lblVoltageSourceMax = new JLabel("Voltage source max EMF (V/m)");
		lblVoltageSourceMax.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVoltageSourceMax.setBounds(6, 143, 208, 16);
		other.add(lblVoltageSourceMax);
		
		max_EMF = new JTextField();
		max_EMF.setColumns(10);
		max_EMF.setBounds(226, 138, 98, 26);
		other.add(max_EMF);
		
		JLabel lblCurrentSourceMax = new JLabel("Current source max (A/m^2)");
		lblCurrentSourceMax.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCurrentSourceMax.setBounds(6, 176, 208, 16);
		other.add(lblCurrentSourceMax);
		
		max_current = new JTextField();
		max_current.setColumns(10);
		max_current.setBounds(226, 171, 98, 26);
		other.add(max_current);
		
		lblDefaultAcFrequency = new JLabel("Default AC frequency (Hz)");
		lblDefaultAcFrequency.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDefaultAcFrequency.setBounds(6, 209, 208, 16);
		other.add(lblDefaultAcFrequency);
		
		default_AC_freq = new JTextField();
		default_AC_freq.setColumns(10);
		default_AC_freq.setBounds(226, 204, 98, 26);
		other.add(default_AC_freq);
		
		btn_apply = new JButton("Apply changes");
		btn_apply.setBounds(255, 372, 144, 29);
		panel.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(549, 372, 144, 29);
		panel.add(btn_cancel);
		
		btn_reset = new JButton("Reset to defaults");
		btn_reset.setBounds(403, 372, 144, 29);
		panel.add(btn_reset);
		
		btn_presets = new JButton("Presets");
		btn_presets.setBounds(10, 372, 144, 29);
		panel.add(btn_presets);
		
		//.add(tabbedPane);
	}
	
	public void initialize(Simulation e) {
		this.e = e;
		btn_apply.addActionListener(this);
		btn_reset.addActionListener(this);
		btn_cancel.addActionListener(this);
		btn_presets.addActionListener(this);
		setVisible(false);
	}
	
	public void applyPreset(Preset p) {

		Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
		JsonObject advsettings = new JsonObject();
		writeAdvancedSettings(gson, advsettings);
		String json = gson.toJson(advsettings);
		
		p.applyPreset(e);
		
		storeAdvancedSettings();
		
		try {
			JsonReader fstr = new JsonReader(new StringReader(json));
			fstr.beginObject();
			this.readAdvancedSettings(gson, fstr);
			fstr.endObject();
		} catch (IOException | RuntimeException e1) {
			e1.printStackTrace();
		}
	}

	public void pickPresets() {
		JList<Preset> tmplist = new JList<>(Preset.values());

		tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		tmplist.setVisibleRowCount(5);
		tmplist.setSelectedValue(Preset.DEFAULT, true);

		JScrollPane scrollPane = new JScrollPane(tmplist);

		int result = JOptionPane.showConfirmDialog(
		null,
		scrollPane,
		"Select preset",
		JOptionPane.OK_CANCEL_OPTION,
		JOptionPane.PLAIN_MESSAGE
		);

		if (result == JOptionPane.OK_OPTION) {
			Preset selected = tmplist.getSelectedValue();

			if (selected != null) {
				applyPreset(selected);
			}
		}
	}

	private JTextField ni_metal_high;
	private JTextField ni_metal_low;
	private JTextField max_EMF;
	private JTextField max_current;
	private JTextField eps_r_semi;
	private JTextField a_factor_n;
	private JTextField a_factor_p;
	private JLabel lblDefaultAcFrequency;
	private JTextField default_AC_freq;

	public void storeAdvancedSettings() {
		width				.setText(Utils.formatDouble(e.default_width				));
		resolution			.setText(Integer.toString(e.default_resolution		));
		depth				.setText(Utils.formatDouble(e.depth						));
		junction_size		.setText(Integer.toString(e.junction_size			));
		dopant_smoothing_distance		.setText(Integer.toString(e.dopant_smoothing_distance			));

		eps0					.setText(Utils.formatDouble(e.eps0						));
		mu0					.setText(Utils.formatDouble(e.mu0							));
		e_charge				.setText(Utils.formatDouble(e.e_charge					));
		T					.setText(Utils.formatDouble(e.T							));
		
		mu_electron			.setText(Utils.formatDouble(e.mu_electron_semi			));
		mu_hole				.setText(Utils.formatDouble(e.mu_hole_semi				));
		ni_semi				.setText(Utils.formatDouble(e.ni_semi						));
		W_semi				.setText(Utils.formatDouble(e.W_semi/e.eVtoJ				));
		E_b_semi				.setText(Utils.formatDouble(e.E_b_semi/e.eVtoJ			));
		k_rad_semi		.setText(Utils.formatDouble(e.k_rad_semi				));
		k_SRH_n_semi		.setText(Utils.formatDouble(e.k_SRH_n_semi			));
		k_SRH_p_semi		.setText(Utils.formatDouble(e.k_SRH_p_semi			));
		k_aug_n_semi		.setText(Utils.formatDouble(e.k_aug_n_semi			));
		k_aug_p_semi		.setText(Utils.formatDouble(e.k_aug_p_semi			));
		
		v_sat_n				.setText(Utils.formatDouble(e.v_sat_n_semi						));
		v_sat_p				.setText(Utils.formatDouble(e.v_sat_p_semi						));
		n_default_doping		.setText(Utils.formatDouble(e.n_default_doping_concentration		));
		p_default_doping		.setText(Utils.formatDouble(e.p_default_doping_concentration		));
		n_light_doping		.setText(Utils.formatDouble(e.n_light_doping_concentration		));
		p_light_doping		.setText(Utils.formatDouble(e.p_light_doping_concentration		));
		n_heavy_doping		.setText(Utils.formatDouble(e.n_heavy_doping_concentration		));
		p_heavy_doping		.setText(Utils.formatDouble(e.p_heavy_doping_concentration		));
		
		ni_metal				.setText(Utils.formatDouble(e.ni_metal					));
		ni_metal_high				.setText(Utils.formatDouble(e.ni_metal_high					));
		ni_metal_low				.setText(Utils.formatDouble(e.ni_metal_low					));
		W_metal				.setText(Utils.formatDouble(e.W_metal_default/e.eVtoJ		));
		E_b_metal			.setText(Utils.formatDouble(e.E_b_metal/e.eVtoJ			));
		W_metal_high			.setText(Utils.formatDouble(e.W_metal_high/e.eVtoJ		));
		W_metal_low			.setText(Utils.formatDouble(e.W_metal_low/e.eVtoJ			));
		recomb_rate_metal	.setText(Utils.formatDouble(e.k_rad_metal					));

		dielectric_eps_r			.setText(Utils.formatDouble(e.dielectric_eps_r		));
		ferromagnet_mu_r			.setText(Utils.formatDouble(e.ferromagnet_mu_r		));
		staticcharge_density		.setText(Utils.formatDouble(e.staticcharge_density	));
		currentsource_mobility	.setText(Utils.formatDouble(e.currentsource_mobility	));
		max_EMF	.setText(Utils.formatDouble(e.max_EMF	));
		max_current	.setText(Utils.formatDouble(e.max_current	));
		default_AC_freq	.setText(Utils.formatDouble(e.default_AC_freq	));

		a_factor_n		.setText(Utils.formatDouble(e.a_factor_n));
		a_factor_p		.setText(Utils.formatDouble(e.a_factor_p));
		eps_r_semi		.setText(Utils.formatDouble(e.eps_r_semi));
	}
	
	public boolean loadAdvancedSettings(boolean show_warning) {

		try {
			int resolution_tmp			= Integer.valueOf(resolution			.getText());

			if (resolution_tmp != e.resolution) {
				int result = JOptionPane.showConfirmDialog(this, "Changing the resolution will delete all materials. Proceed?", "Message", JOptionPane.YES_NO_OPTION);
				if (result != JOptionPane.OK_OPTION)
				{
					return false;
				}
				
				double memory_estimate = 400.0*8.0*(double)resolution_tmp*(double)resolution_tmp;
				if (memory_estimate > 1e9) {
					int result2 = JOptionPane.showConfirmDialog(this, "Warning: This resolution will use approximately " + Units.SI.toString(memory_estimate, Quantity.INFORMATION) + " of memory. Proceed?", "Message", JOptionPane.YES_NO_OPTION);
					if (result2 != JOptionPane.OK_OPTION)
					{
						return false;
					}
				}
			}

			e.default_width				= Double.valueOf(width				.getText());	
			e.default_resolution = resolution_tmp;

			e.depth						= Double.valueOf(depth				.getText());
			e.junction_size				= Integer.valueOf(junction_size		.getText());
			e.dopant_smoothing_distance				= Integer.valueOf(dopant_smoothing_distance		.getText());

			e.T							= Double.valueOf(T					.getText());
			e.eps0						= Double.valueOf(eps0				.getText());
			e.mu0						= Double.valueOf(mu0					.getText());
			e.e_charge					= Double.valueOf(e_charge			.getText());

			e.mu_electron_semi				= Double.valueOf(mu_electron			.getText());
			e.mu_hole_semi					= Double.valueOf(mu_hole				.getText());
			e.ni_semi					= Double.valueOf(ni_semi				.getText());
			e.W_semi					= Double.valueOf(W_semi				.getText())*e.eVtoJ;
			e.E_b_semi					= Double.valueOf(E_b_semi			.getText())*e.eVtoJ;
			e.k_rad_semi			= Double.valueOf(k_rad_semi	.getText());
			e.k_SRH_n_semi			= Double.valueOf(k_SRH_n_semi	.getText());
			e.k_SRH_p_semi			= Double.valueOf(k_SRH_p_semi	.getText());
			e.k_aug_n_semi			= Double.valueOf(k_aug_n_semi	.getText());
			e.k_aug_p_semi			= Double.valueOf(k_aug_p_semi	.getText());
			e.v_sat_n_semi						= Double.valueOf(v_sat_n				.getText());
			e.v_sat_p_semi						= Double.valueOf(v_sat_p				.getText());
			e.n_default_doping_concentration	= Double.valueOf(n_default_doping	.getText());
			e.p_default_doping_concentration	= Double.valueOf(p_default_doping	.getText());
			e.n_light_doping_concentration		= Double.valueOf(n_light_doping	.getText());
			e.p_light_doping_concentration		= Double.valueOf(p_light_doping	.getText());
			e.n_heavy_doping_concentration		= Double.valueOf(n_heavy_doping	.getText());
			e.p_heavy_doping_concentration		= Double.valueOf(p_heavy_doping	.getText());

			e.ni_metal					= Double.valueOf(ni_metal			.getText());
			e.ni_metal_high					= Double.valueOf(ni_metal_high			.getText());
			e.ni_metal_low					= Double.valueOf(ni_metal_low			.getText());
			e.W_metal_default			= Double.valueOf(W_metal				.getText())*e.eVtoJ;
			e.E_b_metal					= Double.valueOf(E_b_metal			.getText())*e.eVtoJ;
			e.W_metal_high				= Double.valueOf(W_metal_high		.getText())*e.eVtoJ;
			e.W_metal_low				= Double.valueOf(W_metal_low			.getText())*e.eVtoJ;
			e.k_rad_metal			= Double.valueOf(recomb_rate_metal	.getText());

			e.dielectric_eps_r			= Double.valueOf(dielectric_eps_r	.getText());
			e.ferromagnet_mu_r			= Double.valueOf(ferromagnet_mu_r	.getText());
			e.staticcharge_density		= Double.valueOf(staticcharge_density	.getText());
			e.currentsource_mobility	= Double.valueOf(currentsource_mobility	.getText());
			e.max_EMF	= Double.valueOf(max_EMF	.getText());
			e.max_current	= Double.valueOf(max_current	.getText());
			e.default_AC_freq	= Double.valueOf(default_AC_freq	.getText());

			
			e.a_factor_n = Double.valueOf(a_factor_n	.getText());
			e.a_factor_p = Double.valueOf(a_factor_p	.getText());
			e.eps_r_semi = Double.valueOf(eps_r_semi	.getText());

			e.calculateDependentConstants();
			e.reset(false, false);
			e.lastsimspeed = -1;
			e.advsettings_tweaked = true;
			return true;
			
		} catch (NumberFormatException e) {
			JOptionPane.showMessageDialog(this, "Invalid number.", "Error", JOptionPane.OK_OPTION);
			return false;
		}
	}

	public void writeAdvancedSettings (Gson gson, JsonObject advsettings) {
		// Width and resolution are handled in SaveManager
		
		advsettings.addProperty("depth", e.depth 						);
		advsettings.addProperty("junction_size", e.junction_size		);
		advsettings.addProperty("dopant_smoothing_distance", e.dopant_smoothing_distance		);

		advsettings.addProperty("T", e.T								);
		advsettings.addProperty("eps0", e.eps0							);
		advsettings.addProperty("mu0", e.mu0							);
		advsettings.addProperty("e_charge", e.e_charge					);
		
		advsettings.addProperty("mu_electron", e.mu_electron_semi			);
		advsettings.addProperty("mu_hole", e.mu_hole_semi					);
		advsettings.addProperty("ni_semi", e.ni_semi					);
		advsettings.addProperty("W_semi", e.W_semi						);
		advsettings.addProperty("E_b_semi", e.E_b_semi					);
		advsettings.addProperty("k_rad_semi", e.k_rad_semi	);
		advsettings.addProperty("k_SRH_n_semi", e.k_SRH_n_semi	);
		advsettings.addProperty("k_SRH_p_semi", e.k_SRH_p_semi	);
		advsettings.addProperty("k_aug_n_semi", e.k_aug_n_semi	);
		advsettings.addProperty("k_aug_p_semi", e.k_aug_p_semi	);
		
		advsettings.addProperty("v_sat_n", e.v_sat_n_semi						);
		advsettings.addProperty("v_sat_p", e.v_sat_p_semi						);
		advsettings.addProperty("n_default_doping", e.n_default_doping_concentration	);
		advsettings.addProperty("p_default_doping", e.p_default_doping_concentration	);
		advsettings.addProperty("n_light_doping", e.n_light_doping_concentration	);
		advsettings.addProperty("p_light_doping", e.p_light_doping_concentration	);
		advsettings.addProperty("n_heavy_doping", e.n_heavy_doping_concentration	);
		advsettings.addProperty("p_heavy_doping", e.p_heavy_doping_concentration	);
		
		advsettings.addProperty("ni_metal", e.ni_metal						);
		advsettings.addProperty("ni_metal_high", e.ni_metal_high						);
		advsettings.addProperty("ni_metal_low", e.ni_metal_low						);
		advsettings.addProperty("W_metal_default", e.W_metal_default		);
		advsettings.addProperty("E_b_metal", e.E_b_metal					);
		advsettings.addProperty("W_metal_high", e.W_metal_high				);
		advsettings.addProperty("W_metal_low", e.W_metal_low				);
		advsettings.addProperty("k_rad_metal", e.k_rad_metal	);

		advsettings.addProperty("dielectric_eps_r", e.dielectric_eps_r	);
		advsettings.addProperty("ferromagnet_mu_r", e.ferromagnet_mu_r	);
		advsettings.addProperty("staticcharge_density", e.staticcharge_density	);
		advsettings.addProperty("currentsource_mobility", e.currentsource_mobility	);
		advsettings.addProperty("max_EMF", e.max_EMF	);
		advsettings.addProperty("max_current", e.max_current	);
		advsettings.addProperty("default_AC_freq", e.default_AC_freq	);

		advsettings.addProperty("a_factor_n", e.a_factor_n	);
		advsettings.addProperty("a_factor_p", e.a_factor_p	);
		advsettings.addProperty("eps_r_semi", e.eps_r_semi	);
	}
	

	public void readAdvancedSettings (Gson gson, JsonReader fstr) throws JsonIOException, JsonSyntaxException, IOException, RuntimeException {
		
		e.setDefaultParameters();

		while (fstr.hasNext()) {
			String name = fstr.nextName();
			switch (name){
			case "depth": e.depth 						= fstr.nextDouble(); break;
			case "junction_size": e.junction_size 		= fstr.nextInt(); break;
			case "dopant_smoothing_distance": e.dopant_smoothing_distance 		= fstr.nextInt(); break;

			case "T": e.T								= fstr.nextDouble(); break;
			case "eps0": e.eps0							= fstr.nextDouble(); break;
			case "mu0": e.mu0							= fstr.nextDouble(); break;
			case "e_charge": e.e_charge					= fstr.nextDouble(); break;
			
			case "mu_electron": e.mu_electron_semi			= fstr.nextDouble(); break;
			case "mu_hole": e.mu_hole_semi					= fstr.nextDouble(); break;
			case "ni_semi": e.ni_semi					= fstr.nextDouble(); break;
			case "W_semi": e.W_semi						= fstr.nextDouble(); break;
			case "E_b_semi": e.E_b_semi					= fstr.nextDouble(); break;
			case "k_rad_semi": e.k_rad_semi		= fstr.nextDouble(); break;
			case "k_SRH_n_semi": e.k_SRH_n_semi		= fstr.nextDouble(); break;
			case "k_SRH_p_semi": e.k_SRH_p_semi		= fstr.nextDouble(); break;
			case "k_aug_n_semi": e.k_aug_n_semi		= fstr.nextDouble(); break;
			case "k_aug_p_semi": e.k_aug_p_semi		= fstr.nextDouble(); break;
			case "v_sat_n": e.v_sat_n_semi					= fstr.nextDouble(); break;
			case "v_sat_p": e.v_sat_p_semi						= fstr.nextDouble(); break;
			case "n_default_doping": e.n_default_doping_concentration	= fstr.nextDouble(); break;
			case "p_default_doping": e.p_default_doping_concentration	= fstr.nextDouble(); break;
			case "n_light_doping": e.n_light_doping_concentration		= fstr.nextDouble(); break;
			case "p_light_doping": e.p_light_doping_concentration		= fstr.nextDouble(); break;
			case "n_heavy_doping": e.n_heavy_doping_concentration		= fstr.nextDouble(); break;
			case "p_heavy_doping": e.p_heavy_doping_concentration		= fstr.nextDouble(); break;
			
			case "ni_metal": e.ni_metal						= fstr.nextDouble(); break;
			case "ni_metal_high": e.ni_metal_high						= fstr.nextDouble(); break;
			case "ni_metal_low": e.ni_metal_low						= fstr.nextDouble(); break;
			case "W_metal_default": e.W_metal_default		= fstr.nextDouble(); break;
			case "E_b_metal": e.E_b_metal					= fstr.nextDouble(); break;
			case "W_metal_high": e.W_metal_high				= fstr.nextDouble(); break;
			case "W_metal_low": e.W_metal_low				= fstr.nextDouble(); break;
			case "k_rad_metal": e.k_rad_metal			= fstr.nextDouble(); break;

			case "dielectric_eps_r": e.dielectric_eps_r		= fstr.nextDouble(); break;
			case "ferromagnet_mu_r": e.ferromagnet_mu_r		= fstr.nextDouble(); break;
			case "staticcharge_density": e.staticcharge_density	= fstr.nextDouble(); break;
			case "currentsource_mobility": e.currentsource_mobility	= fstr.nextDouble(); break;
			case "max_EMF": e.max_EMF	= fstr.nextDouble(); break;
			case "max_current": e.max_current	= fstr.nextDouble(); break;
			case "default_AC_freq": e.default_AC_freq	= fstr.nextDouble(); break;

			case "a_factor_n": e.a_factor_n	= fstr.nextDouble(); break;
			case "a_factor_p": e.a_factor_p	= fstr.nextDouble(); break;
			case "eps_r_semi": e.eps_r_semi	= fstr.nextDouble(); break;

			default: fstr.skipValue(); break; // skip others
			}
		}
	}

	@Override
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == btn_apply) {
			loadAdvancedSettings(true);
			//boolean success = loadAdvancedSettings(true);
			//if (success)
			//	setVisible(false);
		} else if (ev.getSource() == btn_cancel) {
			setVisible(false);
		} else if (ev.getSource() == btn_reset) {
			applyPreset(Preset.DEFAULT);
		} else if (ev.getSource() == btn_presets) {
			pickPresets();
		}
	}
}
