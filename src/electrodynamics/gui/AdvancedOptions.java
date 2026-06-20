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
import com.google.gson.reflect.TypeToken;
import com.google.gson.stream.JsonReader;

import electrodynamics.MaterialType;
import electrodynamics.Preset;
import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import electrodynamics.units.Units;
import electrodynamics.util.Utils;

import javax.swing.JTabbedPane;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.io.StringReader;
import java.util.HashMap;

public class AdvancedOptions extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	HashMap<MaterialType, String> modified_names_tmp;
	private JPanel sim;
	public JTextField ds;
	public JTextField resolution_x;
	public JTextField depth;
	public JTextField mu_electron_semi;
	public JTextField mu_hole_semi;
	public JTextField chi_semi;
	public JTextField Eg_semi;
	public JTextField g_metal;
	public JTextField W_metal;
	public JTextField Eg_metal;
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
	private JTextField v_sat_n_semi;
	private JLabel lblNewLabel_5_9;
	private JTextField v_sat_p_semi;
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

	public AdvancedOptions(Simulation e) {
		this.e = e;
		setResizable(false);
		setTitle("Advanced settings");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 699, 439);

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
		
		JLabel lblNewLabel = new JLabel("Grid spacing [m]");
		lblNewLabel.setToolTipText("Width of simulation domain");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(45, 11, 168, 16);
		sim.add(lblNewLabel);
		
		ds = new JTextField();
		ds.setBounds(225, 7, 98, 23);
		sim.add(ds);
		ds.setColumns(10);
		
		JLabel lblNewLabel_1 = new JLabel("Grid size x");
		lblNewLabel_1.setToolTipText("Number of grid points in x or y direction. Must be power of 2");
		lblNewLabel_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1.setBounds(27, 44, 186, 16);
		sim.add(lblNewLabel_1);
		
		resolution_x = new JTextField();
		resolution_x.setColumns(10);
		resolution_x.setBounds(225, 40, 98, 23);
		sim.add(resolution_x);
		
		JLabel lblNewLabel_2 = new JLabel("Depth [m]");
		lblNewLabel_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2.setBounds(45, 110, 168, 16);
		sim.add(lblNewLabel_2);
		
		depth = new JTextField();
		depth.setColumns(10);
		depth.setBounds(225, 106, 98, 23);
		sim.add(depth);
		
		lblNewLabel_7 = new JLabel("Junction smoothing [px]");
		lblNewLabel_7.setToolTipText("Smooths junction between materials with different chemical potential (ie. metal-semiconductor junctions)");
		lblNewLabel_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_7.setBounds(45, 143, 168, 16);
		sim.add(lblNewLabel_7);
		
		junction_size = new JTextField();
		junction_size.setColumns(10);
		junction_size.setBounds(225, 139, 98, 23);
		sim.add(junction_size);
		
		lblNewLabel_8 = new JLabel("Dopant smoothing [px]");
		lblNewLabel_8.setToolTipText("Smooths dopant density, use when modelling non abrupt PN junctions.");
		lblNewLabel_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_8.setBounds(45, 176, 168, 16);
		sim.add(lblNewLabel_8);
		
		dopant_smoothing_distance = new JTextField();
		dopant_smoothing_distance.setColumns(10);
		dopant_smoothing_distance.setBounds(225, 172, 98, 23);
		sim.add(dopant_smoothing_distance);
		
		JLabel lblNimetal = new JLabel("Effective DOS [1/m^3]");
		lblNimetal.setToolTipText("Metal equilibrium carrier concentration");
		lblNimetal.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNimetal.setBounds(6, 11, 203, 16);
		metal.add(lblNimetal);
		
		g_metal = new JTextField();
		g_metal.setColumns(10);
		g_metal.setBounds(221, 7, 98, 23);
		metal.add(g_metal);
		
		JLabel lblNewLabel_1_1 = new JLabel("Workfunction (default) [eV]");
		lblNewLabel_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1_1.setBounds(33, 77, 176, 16);
		metal.add(lblNewLabel_1_1);
		
		W_metal = new JTextField();
		W_metal.setColumns(10);
		W_metal.setBounds(221, 73, 98, 23);
		metal.add(W_metal);
		
		JLabel lblNewLabel_2_1 = new JLabel("Metal \"bandgap\" [eV]");
		lblNewLabel_2_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2_1.setBounds(33, 44, 176, 16);
		metal.add(lblNewLabel_2_1);
		
		Eg_metal = new JTextField();
		Eg_metal.setColumns(10);
		Eg_metal.setBounds(221, 40, 98, 23);
		metal.add(Eg_metal);
		
		JLabel lblNewLabel_3_1 = new JLabel("Workfunction (High WF metal) [eV]");
		lblNewLabel_3_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3_1.setBounds(330, 11, 232, 16);
		metal.add(lblNewLabel_3_1);
		
		W_metal_high = new JTextField();
		W_metal_high.setColumns(10);
		W_metal_high.setBounds(574, 7, 98, 23);
		metal.add(W_metal_high);
		
		JLabel lblNewLabel_4_1 = new JLabel("Workfunction (Low WF metal) [eV]");
		lblNewLabel_4_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4_1.setBounds(340, 44, 222, 16);
		metal.add(lblNewLabel_4_1);
		
		W_metal_low = new JTextField();
		W_metal_low.setColumns(10);
		W_metal_low.setBounds(574, 40, 98, 23);
		metal.add(W_metal_low);
		
		JLabel lblNewLabel_6_1 = new JLabel("Recomb. rate [m^3/s]");
		lblNewLabel_6_1.setToolTipText("Radiative recombination rate constant in metal");
		lblNewLabel_6_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6_1.setBounds(39, 110, 170, 16);
		metal.add(lblNewLabel_6_1);
		
		recomb_rate_metal = new JTextField();
		recomb_rate_metal.setColumns(10);
		recomb_rate_metal.setBounds(221, 106, 98, 23);
		metal.add(recomb_rate_metal);
		
		JLabel lblCarrierConchigh = new JLabel("Effective DOS (High cond.) [1/m^3]");
		lblCarrierConchigh.setToolTipText("Metal equilibrium carrier concentration");
		lblCarrierConchigh.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCarrierConchigh.setBounds(331, 77, 231, 16);
		metal.add(lblCarrierConchigh);
		
		g_metal_high = new JTextField();
		g_metal_high.setColumns(10);
		g_metal_high.setBounds(574, 73, 98, 23);
		metal.add(g_metal_high);
		
		JLabel lblCarrierConclow = new JLabel("Effective DOS (Low cond.) [1/m^3]");
		lblCarrierConclow.setToolTipText("Metal equilibrium carrier concentration");
		lblCarrierConclow.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCarrierConclow.setBounds(330, 110, 231, 16);
		metal.add(lblCarrierConclow);
		
		g_metal_low = new JTextField();
		g_metal_low.setColumns(10);
		g_metal_low.setBounds(573, 106, 98, 23);
		metal.add(g_metal_low);
		
		
		JLabel lblT = new JLabel("Temperature [K]");
		lblT.setToolTipText("Global temperature");
		lblT.setHorizontalAlignment(SwingConstants.TRAILING);
		lblT.setBounds(39, 110, 176, 16);
		phys.add(lblT);
		
		T = new JTextField();
		T.setColumns(10);
		T.setBounds(226, 107, 98, 23);
		phys.add(T);
		
		lblVacuumPermittivitysi = new JLabel("Vacuum permittivity [SI units]");
		lblVacuumPermittivitysi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVacuumPermittivitysi.setBounds(29, 11, 186, 16);
		phys.add(lblVacuumPermittivitysi);
		
		eps0 = new JTextField();
		eps0.setColumns(10);
		eps0.setBounds(226, 8, 98, 23);
		phys.add(eps0);
		
		lblVacuumPermeabilitysi = new JLabel("Vacuum permeability [SI units]");
		lblVacuumPermeabilitysi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVacuumPermeabilitysi.setBounds(17, 44, 198, 16);
		phys.add(lblVacuumPermeabilitysi);
		
		mu0 = new JTextField();
		mu0.setColumns(10);
		mu0.setBounds(226, 41, 98, 23);
		phys.add(mu0);
		
		lblBoltzmannConstantk = new JLabel("Elementary charge [C]");
		lblBoltzmannConstantk.setToolTipText("Charge of an electron or hole");
		lblBoltzmannConstantk.setHorizontalAlignment(SwingConstants.TRAILING);
		lblBoltzmannConstantk.setBounds(39, 77, 176, 16);
		phys.add(lblBoltzmannConstantk);
		
		e_charge = new JTextField();
		e_charge.setColumns(10);
		e_charge.setBounds(226, 74, 98, 23);
		phys.add(e_charge);
		
		JPanel semi = new JPanel();
		tabbedPane.addTab("Semiconductor", null, semi, null);
		semi.setLayout(null);
		
		mu_electron_semi = new JTextField();
		mu_electron_semi.setColumns(10);
		mu_electron_semi.setBounds(216, 8, 98, 23);
		semi.add(mu_electron_semi);
		
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
		
		mu_hole_semi = new JTextField();
		mu_hole_semi.setColumns(10);
		mu_hole_semi.setBounds(216, 41, 98, 23);
		semi.add(mu_hole_semi);
		
		JLabel lblNewLabel_6 = new JLabel("Electron affinity [eV]");
		lblNewLabel_6.setToolTipText("Semiconductor work function");
		lblNewLabel_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6.setBounds(36, 77, 168, 16);
		semi.add(lblNewLabel_6);
		
		chi_semi = new JTextField();
		chi_semi.setColumns(10);
		chi_semi.setBounds(216, 74, 98, 23);
		semi.add(chi_semi);
		
		JLabel lblEbsemi = new JLabel("Bandgap [eV]");
		lblEbsemi.setToolTipText("Semiconductor band gap");
		lblEbsemi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblEbsemi.setBounds(36, 110, 168, 16);
		semi.add(lblEbsemi);
		
		Eg_semi = new JTextField();
		Eg_semi.setColumns(10);
		Eg_semi.setBounds(216, 107, 98, 23);
		semi.add(Eg_semi);
		
		JLabel lblNewLabel_5_1 = new JLabel("Radiative recomb. rate [m^3/s]");
		lblNewLabel_5_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_1.setBounds(338, 11, 211, 16);
		semi.add(lblNewLabel_5_1);
		
		k_rad_semi = new JTextField();
		k_rad_semi.setColumns(10);
		k_rad_semi.setBounds(561, 8, 98, 23);
		semi.add(k_rad_semi);

		
		other = new JPanel();
		tabbedPane.addTab("Material parameters", null, other, null);
		other.setLayout(null);
		
		lblNewLabel_5_2 = new JLabel("n-type default doping conc. [1/m^3]");
		lblNewLabel_5_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_2.setBounds(326, 11, 236, 16);
		other.add(lblNewLabel_5_2);
		
		n_default_doping = new JTextField();
		n_default_doping.setColumns(10);
		n_default_doping.setBounds(574, 8, 98, 23);
		other.add(n_default_doping);
		
		lblNewLabel_5_3 = new JLabel("p-type default doping conc. [1/m^3]");
		lblNewLabel_5_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_3.setBounds(326, 44, 236, 16);
		other.add(lblNewLabel_5_3);
		
		p_default_doping = new JTextField();
		p_default_doping.setColumns(10);
		p_default_doping.setBounds(574, 41, 98, 23);
		other.add(p_default_doping);
		
		lblNewLabel_5_4 = new JLabel("n-type light doping conc. [1/m^3]");
		lblNewLabel_5_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_4.setBounds(326, 77, 236, 16);
		other.add(lblNewLabel_5_4);
		
		n_light_doping = new JTextField();
		n_light_doping.setColumns(10);
		n_light_doping.setBounds(574, 74, 98, 23);
		other.add(n_light_doping);
		
		lblNewLabel_5_5 = new JLabel("p-type light doping conc. [1/m^3]");
		lblNewLabel_5_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_5.setBounds(326, 110, 236, 16);
		other.add(lblNewLabel_5_5);
		
		p_light_doping = new JTextField();
		p_light_doping.setColumns(10);
		p_light_doping.setBounds(574, 107, 98, 23);
		other.add(p_light_doping);
		
		lblNewLabel_5_6 = new JLabel("n-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_6.setBounds(326, 143, 236, 16);
		other.add(lblNewLabel_5_6);
		
		n_heavy_doping = new JTextField();
		n_heavy_doping.setColumns(10);
		n_heavy_doping.setBounds(574, 140, 98, 23);
		other.add(n_heavy_doping);
		
		lblNewLabel_5_7 = new JLabel("p-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_7.setBounds(326, 176, 236, 16);
		other.add(lblNewLabel_5_7);
		
		p_heavy_doping = new JTextField();
		p_heavy_doping.setColumns(10);
		p_heavy_doping.setBounds(574, 173, 98, 23);
		other.add(p_heavy_doping);
		
		lblNewLabel_5_8 = new JLabel("Electron sat. velocity [m/s]");
		lblNewLabel_5_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_8.setBounds(12, 209, 192, 16);
		semi.add(lblNewLabel_5_8);
		
		v_sat_n_semi = new JTextField();
		v_sat_n_semi.setColumns(10);
		v_sat_n_semi.setBounds(216, 206, 98, 23);
		semi.add(v_sat_n_semi);
		
		lblNewLabel_5_9 = new JLabel("Hole sat. velocity [m/s]");
		lblNewLabel_5_9.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9.setBounds(12, 242, 192, 16);
		semi.add(lblNewLabel_5_9);
		
		v_sat_p_semi = new JTextField();
		v_sat_p_semi.setColumns(10);
		v_sat_p_semi.setBounds(216, 239, 98, 23);
		semi.add(v_sat_p_semi);
		
		lblNewLabel_5_10 = new JLabel("SRH recomb. rate n [1/s]");
		lblNewLabel_5_10.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_10.setBounds(357, 44, 192, 16);
		semi.add(lblNewLabel_5_10);
		
		k_SRH_n_semi = new JTextField();
		k_SRH_n_semi.setColumns(10);
		k_SRH_n_semi.setBounds(561, 41, 98, 23);
		semi.add(k_SRH_n_semi);
		
		lblNewLabel_5_11 = new JLabel("SRH recomb. rate p [1/s]");
		lblNewLabel_5_11.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_11.setBounds(367, 77, 182, 16);
		semi.add(lblNewLabel_5_11);
		
		k_SRH_p_semi = new JTextField();
		k_SRH_p_semi.setColumns(10);
		k_SRH_p_semi.setBounds(561, 74, 98, 23);
		semi.add(k_SRH_p_semi);
		
		lblNewLabel_5_12 = new JLabel("Auger recomb. rate n [m^6/s]");
		lblNewLabel_5_12.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12.setBounds(357, 110, 192, 16);
		semi.add(lblNewLabel_5_12);
		
		k_aug_n_semi = new JTextField();
		k_aug_n_semi.setColumns(10);
		k_aug_n_semi.setBounds(561, 107, 98, 23);
		semi.add(k_aug_n_semi);
		
		lblNewLabel_5_13 = new JLabel("Auger recomb. rate p [m^6/s]");
		lblNewLabel_5_13.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_13.setBounds(357, 143, 192, 16);
		semi.add(lblNewLabel_5_13);
		
		k_aug_p_semi = new JTextField();
		k_aug_p_semi.setColumns(10);
		k_aug_p_semi.setBounds(561, 140, 98, 23);
		semi.add(k_aug_p_semi);
		
		JLabel lblNewLabel_5_9_1 = new JLabel("Dielectric constant");
		lblNewLabel_5_9_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9_1.setBounds(357, 242, 192, 16);
		semi.add(lblNewLabel_5_9_1);
		
		eps_r_semi = new JTextField();
		eps_r_semi.setColumns(10);
		eps_r_semi.setBounds(561, 239, 98, 23);
		semi.add(eps_r_semi);
		
		JLabel lblNewLabel_5_12_1 = new JLabel("Doping dep. mob. factor n [1/m^3]");
		lblNewLabel_5_12_1.setToolTipText("Density at which n mobility drops to 1/2 of original value");
		lblNewLabel_5_12_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12_1.setBounds(317, 176, 232, 16);
		semi.add(lblNewLabel_5_12_1);
		
		d_crit_n = new JTextField();
		d_crit_n.setColumns(10);
		d_crit_n.setBounds(561, 173, 98, 23);
		semi.add(d_crit_n);
		
		JLabel lblNewLabel_5_12_2 = new JLabel("Doping dep. mob. factor p [1/m^3]");
		lblNewLabel_5_12_2.setToolTipText("Density at which p mobility drops to 1/2 of original value");
		lblNewLabel_5_12_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12_2.setBounds(317, 209, 232, 16);
		semi.add(lblNewLabel_5_12_2);
		
		d_crit_p = new JTextField();
		d_crit_p.setColumns(10);
		d_crit_p.setBounds(561, 206, 98, 23);
		semi.add(d_crit_p);
		
		lblNewLabel_9 = new JLabel("Eff. cond. band DOS [1/m^3]");
		lblNewLabel_9.setToolTipText("Semiconductor work function");
		lblNewLabel_9.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_9.setBounds(12, 143, 192, 16);
		semi.add(lblNewLabel_9);
		
		gc_semi = new JTextField();
		gc_semi.setColumns(10);
		gc_semi.setBounds(216, 140, 98, 23);
		semi.add(gc_semi);
		
		lblNewLabel_10 = new JLabel("Eff. valence band DOS [1/m^3]");
		lblNewLabel_10.setToolTipText("Semiconductor work function");
		lblNewLabel_10.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_10.setBounds(6, 176, 198, 16);
		semi.add(lblNewLabel_10);
		
		gv_semi = new JTextField();
		gv_semi.setColumns(10);
		gv_semi.setBounds(216, 173, 98, 23);
		semi.add(gv_semi);
		
		lblDielectricRelPermittivity = new JLabel("Dielectric rel. permittivity");
		lblDielectricRelPermittivity.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDielectricRelPermittivity.setBounds(6, 11, 208, 16);
		other.add(lblDielectricRelPermittivity);
		
		dielectric_eps_r = new JTextField();
		dielectric_eps_r.setColumns(10);
		dielectric_eps_r.setBounds(226, 8, 98, 23);
		other.add(dielectric_eps_r);
		
		lblFerromagnetRelPermeability = new JLabel("Ferromagnet rel. permeability");
		lblFerromagnetRelPermeability.setHorizontalAlignment(SwingConstants.TRAILING);
		lblFerromagnetRelPermeability.setBounds(6, 44, 208, 16);
		other.add(lblFerromagnetRelPermeability);
		
		ferromagnet_mu_r = new JTextField();
		ferromagnet_mu_r.setColumns(10);
		ferromagnet_mu_r.setBounds(226, 41, 98, 23);
		other.add(ferromagnet_mu_r);
		
		lblStaticChargeDensity = new JLabel("Static charge density [C/m^3]");
		lblStaticChargeDensity.setHorizontalAlignment(SwingConstants.TRAILING);
		lblStaticChargeDensity.setBounds(6, 77, 208, 16);
		other.add(lblStaticChargeDensity);
		
		staticcharge_density = new JTextField();
		staticcharge_density.setColumns(10);
		staticcharge_density.setBounds(226, 74, 98, 23);
		other.add(staticcharge_density);
		
		lblCurrentSourceMobilit = new JLabel("Current source rel. mobility");
		lblCurrentSourceMobilit.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCurrentSourceMobilit.setBounds(6, 110, 208, 16);
		other.add(lblCurrentSourceMobilit);
		
		currentsource_mobility = new JTextField();
		currentsource_mobility.setColumns(10);
		currentsource_mobility.setBounds(226, 107, 98, 23);
		other.add(currentsource_mobility);
		
		JLabel lblVoltageSourceMax = new JLabel("Voltage source max EMF [V/m]");
		lblVoltageSourceMax.setHorizontalAlignment(SwingConstants.TRAILING);
		lblVoltageSourceMax.setBounds(6, 143, 208, 16);
		other.add(lblVoltageSourceMax);
		
		max_EMF = new JTextField();
		max_EMF.setColumns(10);
		max_EMF.setBounds(226, 140, 98, 23);
		other.add(max_EMF);
		
		JLabel lblCurrentSourceMax = new JLabel("Current source max [A/m^2]");
		lblCurrentSourceMax.setHorizontalAlignment(SwingConstants.TRAILING);
		lblCurrentSourceMax.setBounds(6, 176, 208, 16);
		other.add(lblCurrentSourceMax);
		
		max_current = new JTextField();
		max_current.setColumns(10);
		max_current.setBounds(226, 173, 98, 23);
		other.add(max_current);
		
		lblDefaultAcFrequency = new JLabel("Default AC frequency [Hz]");
		lblDefaultAcFrequency.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDefaultAcFrequency.setBounds(6, 209, 208, 16);
		other.add(lblDefaultAcFrequency);
		
		default_AC_freq = new JTextField();
		default_AC_freq.setColumns(10);
		default_AC_freq.setBounds(226, 206, 98, 23);
		other.add(default_AC_freq);
		
		btn_apply = new JButton("Apply changes");
		btn_apply.setBounds(255, 372, 144, 23);
		panel.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(539, 372, 134, 23);
		panel.add(btn_cancel);
		
		btn_reset = new JButton("Reset to defaults");
		btn_reset.setBounds(403, 372, 133, 23);
		panel.add(btn_reset);
		
		btn_presets = new JButton("Presets");
		btn_presets.setBounds(10, 372, 144, 23);
		panel.add(btn_presets);
		
		addEnterKey(phys);
		addEnterKey(semi);
		addEnterKey(metal);
		
		JLabel lblNewLabel_3_2 = new JLabel("Electron mobility [m^2/(V s)]");
		lblNewLabel_3_2.setToolTipText("Electron mobility");
		lblNewLabel_3_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3_2.setBounds(17, 143, 192, 16);
		metal.add(lblNewLabel_3_2);
		
		mu_electron_metal = new JTextField();
		mu_electron_metal.setColumns(10);
		mu_electron_metal.setBounds(221, 139, 98, 23);
		metal.add(mu_electron_metal);
		
		JLabel lblNewLabel_4_2 = new JLabel("Hole mobility [m^2/(V s)]");
		lblNewLabel_4_2.setToolTipText("Hole mobility");
		lblNewLabel_4_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4_2.setBounds(41, 176, 168, 16);
		metal.add(lblNewLabel_4_2);
		
		mu_hole_metal = new JTextField();
		mu_hole_metal.setColumns(10);
		mu_hole_metal.setBounds(221, 172, 98, 23);
		metal.add(mu_hole_metal);
		
		JLabel lblNewLabel_5_8_1 = new JLabel("Electron sat. velocity [m/s]");
		lblNewLabel_5_8_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_8_1.setBounds(17, 208, 192, 16);
		metal.add(lblNewLabel_5_8_1);
		
		v_sat_n_metal = new JTextField();
		v_sat_n_metal.setColumns(10);
		v_sat_n_metal.setBounds(221, 204, 98, 23);
		metal.add(v_sat_n_metal);
		
		JLabel lblNewLabel_5_9_2 = new JLabel("Hole sat. velocity [m/s]");
		lblNewLabel_5_9_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9_2.setBounds(17, 241, 192, 16);
		metal.add(lblNewLabel_5_9_2);
		
		v_sat_p_metal = new JTextField();
		v_sat_p_metal.setColumns(10);
		v_sat_p_metal.setBounds(221, 237, 98, 23);
		metal.add(v_sat_p_metal);
		
		JLabel lblNewLabel_5_9_1_1 = new JLabel("Dielectric constant");
		lblNewLabel_5_9_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9_1_1.setBounds(370, 143, 192, 16);
		metal.add(lblNewLabel_5_9_1_1);
		
		eps_r_metal = new JTextField();
		eps_r_metal.setColumns(10);
		eps_r_metal.setBounds(574, 139, 98, 23);
		metal.add(eps_r_metal);
		addEnterKey(other);
		
		lblOpenSwitchMobility = new JLabel("Switch open mobility");
		lblOpenSwitchMobility.setHorizontalAlignment(SwingConstants.TRAILING);
		lblOpenSwitchMobility.setBounds(6, 242, 208, 16);
		other.add(lblOpenSwitchMobility);
		
		switch_mobility = new JTextField();
		switch_mobility.setColumns(10);
		switch_mobility.setBounds(226, 239, 98, 23);
		other.add(switch_mobility);
		addEnterKey(sim);
		
		lblNewLabel_5 = new JLabel("Grid size y");
		lblNewLabel_5.setToolTipText("Number of grid points in x or y direction. Must be power of 2");
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5.setBounds(27, 77, 186, 16);
		sim.add(lblNewLabel_5);
		
		resolution_y = new JTextField();
		resolution_y.setColumns(10);
		resolution_y.setBounds(225, 73, 98, 23);
		sim.add(resolution_y);
	}
	
	public void initialize() {
		btn_apply.addActionListener(this);
		btn_reset.addActionListener(this);
		btn_cancel.addActionListener(this);
		btn_presets.addActionListener(this);
		setLocationRelativeTo(null);
		setVisible(false);
	}
	
	public void addEnterKey(JPanel panel) {
		for (Component c : panel.getComponents()) {
			if (c instanceof JTextField)
				((JTextField) c).addActionListener(this);
		}
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
			this.readAdvancedSettings(gson, fstr, Preset.DEFAULT);
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
        scrollPane.setPreferredSize(new Dimension(300, 250));

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

	private JTextField g_metal_high;
	private JTextField g_metal_low;
	private JTextField max_EMF;
	private JTextField max_current;
	private JTextField eps_r_semi;
	private JTextField d_crit_n;
	private JTextField d_crit_p;
	private JLabel lblDefaultAcFrequency;
	private JTextField default_AC_freq;
	private JLabel lblNewLabel_9;
	private JTextField gc_semi;
	private JLabel lblNewLabel_10;
	private JTextField gv_semi;
	private JTextField mu_electron_metal;
	private JTextField mu_hole_metal;
	private JTextField v_sat_n_metal;
	private JTextField v_sat_p_metal;
	private JTextField eps_r_metal;
	private JLabel lblOpenSwitchMobility;
	private JTextField switch_mobility;
	private JLabel lblNewLabel_5;
	private JTextField resolution_y;

	public void storeAdvancedSettings() {
		ds				.setText(Utils.formatDouble(e.ds				));
		resolution_x			.setText(Integer.toString(e.default_resolution_x		));
		resolution_y			.setText(Integer.toString(e.default_resolution_y		));
		depth				.setText(Utils.formatDouble(e.depth						));
		junction_size		.setText(Integer.toString(e.junction_size			));
		dopant_smoothing_distance		.setText(Integer.toString(e.dopant_smoothing_distance			));

		eps0				.setText(Utils.formatDouble(e.eps0						));
		mu0					.setText(Utils.formatDouble(e.mu0							));
		e_charge			.setText(Utils.formatDouble(e.e_charge					));
		T					.setText(Utils.formatDouble(e.T							));
		
		mu_electron_semi			.setText(Utils.formatDouble(e.mu_electron_semi			));
		mu_hole_semi				.setText(Utils.formatDouble(e.mu_hole_semi				));
		gc_semi				.setText(Utils.formatDouble(e.gc_semi						));
		gv_semi				.setText(Utils.formatDouble(e.gv_semi						));
		chi_semi				.setText(Utils.formatDouble(e.chi_semi/e.eVtoJ				));
		Eg_semi				.setText(Utils.formatDouble(e.Eg_semi/e.eVtoJ			));
		k_rad_semi		.setText(Utils.formatDouble(e.k_rad_semi				));
		k_SRH_n_semi		.setText(Utils.formatDouble(e.k_SRH_n_semi			));
		k_SRH_p_semi		.setText(Utils.formatDouble(e.k_SRH_p_semi			));
		k_aug_n_semi		.setText(Utils.formatDouble(e.k_aug_n_semi			));
		k_aug_p_semi		.setText(Utils.formatDouble(e.k_aug_p_semi			));
		
		v_sat_n_semi			.setText(Utils.formatDouble(e.v_sat_n_semi						));
		v_sat_p_semi			.setText(Utils.formatDouble(e.v_sat_p_semi						));
		n_default_doping		.setText(Utils.formatDouble(e.n_default_doping_concentration		));
		p_default_doping		.setText(Utils.formatDouble(e.p_default_doping_concentration		));
		n_light_doping		.setText(Utils.formatDouble(e.n_light_doping_concentration		));
		p_light_doping		.setText(Utils.formatDouble(e.p_light_doping_concentration		));
		n_heavy_doping		.setText(Utils.formatDouble(e.n_heavy_doping_concentration		));
		p_heavy_doping		.setText(Utils.formatDouble(e.p_heavy_doping_concentration		));
		
		g_metal				.setText(Utils.formatDouble(e.g_metal					));
		g_metal_high				.setText(Utils.formatDouble(e.g_metal_high					));
		g_metal_low				.setText(Utils.formatDouble(e.g_metal_low					));
		W_metal				.setText(Utils.formatDouble(e.W_metal_default/e.eVtoJ		));
		Eg_metal			.setText(Utils.formatDouble(e.Eg_metal/e.eVtoJ			));
		W_metal_high			.setText(Utils.formatDouble(e.W_metal_high/e.eVtoJ		));
		W_metal_low			.setText(Utils.formatDouble(e.W_metal_low/e.eVtoJ			));
		recomb_rate_metal	.setText(Utils.formatDouble(e.k_rad_metal					));

		mu_electron_metal			.setText(Utils.formatDouble(e.mu_electron_metal			));
		mu_hole_metal				.setText(Utils.formatDouble(e.mu_hole_metal				));
		v_sat_n_metal			.setText(Utils.formatDouble(e.v_sat_n_metal						));
		v_sat_p_metal			.setText(Utils.formatDouble(e.v_sat_p_metal						));
		eps_r_metal		.setText(Utils.formatDouble(e.eps_r_metal));

		dielectric_eps_r			.setText(Utils.formatDouble(e.dielectric_eps_r		));
		ferromagnet_mu_r			.setText(Utils.formatDouble(e.ferromagnet_mu_r		));
		staticcharge_density		.setText(Utils.formatDouble(e.staticcharge_density	));
		currentsource_mobility	.setText(Utils.formatDouble(e.currentsource_mobility	));
		max_EMF	.setText(Utils.formatDouble(e.max_EMF	));
		max_current	.setText(Utils.formatDouble(e.max_current	));
		default_AC_freq	.setText(Utils.formatDouble(e.default_AC_freq	));
		switch_mobility	.setText(Utils.formatDouble(e.switch_open_mobility	));

		d_crit_n		.setText(Utils.formatDouble(e.d_crit_n));
		d_crit_p		.setText(Utils.formatDouble(e.d_crit_p));
		eps_r_semi		.setText(Utils.formatDouble(e.eps_r_semi));

		modified_names_tmp = new HashMap<>(e.modified_names);
	}
	
	public boolean loadAdvancedSettings(boolean show_warning) {

		try {
			int resolution_tmp_x			= Integer.valueOf(resolution_x			.getText());
			int resolution_tmp_y			= Integer.valueOf(resolution_y			.getText());

			if (resolution_tmp_x != e.nx || resolution_tmp_y != e.ny) {
				int result = JOptionPane.showConfirmDialog(this, "Changing the resolution will delete all materials. Proceed?", "Message", JOptionPane.YES_NO_OPTION);
				if (result != JOptionPane.OK_OPTION)
				{
					return false;
				}
				
				if(resolution_tmp_x < 4) {
					JOptionPane.showMessageDialog(this, "X resolution will be set to the minimum of 4.", "Message", JOptionPane.OK_OPTION);
					resolution_tmp_x = 4;
					resolution_x.setText(Integer.toString(resolution_tmp_x));
				}
				
				if(resolution_tmp_y < 4) {
					JOptionPane.showMessageDialog(this, "Y resolution will be set to the minimum of 4.", "Message", JOptionPane.OK_OPTION);
					resolution_tmp_y = 4;
					resolution_y.setText(Integer.toString(resolution_tmp_y));
				}

				int log2_resolution_x = (int) Math.round(Math.log(resolution_tmp_x)/Math.log(2));
				if(1 << log2_resolution_x != resolution_tmp_x) {
					JOptionPane.showMessageDialog(this, "X resolution will be rounded to the nearest power of 2.", "Message", JOptionPane.OK_OPTION);
					resolution_tmp_x = 1 << log2_resolution_x;
					resolution_x.setText(Integer.toString(resolution_tmp_x));
				}
				

				int log2_resolution_y = (int) Math.round(Math.log(resolution_tmp_y)/Math.log(2));
				if(1 << log2_resolution_y != resolution_tmp_y) {
					JOptionPane.showMessageDialog(this, "Y resolution will be rounded to the nearest power of 2.", "Message", JOptionPane.OK_OPTION);
					resolution_tmp_y = 1 << log2_resolution_y;
					resolution_y.setText(Integer.toString(resolution_tmp_y));
				}
				
				double memory_estimate = 400.0*8.0*(double)resolution_tmp_x*(double)resolution_tmp_y;
				if (memory_estimate > 1e9) {
					int result2 = JOptionPane.showConfirmDialog(this, "Warning: This resolution will use approximately " + Units.SI.toString(memory_estimate, Quantity.INFORMATION) + " of memory. Proceed?", "Message", JOptionPane.YES_NO_OPTION);
					if (result2 != JOptionPane.OK_OPTION)
					{
						return false;
					}
				}
			}

			e.ds				= Double.valueOf(ds				.getText());	
			e.default_resolution_x = resolution_tmp_x;
			e.default_resolution_y = resolution_tmp_y;

			e.depth						= Double.valueOf(depth				.getText());
			e.junction_size				= Integer.valueOf(junction_size		.getText());
			e.dopant_smoothing_distance	= Integer.valueOf(dopant_smoothing_distance		.getText());

			e.T							= Double.valueOf(T					.getText());
			e.eps0						= Double.valueOf(eps0				.getText());
			e.mu0						= Double.valueOf(mu0					.getText());
			e.e_charge					= Double.valueOf(e_charge			.getText());

			e.mu_electron_semi				= Double.valueOf(mu_electron_semi			.getText());
			e.mu_hole_semi					= Double.valueOf(mu_hole_semi				.getText());
			e.gc_semi					= Double.valueOf(gc_semi				.getText());
			e.gv_semi					= Double.valueOf(gv_semi				.getText());
			e.chi_semi					= Double.valueOf(chi_semi				.getText())*e.eVtoJ;
			e.Eg_semi					= Double.valueOf(Eg_semi			.getText())*e.eVtoJ;
			e.k_rad_semi			= Double.valueOf(k_rad_semi	.getText());
			e.k_SRH_n_semi			= Double.valueOf(k_SRH_n_semi	.getText());
			e.k_SRH_p_semi			= Double.valueOf(k_SRH_p_semi	.getText());
			e.k_aug_n_semi			= Double.valueOf(k_aug_n_semi	.getText());
			e.k_aug_p_semi			= Double.valueOf(k_aug_p_semi	.getText());
			e.v_sat_n_semi						= Double.valueOf(v_sat_n_semi				.getText());
			e.v_sat_p_semi						= Double.valueOf(v_sat_p_semi				.getText());
			e.n_default_doping_concentration	= Double.valueOf(n_default_doping	.getText());
			e.p_default_doping_concentration	= Double.valueOf(p_default_doping	.getText());
			e.n_light_doping_concentration		= Double.valueOf(n_light_doping	.getText());
			e.p_light_doping_concentration		= Double.valueOf(p_light_doping	.getText());
			e.n_heavy_doping_concentration		= Double.valueOf(n_heavy_doping	.getText());
			e.p_heavy_doping_concentration		= Double.valueOf(p_heavy_doping	.getText());

			e.g_metal					= Double.valueOf(g_metal			.getText());
			e.g_metal_high					= Double.valueOf(g_metal_high			.getText());
			e.g_metal_low					= Double.valueOf(g_metal_low			.getText());
			e.W_metal_default			= Double.valueOf(W_metal				.getText())*e.eVtoJ;
			e.Eg_metal					= Double.valueOf(Eg_metal			.getText())*e.eVtoJ;
			e.W_metal_high				= Double.valueOf(W_metal_high		.getText())*e.eVtoJ;
			e.W_metal_low				= Double.valueOf(W_metal_low			.getText())*e.eVtoJ;
			e.k_rad_metal			= Double.valueOf(recomb_rate_metal	.getText());

			e.mu_electron_metal				= Double.valueOf(mu_electron_metal			.getText());
			e.mu_hole_metal					= Double.valueOf(mu_hole_metal				.getText());
			e.v_sat_n_metal						= Double.valueOf(v_sat_n_metal				.getText());
			e.v_sat_p_metal						= Double.valueOf(v_sat_p_metal				.getText());
			e.eps_r_metal = Double.valueOf(eps_r_metal	.getText());

			e.dielectric_eps_r			= Double.valueOf(dielectric_eps_r	.getText());
			e.ferromagnet_mu_r			= Double.valueOf(ferromagnet_mu_r	.getText());
			e.staticcharge_density		= Double.valueOf(staticcharge_density	.getText());
			e.currentsource_mobility	= Double.valueOf(currentsource_mobility	.getText());
			e.max_EMF	= Double.valueOf(max_EMF	.getText());
			e.max_current	= Double.valueOf(max_current	.getText());
			e.default_AC_freq	= Double.valueOf(default_AC_freq	.getText());
			e.switch_open_mobility	= Double.valueOf(switch_mobility	.getText());

			
			e.d_crit_n = Double.valueOf(d_crit_n	.getText());
			e.d_crit_p = Double.valueOf(d_crit_p	.getText());
			e.eps_r_semi = Double.valueOf(eps_r_semi	.getText());
			
			e.modified_names = new HashMap<>(modified_names_tmp);

			e.calculateDependentConstants();
			e.reset(false, null);
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
		advsettings.addProperty("gc_semi", e.gc_semi					);
		advsettings.addProperty("gv_semi", e.gv_semi					);
		advsettings.addProperty("chi_semi", e.chi_semi						);
		advsettings.addProperty("Eg_semi", e.Eg_semi					);
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
		
		advsettings.addProperty("g_metal", e.g_metal						);
		advsettings.addProperty("g_metal_high", e.g_metal_high						);
		advsettings.addProperty("g_metal_low", e.g_metal_low						);
		advsettings.addProperty("W_metal_default", e.W_metal_default		);
		advsettings.addProperty("Eg_metal", e.Eg_metal					);
		advsettings.addProperty("W_metal_high", e.W_metal_high				);
		advsettings.addProperty("W_metal_low", e.W_metal_low				);
		advsettings.addProperty("k_rad_metal", e.k_rad_metal	);

		advsettings.addProperty("mu_electron_metal", e.mu_electron_metal			);
		advsettings.addProperty("mu_hole_metal", e.mu_hole_metal					);
		advsettings.addProperty("v_sat_n_metal", e.v_sat_n_metal						);
		advsettings.addProperty("v_sat_p_metal", e.v_sat_p_metal						);
		advsettings.addProperty("eps_r_metal", e.eps_r_metal	);

		advsettings.addProperty("dielectric_eps_r", e.dielectric_eps_r	);
		advsettings.addProperty("ferromagnet_mu_r", e.ferromagnet_mu_r	);
		advsettings.addProperty("staticcharge_density", e.staticcharge_density	);
		advsettings.addProperty("currentsource_mobility", e.currentsource_mobility	);
		advsettings.addProperty("max_EMF", e.max_EMF	);
		advsettings.addProperty("max_current", e.max_current	);
		advsettings.addProperty("default_AC_freq", e.default_AC_freq	);
		advsettings.addProperty("switch_open_mobility", e.switch_open_mobility	);

		advsettings.addProperty("a_factor_n", e.d_crit_n	);
		advsettings.addProperty("a_factor_p", e.d_crit_p	);
		advsettings.addProperty("eps_r_semi", e.eps_r_semi	);

		advsettings.add("modified_names", gson.toJsonTree(e.modified_names));
	}
	

	@SuppressWarnings("unchecked")
	public void readAdvancedSettings (Gson gson, JsonReader fstr, Preset p) throws JsonIOException, JsonSyntaxException, IOException, RuntimeException {
		
		if (p == null)
			e.setDefaultParameters();
		else
			p.applyPreset(e);

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
			case "gc_semi": e.gc_semi					= fstr.nextDouble(); break;
			case "gv_semi": e.gv_semi					= fstr.nextDouble(); break;
			case "chi_semi": e.chi_semi						= fstr.nextDouble(); break;
			case "Eg_semi": e.Eg_semi					= fstr.nextDouble(); break;
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
			
			case "g_metal": e.g_metal						= fstr.nextDouble(); break;
			case "g_metal_high": e.g_metal_high						= fstr.nextDouble(); break;
			case "g_metal_low": e.g_metal_low						= fstr.nextDouble(); break;
			case "W_metal_default": e.W_metal_default		= fstr.nextDouble(); break;
			case "Eg_metal": e.Eg_metal					= fstr.nextDouble(); break;
			case "W_metal_high": e.W_metal_high				= fstr.nextDouble(); break;
			case "W_metal_low": e.W_metal_low				= fstr.nextDouble(); break;
			case "k_rad_metal": e.k_rad_metal			= fstr.nextDouble(); break;

			case "mu_electron_metal": e.mu_electron_metal			= fstr.nextDouble(); break;
			case "mu_hole_metal": e.mu_hole_metal					= fstr.nextDouble(); break;
			case "v_sat_n_metal": e.v_sat_n_metal					= fstr.nextDouble(); break;
			case "v_sat_p_metal": e.v_sat_p_metal						= fstr.nextDouble(); break;
			case "eps_r_metal": e.eps_r_metal	= fstr.nextDouble(); break;

			case "dielectric_eps_r": e.dielectric_eps_r		= fstr.nextDouble(); break;
			case "ferromagnet_mu_r": e.ferromagnet_mu_r		= fstr.nextDouble(); break;
			case "staticcharge_density": e.staticcharge_density	= fstr.nextDouble(); break;
			case "currentsource_mobility": e.currentsource_mobility	= fstr.nextDouble(); break;
			case "max_EMF": e.max_EMF	= fstr.nextDouble(); break;
			case "max_current": e.max_current	= fstr.nextDouble(); break;
			case "default_AC_freq": e.default_AC_freq	= fstr.nextDouble(); break;
			case "switch_open_mobility": e.switch_open_mobility	= fstr.nextDouble(); break;

			case "a_factor_n": e.d_crit_n	= fstr.nextDouble(); break;
			case "a_factor_p": e.d_crit_p	= fstr.nextDouble(); break;
			case "eps_r_semi": e.eps_r_semi	= fstr.nextDouble(); break;

			case "modified_names": e.modified_names = (HashMap<MaterialType, String>) gson.fromJson(fstr, new TypeToken<HashMap<MaterialType, String>>(){}.getType()); break;

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
		} else if (ev.getSource() instanceof JTextField) {
			loadAdvancedSettings(true);
		}
	}
}
