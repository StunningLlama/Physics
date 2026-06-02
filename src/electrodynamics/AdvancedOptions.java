// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.JsonObject;
import com.google.gson.JsonSyntaxException;
import com.google.gson.stream.JsonReader;

import javax.swing.JTabbedPane;
import java.io.IOException;
import java.io.StringReader;

public class AdvancedOptions extends JFrame {

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
	public JTextField recomb_rate_semi;
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
		lblNewLabel_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_7.setBounds(45, 110, 168, 16);
		sim.add(lblNewLabel_7);
		
		junction_size = new JTextField();
		junction_size.setColumns(10);
		junction_size.setBounds(225, 105, 98, 26);
		sim.add(junction_size);
		
		JLabel lblNimetal = new JLabel("Intrinsic carrier conc. [1/m^3]");
		lblNimetal.setToolTipText("Metal equilibrium carrier concentration");
		lblNimetal.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNimetal.setBounds(35, 11, 203, 16);
		metal.add(lblNimetal);
		
		ni_metal = new JTextField();
		ni_metal.setColumns(10);
		ni_metal.setBounds(250, 6, 98, 26);
		metal.add(ni_metal);
		
		JLabel lblNewLabel_1_1 = new JLabel("Workfunction (default) [eV]");
		lblNewLabel_1_1.setToolTipText("Metal work function");
		lblNewLabel_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_1_1.setBounds(62, 77, 176, 16);
		metal.add(lblNewLabel_1_1);
		
		W_metal = new JTextField();
		W_metal.setColumns(10);
		W_metal.setBounds(250, 72, 98, 26);
		metal.add(W_metal);
		
		JLabel lblNewLabel_2_1 = new JLabel("Metal \"bandgap\" [eV]");
		lblNewLabel_2_1.setToolTipText("Metal \"band gap\"");
		lblNewLabel_2_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_2_1.setBounds(62, 44, 176, 16);
		metal.add(lblNewLabel_2_1);
		
		E_b_metal = new JTextField();
		E_b_metal.setColumns(10);
		E_b_metal.setBounds(250, 39, 98, 26);
		metal.add(E_b_metal);
		
		JLabel lblNewLabel_3_1 = new JLabel("Workfunction (High WF metal) [eV]");
		lblNewLabel_3_1.setToolTipText("Workfunction of low-workfunction metal");
		lblNewLabel_3_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3_1.setBounds(6, 110, 232, 16);
		metal.add(lblNewLabel_3_1);
		
		W_metal_high = new JTextField();
		W_metal_high.setColumns(10);
		W_metal_high.setBounds(250, 105, 98, 26);
		metal.add(W_metal_high);
		
		JLabel lblNewLabel_4_1 = new JLabel("Workfunction (Low WF metal) [eV]");
		lblNewLabel_4_1.setToolTipText("Workfunction of high-workfunction metal");
		lblNewLabel_4_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4_1.setBounds(16, 143, 222, 16);
		metal.add(lblNewLabel_4_1);
		
		W_metal_low = new JTextField();
		W_metal_low.setColumns(10);
		W_metal_low.setBounds(250, 138, 98, 26);
		metal.add(W_metal_low);
		
		JLabel lblNewLabel_6_1 = new JLabel("Recomb. rate [m^3/s]");
		lblNewLabel_6_1.setToolTipText("Metal carrier recombination rate (in number density units)");
		lblNewLabel_6_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6_1.setBounds(68, 176, 170, 16);
		metal.add(lblNewLabel_6_1);
		
		recomb_rate_metal = new JTextField();
		recomb_rate_metal.setColumns(10);
		recomb_rate_metal.setBounds(250, 171, 98, 26);
		metal.add(recomb_rate_metal);
		
		
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
		
		JLabel lblNewLabel_5_1 = new JLabel("Recomb. rate [m^3/s]");
		lblNewLabel_5_1.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_1.setBounds(36, 176, 168, 16);
		semi.add(lblNewLabel_5_1);
		
		recomb_rate_semi = new JTextField();
		recomb_rate_semi.setColumns(10);
		recomb_rate_semi.setBounds(216, 171, 98, 26);
		semi.add(recomb_rate_semi);
		
		lblNewLabel_5_2 = new JLabel("n-type default doping conc. [1/m^3]");
		lblNewLabel_5_2.setToolTipText("Doping concentration");
		lblNewLabel_5_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_2.setBounds(326, 11, 236, 16);
		semi.add(lblNewLabel_5_2);
		
		n_default_doping = new JTextField();
		n_default_doping.setColumns(10);
		n_default_doping.setBounds(574, 6, 98, 26);
		semi.add(n_default_doping);
		
		lblNewLabel_5_3 = new JLabel("p-type default doping conc. [1/m^3]");
		lblNewLabel_5_3.setToolTipText("Doping concentration");
		lblNewLabel_5_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_3.setBounds(326, 44, 236, 16);
		semi.add(lblNewLabel_5_3);
		
		p_default_doping = new JTextField();
		p_default_doping.setColumns(10);
		p_default_doping.setBounds(574, 39, 98, 26);
		semi.add(p_default_doping);
		
		lblNewLabel_5_4 = new JLabel("n-type light doping conc. [1/m^3]");
		lblNewLabel_5_4.setToolTipText("Doping concentration");
		lblNewLabel_5_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_4.setBounds(326, 77, 236, 16);
		semi.add(lblNewLabel_5_4);
		
		n_light_doping = new JTextField();
		n_light_doping.setColumns(10);
		n_light_doping.setBounds(574, 72, 98, 26);
		semi.add(n_light_doping);
		
		lblNewLabel_5_5 = new JLabel("p-type light doping conc. [1/m^3]");
		lblNewLabel_5_5.setToolTipText("Doping concentration");
		lblNewLabel_5_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_5.setBounds(326, 110, 236, 16);
		semi.add(lblNewLabel_5_5);
		
		p_light_doping = new JTextField();
		p_light_doping.setColumns(10);
		p_light_doping.setBounds(574, 105, 98, 26);
		semi.add(p_light_doping);
		
		lblNewLabel_5_6 = new JLabel("n-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_6.setToolTipText("Doping concentration");
		lblNewLabel_5_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_6.setBounds(326, 143, 236, 16);
		semi.add(lblNewLabel_5_6);
		
		n_heavy_doping = new JTextField();
		n_heavy_doping.setColumns(10);
		n_heavy_doping.setBounds(574, 138, 98, 26);
		semi.add(n_heavy_doping);
		
		lblNewLabel_5_7 = new JLabel("p-type heavy doping conc. [1/m^3]");
		lblNewLabel_5_7.setToolTipText("Doping concentration");
		lblNewLabel_5_7.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_7.setBounds(326, 176, 236, 16);
		semi.add(lblNewLabel_5_7);
		
		p_heavy_doping = new JTextField();
		p_heavy_doping.setColumns(10);
		p_heavy_doping.setBounds(574, 171, 98, 26);
		semi.add(p_heavy_doping);
		
		lblNewLabel_5_8 = new JLabel("Electron sat. velocity [m/s]");
		lblNewLabel_5_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_8.setBounds(12, 209, 192, 16);
		semi.add(lblNewLabel_5_8);
		
		v_sat_n = new JTextField();
		v_sat_n.setColumns(10);
		v_sat_n.setBounds(216, 204, 98, 26);
		semi.add(v_sat_n);
		
		lblNewLabel_5_9 = new JLabel("Hole sat. velocity [m/s]");
		lblNewLabel_5_9.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9.setBounds(12, 242, 192, 16);
		semi.add(lblNewLabel_5_9);
		
		v_sat_p = new JTextField();
		v_sat_p.setColumns(10);
		v_sat_p.setBounds(216, 237, 98, 26);
		semi.add(v_sat_p);
		
		other = new JPanel();
		tabbedPane.addTab("Other materials", null, other, null);
		other.setLayout(null);
		
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
		
		btn_apply = new JButton("Apply changes");
		btn_apply.setBounds(255, 372, 144, 29);
		panel.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(549, 372, 144, 29);
		panel.add(btn_cancel);
		
		btn_reset = new JButton("Reset to defaults");
		btn_reset.setBounds(403, 372, 144, 29);
		panel.add(btn_reset);
		
		//.add(tabbedPane);
	}
	
	public void initialize(Simulation e) {
		this.e = e;
		btn_apply.addActionListener(e.controls);
		btn_reset.addActionListener(e.controls);
		btn_cancel.addActionListener(e.controls);
		setVisible(false);
	}
	
	public void setInputsToDefault() {

		Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
		JsonObject advsettings = new JsonObject();
		writeAdvancedSettings(gson, advsettings);
		String json = gson.toJson(advsettings);
		
		e.setDefaultParameters();
		
		storeAdvancedSettings();
		e.adv_opts.width				.setText(formatDouble(e.default_width						));
		e.adv_opts.resolution			.setText(Integer.toString(e.default_resolution				));

		try {
			JsonReader fstr = new JsonReader(new StringReader(json));
			fstr.beginObject();
			this.readAdvancedSettings(gson, fstr);
			fstr.endObject();
		} catch (IOException | RuntimeException e1) {
			e1.printStackTrace();
		}
	}


	public String formatDouble(double d) {
		return Double.toString(d);
	}

	public void storeAdvancedSettings() {
		e.adv_opts.width				.setText(formatDouble(e.width						));
		e.adv_opts.resolution			.setText(Integer.toString(e.resolution				));
		e.adv_opts.depth				.setText(formatDouble(e.depth						));
		e.adv_opts.junction_size		.setText(Integer.toString(e.junction_size			));

		e.adv_opts.eps0					.setText(formatDouble(e.eps0						));
		e.adv_opts.mu0					.setText(formatDouble(e.mu0							));
		e.adv_opts.e_charge				.setText(formatDouble(e.e_charge					));
		e.adv_opts.T					.setText(formatDouble(e.T							));
		
		e.adv_opts.mu_electron			.setText(formatDouble(e.mu_electron					));
		e.adv_opts.mu_hole				.setText(formatDouble(e.mu_hole						));
		e.adv_opts.ni_semi				.setText(formatDouble(e.ni_semi						));
		e.adv_opts.W_semi				.setText(formatDouble(e.W_semi/e.eVtoJ				));
		e.adv_opts.E_b_semi				.setText(formatDouble(e.E_b_semi/e.eVtoJ			));
		e.adv_opts.recomb_rate_semi		.setText(formatDouble(e.recomb_rate_semi			));
		e.adv_opts.v_sat_n				.setText(formatDouble(e.v_sat_n			));
		e.adv_opts.v_sat_p				.setText(formatDouble(e.v_sat_p			));
		e.adv_opts.n_default_doping		.setText(formatDouble(e.n_default_doping_concentration		));
		e.adv_opts.p_default_doping		.setText(formatDouble(e.p_default_doping_concentration		));
		e.adv_opts.n_light_doping		.setText(formatDouble(e.n_light_doping_concentration		));
		e.adv_opts.p_light_doping		.setText(formatDouble(e.p_light_doping_concentration		));
		e.adv_opts.n_heavy_doping		.setText(formatDouble(e.n_heavy_doping_concentration		));
		e.adv_opts.p_heavy_doping		.setText(formatDouble(e.p_heavy_doping_concentration		));
		
		e.adv_opts.ni_metal				.setText(formatDouble(e.ni_metal					));
		e.adv_opts.W_metal				.setText(formatDouble(e.W_metal_default/e.eVtoJ		));
		e.adv_opts.E_b_metal			.setText(formatDouble(e.E_b_metal/e.eVtoJ			));
		e.adv_opts.W_metal_high			.setText(formatDouble(e.W_metal_high/e.eVtoJ		));
		e.adv_opts.W_metal_low			.setText(formatDouble(e.W_metal_low/e.eVtoJ			));
		e.adv_opts.recomb_rate_metal	.setText(formatDouble(e.recomb_rate_metal			));

		e.adv_opts.dielectric_eps_r			.setText(formatDouble(e.dielectric_eps_r			));
		e.adv_opts.ferromagnet_mu_r			.setText(formatDouble(e.ferromagnet_mu_r			));
		e.adv_opts.staticcharge_density		.setText(formatDouble(e.staticcharge_density			));
		e.adv_opts.currentsource_mobility	.setText(formatDouble(e.currentsource_mobility			));
	}
	
	public boolean loadAdvancedSettings(boolean show_warning) {

		try {
			double width_tmp			= Double.valueOf(e.adv_opts.width				.getText());	
			int resolution_tmp			= Integer.valueOf(e.adv_opts.resolution			.getText());

			if (resolution_tmp != e.resolution) {
				int result = JOptionPane.showConfirmDialog(this, "Changing the resolution will delete all materials. Proceed?", "Message", JOptionPane.YES_NO_OPTION);
				if (result != JOptionPane.OK_OPTION)
				{
					return false;
				}
			}

			e.depth						= Double.valueOf(e.adv_opts.depth				.getText());
			e.junction_size				= Integer.valueOf(e.adv_opts.junction_size		.getText());

			e.T							= Double.valueOf(e.adv_opts.T					.getText());
			e.eps0						= Double.valueOf(e.adv_opts.eps0				.getText());
			e.mu0						= Double.valueOf(e.adv_opts.mu0					.getText());
			e.e_charge					= Double.valueOf(e.adv_opts.e_charge			.getText());

			e.mu_electron				= Double.valueOf(e.adv_opts.mu_electron			.getText());
			e.mu_hole					= Double.valueOf(e.adv_opts.mu_hole				.getText());
			e.ni_semi					= Double.valueOf(e.adv_opts.ni_semi				.getText());
			e.W_semi					= Double.valueOf(e.adv_opts.W_semi				.getText())*e.eVtoJ;
			e.E_b_semi					= Double.valueOf(e.adv_opts.E_b_semi			.getText())*e.eVtoJ;
			e.recomb_rate_semi			= Double.valueOf(e.adv_opts.recomb_rate_semi	.getText());
			e.v_sat_n						= Double.valueOf(e.adv_opts.v_sat_n				.getText());
			e.v_sat_p						= Double.valueOf(e.adv_opts.v_sat_p				.getText());
			e.n_default_doping_concentration	= Double.valueOf(e.adv_opts.n_default_doping	.getText());
			e.p_default_doping_concentration	= Double.valueOf(e.adv_opts.p_default_doping	.getText());
			e.n_light_doping_concentration		= Double.valueOf(e.adv_opts.n_light_doping	.getText());
			e.p_light_doping_concentration		= Double.valueOf(e.adv_opts.p_light_doping	.getText());
			e.n_heavy_doping_concentration		= Double.valueOf(e.adv_opts.n_heavy_doping	.getText());
			e.p_heavy_doping_concentration		= Double.valueOf(e.adv_opts.p_heavy_doping	.getText());

			e.ni_metal					= Double.valueOf(e.adv_opts.ni_metal			.getText());
			e.W_metal_default			= Double.valueOf(e.adv_opts.W_metal				.getText())*e.eVtoJ;
			e.E_b_metal					= Double.valueOf(e.adv_opts.E_b_metal			.getText())*e.eVtoJ;
			e.W_metal_high				= Double.valueOf(e.adv_opts.W_metal_high		.getText())*e.eVtoJ;
			e.W_metal_low				= Double.valueOf(e.adv_opts.W_metal_low			.getText())*e.eVtoJ;
			e.recomb_rate_metal			= Double.valueOf(e.adv_opts.recomb_rate_metal	.getText());

			e.dielectric_eps_r			= Double.valueOf(e.adv_opts.dielectric_eps_r	.getText());
			e.ferromagnet_mu_r			= Double.valueOf(e.adv_opts.ferromagnet_mu_r	.getText());
			e.staticcharge_density		= Double.valueOf(e.adv_opts.staticcharge_density	.getText());
			e.currentsource_mobility	= Double.valueOf(e.adv_opts.currentsource_mobility	.getText());

			e.calculateDependentConstants();
			e.setSize(resolution_tmp, width_tmp);
			e.resetFields(false);
			e.lastsimspeed = -1;
			e.advsettings_tweaked = true;
			return true;
			
		} catch (NumberFormatException e) {
			JOptionPane.showMessageDialog(this, "Invalid number.", "Error", JOptionPane.OK_OPTION);
			return false;
		}
	}

	public void writeAdvancedSettings (Gson gson, JsonObject advsettings) {
		advsettings.addProperty("depth", e.depth 						);
		advsettings.addProperty("junction_size", e.junction_size		);

		advsettings.addProperty("T", e.T								);
		advsettings.addProperty("eps0", e.eps0							);
		advsettings.addProperty("mu0", e.mu0							);
		advsettings.addProperty("e_charge", e.e_charge					);
		
		advsettings.addProperty("mu_electron", e.mu_electron			);
		advsettings.addProperty("mu_hole", e.mu_hole					);
		advsettings.addProperty("ni_semi", e.ni_semi					);
		advsettings.addProperty("W_semi", e.W_semi						);
		advsettings.addProperty("E_b_semi", e.E_b_semi					);
		advsettings.addProperty("recomb_rate_semi", e.recomb_rate_semi	);
		advsettings.addProperty("v_sat_n", e.v_sat_n						);
		advsettings.addProperty("v_sat_p", e.v_sat_p						);
		advsettings.addProperty("n_default_doping", e.n_default_doping_concentration	);
		advsettings.addProperty("p_default_doping", e.p_default_doping_concentration	);
		advsettings.addProperty("n_light_doping", e.n_light_doping_concentration	);
		advsettings.addProperty("p_light_doping", e.p_light_doping_concentration	);
		advsettings.addProperty("n_heavy_doping", e.n_heavy_doping_concentration	);
		advsettings.addProperty("p_heavy_doping", e.p_heavy_doping_concentration	);
		
		advsettings.addProperty("ni_metal", e.ni_metal						);
		advsettings.addProperty("W_metal_default", e.W_metal_default		);
		advsettings.addProperty("E_b_metal", e.E_b_metal					);
		advsettings.addProperty("W_metal_high", e.W_metal_high				);
		advsettings.addProperty("W_metal_low", e.W_metal_low				);
		advsettings.addProperty("recomb_rate_metal", e.recomb_rate_metal	);

		advsettings.addProperty("dielectric_eps_r", e.dielectric_eps_r	);
		advsettings.addProperty("ferromagnet_mu_r", e.ferromagnet_mu_r	);
		advsettings.addProperty("staticcharge_density", e.staticcharge_density	);
		advsettings.addProperty("currentsource_mobility", e.currentsource_mobility	);
	}
	

	public void readAdvancedSettings (Gson gson, JsonReader fstr) throws JsonIOException, JsonSyntaxException, IOException, RuntimeException {
		
		e.setDefaultParameters();

		while (fstr.hasNext()) {
			String name = fstr.nextName();
			switch (name){
			case "depth": e.depth 						= fstr.nextDouble(); break;
			case "junction_size": e.junction_size 		= fstr.nextInt(); break;

			case "T": e.T								= fstr.nextDouble(); break;
			case "eps0": e.eps0							= fstr.nextDouble(); break;
			case "mu0": e.mu0							= fstr.nextDouble(); break;
			case "e_charge": e.e_charge					= fstr.nextDouble(); break;
			
			case "mu_electron": e.mu_electron			= fstr.nextDouble(); break;
			case "mu_hole": e.mu_hole					= fstr.nextDouble(); break;
			case "ni_semi": e.ni_semi					= fstr.nextDouble(); break;
			case "W_semi": e.W_semi						= fstr.nextDouble(); break;
			case "E_b_semi": e.E_b_semi					= fstr.nextDouble(); break;
			case "recomb_rate_semi": e.recomb_rate_semi		= fstr.nextDouble(); break;
			case "v_sat_n": e.v_sat_n						= fstr.nextDouble(); break;
			case "v_sat_p": e.v_sat_p						= fstr.nextDouble(); break;
			case "n_default_doping": e.n_default_doping_concentration	= fstr.nextDouble(); break;
			case "p_default_doping": e.p_default_doping_concentration	= fstr.nextDouble(); break;
			case "n_light_doping": e.n_light_doping_concentration		= fstr.nextDouble(); break;
			case "p_light_doping": e.p_light_doping_concentration		= fstr.nextDouble(); break;
			case "n_heavy_doping": e.n_heavy_doping_concentration		= fstr.nextDouble(); break;
			case "p_heavy_doping": e.p_heavy_doping_concentration		= fstr.nextDouble(); break;
			
			case "ni_metal": e.ni_metal						= fstr.nextDouble(); break;
			case "W_metal_default": e.W_metal_default		= fstr.nextDouble(); break;
			case "E_b_metal": e.E_b_metal					= fstr.nextDouble(); break;
			case "W_metal_high": e.W_metal_high				= fstr.nextDouble(); break;
			case "W_metal_low": e.W_metal_low				= fstr.nextDouble(); break;
			case "recomb_rate_metal": e.recomb_rate_metal	= fstr.nextDouble(); break;

			case "dielectric_eps_r": e.dielectric_eps_r		= fstr.nextDouble(); break;
			case "ferromagnet_mu_r": e.ferromagnet_mu_r		= fstr.nextDouble(); break;
			case "staticcharge_density": e.staticcharge_density	= fstr.nextDouble(); break;
			case "currentsource_mobility": e.currentsource_mobility	= fstr.nextDouble(); break;

			default: fstr.skipValue(); break; // skip others
			}
		}
	}
}
