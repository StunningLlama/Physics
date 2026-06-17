package electrodynamics.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import electrodynamics.GeneralMaterialType;
import electrodynamics.Material;
import electrodynamics.MaterialType;
import electrodynamics.Simulation;
import electrodynamics.units.Quantity;
import java.awt.BorderLayout;
import javax.swing.BoxLayout;

public class MaterialViewer extends JFrame implements ActionListener, ListSelectionListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	private JPanel contentPane;
	public JTextField width;
	public JTextField resolution;
	public JTextField depth;
	public JTextField ni_metal;
	public JTextField W_metal;
	public JTextField E_b_metal;
	public JTextField W_metal_high;
	public JTextField W_metal_low;
	public JTextField recomb_rate_metal;
	public JTextField T;
	
	private JList<GeneralMaterialType> list;
	private JButton btn_cancel;
	private JButton btn_refresh;
	private JTextPane textPane;
	private JPanel panel;
	private JPanel panel_1;
	private JPanel panel_2;
	private JPanel panel_3;
	private JScrollPane scrollPane_1;

	public MaterialViewer(Simulation e) {
		this.e = e;
		setTitle("Material property viewer");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 705, 595);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new BorderLayout(0, 0));
		
		panel = new JPanel();
		panel.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel, BorderLayout.WEST);
		panel.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane = new JScrollPane();
		panel.add(scrollPane);
		
		list = new JList<>();
		scrollPane.setViewportView(list);
		
		panel_1 = new JPanel();
		panel_1.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel_1, BorderLayout.CENTER);
		panel_1.setLayout(new BorderLayout(0, 0));
		
		scrollPane_1 = new JScrollPane();
		panel_1.add(scrollPane_1);
		
		textPane = new JTextPane();
		scrollPane_1.setViewportView(textPane);
		
		panel_3 = new JPanel();
		contentPane.add(panel_3, BorderLayout.SOUTH);
		panel_3.setLayout(new BorderLayout(0, 0));
		
		panel_2 = new JPanel();
		panel_3.add(panel_2, BorderLayout.EAST);
		panel_2.setLayout(new BoxLayout(panel_2, BoxLayout.X_AXIS));
		
		btn_refresh = new JButton("Refresh");
		panel_2.add(btn_refresh);
		
		btn_cancel = new JButton("Close");
		panel_2.add(btn_cancel);
	}
	
	public void initialize() {
		this.btn_cancel.addActionListener(this);
		this.btn_refresh.addActionListener(this);
		this.list.addListSelectionListener(this);
		setLocationRelativeTo(null);

		textPane.setContentType("text/html");
		updateUI();
	}
	
	public void updateUI() {
		list.setModel(new DefaultComboBoxModel<>(makelist()));
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == btn_cancel) {
			this.setVisible(false);
		} else if (e.getSource() == btn_refresh) {
			GeneralMaterialType type = list.getSelectedValue();
			updateUI();
			list.setSelectedValue(type, true);
		}
	}

	@Override
	public void valueChanged(ListSelectionEvent ev) {
		GeneralMaterialType type = list.getSelectedValue();
		if (type != null) {
			Material mat = new Material();
			e.initializeMaterial(mat, type);
			
			double rho_n = mat.calcEquilibriumElectronCharge(e.e_charge, e.k*e.T);
			double rho_p = mat.calcEquilibriumHoleCharge(e.e_charge, e.k*e.T);
			double sigma = mat.calcConductivity(e.e_charge, e.k*e.T);
			

			double n = rho_n/e.q_n;
			double p = rho_p/e.q_p;
			double ni = Math.sqrt(n*p);
			double rate_const = (mat.k_aug_n*n+mat.k_aug_p*p)
				+ (mat.k_SRH_n*mat.k_SRH_p)/(mat.k_SRH_n*(n+ni) + mat.k_SRH_p*(p+ni) + Double.MIN_VALUE)
				+ mat.k_rad;
			double n_time = 1/(rate_const*p);
			double p_time = 1/(rate_const*n);
			
			String str = "<html>";
			str += "Name: " + mat.getDisplayName() + "<br>";
			str += "<br>";
			str += "Dielectric const.: <b>" + e.units.toString(mat.eps_r, Quantity.DIMENSIONLESS) + "</b><br>";
			str += "Relative permeability: <b>" + e.units.toString(mat.mu_r, Quantity.DIMENSIONLESS) + "</b><br>";
			str += "Speed of light in material: <b>" + e.units.toString(1/Math.sqrt(e.eps0*mat.eps_r * e.mu0*mat.mu_r), Quantity.VELOCITY) + "</b><br>";
			str += "<br>";
			str += "Electron affinity: <b>" + e.units.toString(-mat.Ec/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Conduction band energy: <b>" + e.units.toString(mat.Ec/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Valence band energy: <b>" + e.units.toString(mat.Ev/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Band gap: <b>" + e.units.toString((mat.Ec - mat.Ev)/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Workfunction: <b>" + e.units.toString(mat.calcPhi(e.k*e.T)/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Effective conduction band DOS: <b>" + e.units.toString(mat.gc, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Effective valence band DOS: <b>" + e.units.toString(mat.gv, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Dopant concentration: <b>" + e.units.toString(Math.abs(mat.rho_back/e.e_charge), Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Eq. carrier concentration (ni): <b>" + e.units.toString(mat.calc_ni(e.k*e.T) , Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Eq. electron density: <b>" + e.units.toString(-rho_n/e.e_charge, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Eq. hole density: <b>" + e.units.toString(rho_p/e.e_charge, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "<br>";
			str += "Electron mobility: <b>" + e.units.toString(mat.D_n*e.beta*e.e_charge, Quantity.ELECTRIC_MOBILITY) + "</b><br>";
			str += "Electron diffusivity: <b>" + e.units.toString(mat.D_n, Quantity.DIFFUSIVITY) + "</b><br>";
			str += "Electron saturation velocity: <b>" + e.units.toString(mat.v_sat_n, Quantity.VELOCITY) + "</b><br>";
			str += "Electron saturation field: <b>" + e.units.toString(mat.v_sat_n/(mat.D_n*e.beta*e.e_charge), Quantity.ELECTRIC_FIELD) + "</b><br>";
			str += "Hole mobility: <b>" + e.units.toString(mat.D_p*e.beta*e.e_charge, Quantity.ELECTRIC_MOBILITY) + "</b><br>";
			str += "Hole diffusivity: <b>" + e.units.toString(mat.D_p, Quantity.DIFFUSIVITY) + "</b><br>";
			str += "Hole saturation velocity: <b>" + e.units.toString(mat.v_sat_p, Quantity.VELOCITY) + "</b><br>";
			str += "Hole saturation field: <b>" + e.units.toString(mat.v_sat_p/(mat.D_p*e.beta*e.e_charge), Quantity.ELECTRIC_FIELD) + "</b><br>";
			str += "<br>";
			str += "Radiative recombination rate coefficient: <b>" + e.units.toString(mat.k_rad, Quantity.RATE_DENSITY) + "</b><br>";
			str += "SRH rate coefficient n: <b>" + e.units.toString(mat.k_SRH_n, Quantity.RATE) + "</b><br>";
			str += "SRH rate coefficient p: <b>" + e.units.toString(mat.k_SRH_p, Quantity.RATE) + "</b><br>";
			str += "Auger rate coefficient n: <b>" + e.units.toString(mat.k_aug_n, Quantity.RATE_DENSITY_SQUARED) + "</b><br>";
			str += "Auger rate coefficient p: <b>" + e.units.toString(mat.k_aug_p, Quantity.RATE_DENSITY_SQUARED) + "</b><br>";
			str += "<br>";
			str += "Excess electron lifetime: <b>" + e.units.toString(n_time, Quantity.TIME) + "</b><br>";
			str += "Electron diffusion length: <b>" + e.units.toString(Math.sqrt(n_time*mat.D_n), Quantity.LENGTH) + "</b><br>";
			str += "Excess hole lifetime: <b>" + e.units.toString(p_time, Quantity.TIME) + "</b><br>";
			str += "Hole diffusion length: <b>" + e.units.toString(Math.sqrt(p_time*mat.D_p), Quantity.LENGTH) + "</b><br>";
			str += "Conductivity: <b>" + e.units.toString(sigma, Quantity.CONDUCTIVITY) + "</b><br>";
			str += "Resistivity: <b>" + e.units.toString(1/sigma, Quantity.RESISTIVITY) + "</b><br>";
			str += "</html>";
			textPane.setText(str);
			textPane.setCaretPosition(0);
			//this.scrollPane_1.getVerticalScrollBar().setValue(0);
		}
	}
	
	public GeneralMaterialType[] makelist() {
		ArrayList<GeneralMaterialType> matlist = new ArrayList<GeneralMaterialType>();
		for (MaterialType t : MaterialType.values()) {
			if (t != MaterialType.CUSTOM)
				matlist.add(new GeneralMaterialType(t));
		}
		
		for (int i : e.materialmanager.mat_map.keySet()) {
			matlist.add(new GeneralMaterialType(i, e.materialmanager.mat_map.get(i).name));
		}
		
		return matlist.toArray(new GeneralMaterialType[0]);
	}
}
