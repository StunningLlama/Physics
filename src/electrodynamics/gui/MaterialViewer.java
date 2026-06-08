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

	public MaterialViewer() {
		setTitle("Material property viewer");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 705, 562);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 6, 235, 522);
		contentPane.add(scrollPane);
		
		list = new JList<>();
		scrollPane.setViewportView(list);
		
		btn_cancel = new JButton("Close");
		btn_cancel.setBounds(575, 499, 124, 29);
		contentPane.add(btn_cancel);
		
		btn_refresh = new JButton("Refresh");
		btn_refresh.setBounds(452, 499, 124, 29);
		contentPane.add(btn_refresh);
		
		JScrollPane scrollPane_1 = new JScrollPane();
		scrollPane_1.setBounds(253, 6, 435, 480);
		contentPane.add(scrollPane_1);
		
		textPane = new JTextPane();
		scrollPane_1.setViewportView(textPane);
	}
	
	public void initialize(Simulation e) {
		this.e = e;
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
			

			double rho_n = e.calcEquilibriumElectronCharge(mat.rho_back, mat.ni*mat.ni);
			double rho_p = e.calcEquilibriumHoleCharge(mat.rho_back, mat.ni*mat.ni);
			double sigma = e.e_charge*(-mat.D_n*e.beta*rho_n + mat.D_p*e.beta*rho_p);
			
			String str = "<html>";
			str += "Name: " + mat.toString() + "<br>";
			str += "<br>";
			str += "Dielectric const.: <b>" + e.units.toString(mat.eps_r, Quantity.DIMENSIONLESS) + "</b><br>";
			str += "Relative permeability: <b>" + e.units.toString(mat.mu_r, Quantity.DIMENSIONLESS) + "</b><br>";
			str += "Speed of light in material: <b>" + e.units.toString(1/Math.sqrt(e.eps0*mat.eps_r * e.mu0*mat.mu_r), Quantity.VELOCITY) + "</b><br>";
			str += "<br>";
			str += "Electron affinity: <b>" + e.units.toString((mat.W - 0.5*mat.Eb)/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Bandgap: <b>" + e.units.toString(mat.Eb/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Workfunction: <b>" + e.units.toString(mat.W/e.eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "</b><br>";
			str += "Eq. carrier concentration (ni): <b>" + e.units.toString(mat.ni, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Eq. electron density: <b>" + e.units.toString(-rho_n/e.e_charge, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Eq. hole density: <b>" + e.units.toString(rho_p/e.e_charge, Quantity.NUMBER_DENSITY) + "</b><br>";
			str += "Conductivity: <b>" + e.units.toString(sigma, Quantity.CONDUCTIVITY) + "</b><br>";
			str += "<br>";
			str += "Electron mobility: <b>" + e.units.toString(mat.D_n*e.beta*e.e_charge, Quantity.ELECTRIC_MOBILITY) + "</b><br>";
			str += "Electron diffusivity: <b>" + e.units.toString(mat.D_n, Quantity.DIFFUSIVITY) + "</b><br>";
			str += "Electron saturation velocity: <b>" + e.units.toString(mat.v_sat_n, Quantity.VELOCITY) + "</b><br>";
			str += "Electron saturation field: <b>" + e.units.toString(mat.v_sat_n/(mat.D_n*e.beta*e.e_charge), Quantity.ELECTRIC_FIELD) + "</b><br>";
			str += "Hole mobility: <b>" + e.units.toString(mat.D_p*e.beta*e.e_charge, Quantity.ELECTRIC_MOBILITY) + "</b><br>";
			str += "Hole diffusivity: <b>" + e.units.toString(mat.D_p, Quantity.DIFFUSIVITY) + "</b><br>";
			str += "Hole saturation velocity: <b>" + e.units.toString(mat.v_sat_p, Quantity.VELOCITY) + "</b><br>";
			str += "Hole saturation field: <b>" + e.units.toString(mat.v_sat_p/(mat.D_p*e.beta*e.e_charge), Quantity.ELECTRIC_FIELD) + "</b><br>";
			str += "</html>";
			textPane.setText(str);
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
