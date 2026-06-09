package electrodynamics.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import electrodynamics.GeneralMaterialType;
import electrodynamics.Material;
import electrodynamics.MaterialType;
import electrodynamics.Simulation;
import electrodynamics.util.Utils;

public class MaterialManager extends JFrame implements ActionListener, ListSelectionListener, ItemListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	private JPanel contentPane;
	public JTextField width;
	public JTextField resolution;
	public JTextField depth;
	public JTextField mu_electron;
	public JTextField mu_hole;
	public JTextField ni;
	public JTextField W;
	public JTextField Eb;
	public JTextField ni_metal;
	public JTextField W_metal;
	public JTextField E_b_metal;
	public JTextField W_metal_high;
	public JTextField W_metal_low;
	public JTextField k_rad;
	public JTextField recomb_rate_metal;
	public JTextField T;
	private JLabel lbl_v_sat_n;
	private JTextField v_sat_n;
	private JLabel lbl_v_sat_p;
	private JTextField v_sat_p;
	private JLabel lbl_k_SRH_n;
	private JTextField k_SRH_n;
	private JLabel lbl_k_SRH_p;
	private JTextField k_SRH_p;
	private JLabel lbl_k_aug_n;
	private JTextField k_aug_n;
	private JLabel lbl_k_aug_p;
	private JTextField k_aug_p;
	private JTextField eps_r;
	private JLabel lbl_mu_r;
	private JTextField mu_r;
	private JLabel lbl_rho_back;
	private JTextField rho_back;
	private JLabel lblNewLabel;
	private JTextField name;
	private JComboBox<MaterialClass> type;
	private JLabel lblType;
	private JButton btn_delete;
	private JButton btn_add;
	
	public HashMap<Integer, Material> mat_map = new HashMap<Integer, Material>();
	public int id_counter = 0;
	private JButton btn_apply;
	private JList<Material> list;
	private JLabel lbl_mu_electron;
	private JLabel lbl_mu_hole;
	private JLabel lbl_ni;
	private JLabel lbl_W;
	private JLabel lbl_Eb;
	private JLabel lbl_eps_r;
	private JLabel lbl_k_rad;
	private JButton btn_cancel;

	public MaterialManager(Simulation e) {
		this.e = e;
		setResizable(false);
		setTitle("Material editor");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 905, 363);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 6, 235, 283);
		contentPane.add(scrollPane);
		
		list = new JList<>();
		scrollPane.setViewportView(list);
		
		btn_add = new JButton("New material");
		btn_add.setBounds(6, 301, 117, 29);
		contentPane.add(btn_add);
		
		btn_delete = new JButton("Delete");
		btn_delete.setBounds(124, 301, 117, 29);
		contentPane.add(btn_delete);

		mu_electron = new JTextField();
		mu_electron.setColumns(10);
		mu_electron.setBounds(477, 66, 98, 26);
		contentPane.add(mu_electron);
		
		lbl_mu_electron = new JLabel("Electron mobility [m^2/(V s)]");
		lbl_mu_electron.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_electron.setBounds(273, 71, 192, 16);
		contentPane.add(lbl_mu_electron);
		
		lbl_mu_hole = new JLabel("Hole mobility [m^2/(V s)]");
		lbl_mu_hole.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_hole.setBounds(297, 104, 168, 16);
		contentPane.add(lbl_mu_hole);
		
		mu_hole = new JTextField();
		mu_hole.setColumns(10);
		mu_hole.setBounds(477, 99, 98, 26);
		contentPane.add(mu_hole);
		
		lbl_ni = new JLabel("Intrinsic carrier conc. [1/m^3]");
		lbl_ni.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_ni.setBounds(267, 137, 198, 16);
		contentPane.add(lbl_ni);
		
		ni = new JTextField();
		ni.setColumns(10);
		ni.setBounds(477, 132, 98, 26);
		contentPane.add(ni);
		
		lbl_W = new JLabel("Workfunction [eV]");
		lbl_W.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_W.setBounds(297, 170, 168, 16);
		contentPane.add(lbl_W);
		
		W = new JTextField();
		W.setColumns(10);
		W.setBounds(477, 165, 98, 26);
		contentPane.add(W);
		
		lbl_Eb = new JLabel("Bandgap [eV]");
		lbl_Eb.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_Eb.setBounds(297, 203, 168, 16);
		contentPane.add(lbl_Eb);
		
		Eb = new JTextField();
		Eb.setColumns(10);
		Eb.setBounds(477, 198, 98, 26);
		contentPane.add(Eb);
		
		lbl_k_rad = new JLabel("Radiative recomb. rate [m^3/s]");
		lbl_k_rad.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_rad.setBounds(578, 5, 211, 16);
		contentPane.add(lbl_k_rad);
		
		k_rad = new JTextField();
		k_rad.setColumns(10);
		k_rad.setBounds(801, 0, 98, 26);
		contentPane.add(k_rad);
		
		lbl_v_sat_n = new JLabel("Electron sat. velocity [m/s]");
		lbl_v_sat_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_v_sat_n.setBounds(597, 170, 192, 16);
		contentPane.add(lbl_v_sat_n);
		
		v_sat_n = new JTextField();
		v_sat_n.setColumns(10);
		v_sat_n.setBounds(801, 165, 98, 26);
		contentPane.add(v_sat_n);
		
		lbl_v_sat_p = new JLabel("Hole sat. velocity [m/s]");
		lbl_v_sat_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_v_sat_p.setBounds(597, 203, 192, 16);
		contentPane.add(lbl_v_sat_p);
		
		v_sat_p = new JTextField();
		v_sat_p.setColumns(10);
		v_sat_p.setBounds(801, 198, 98, 26);
		contentPane.add(v_sat_p);
		
		lbl_k_SRH_n = new JLabel("SRH recomb. rate n [1/s]");
		lbl_k_SRH_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_SRH_n.setBounds(597, 38, 192, 16);
		contentPane.add(lbl_k_SRH_n);
		
		k_SRH_n = new JTextField();
		k_SRH_n.setColumns(10);
		k_SRH_n.setBounds(801, 33, 98, 26);
		contentPane.add(k_SRH_n);
		
		lbl_k_SRH_p = new JLabel("SRH recomb. rate p [1/s]");
		lbl_k_SRH_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_SRH_p.setBounds(607, 71, 182, 16);
		contentPane.add(lbl_k_SRH_p);
		
		k_SRH_p = new JTextField();
		k_SRH_p.setColumns(10);
		k_SRH_p.setBounds(801, 66, 98, 26);
		contentPane.add(k_SRH_p);
		
		lbl_k_aug_n = new JLabel("Auger recomb. rate n [m^6/s]");
		lbl_k_aug_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_aug_n.setBounds(597, 104, 192, 16);
		contentPane.add(lbl_k_aug_n);
		
		k_aug_n = new JTextField();
		k_aug_n.setColumns(10);
		k_aug_n.setBounds(801, 99, 98, 26);
		contentPane.add(k_aug_n);
		
		lbl_k_aug_p = new JLabel("Auger recomb. rate p [m^6/s]");
		lbl_k_aug_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_aug_p.setBounds(597, 137, 192, 16);
		contentPane.add(lbl_k_aug_p);
		
		k_aug_p = new JTextField();
		k_aug_p.setColumns(10);
		k_aug_p.setBounds(801, 132, 98, 26);
		contentPane.add(k_aug_p);
		
		lbl_eps_r = new JLabel("Dielectric constant");
		lbl_eps_r.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_eps_r.setBounds(297, 235, 168, 16);
		contentPane.add(lbl_eps_r);
		
		eps_r = new JTextField();
		eps_r.setColumns(10);
		eps_r.setBounds(477, 230, 98, 26);
		contentPane.add(eps_r);
		
		lbl_mu_r = new JLabel("Rel. permeability");
		lbl_mu_r.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_r.setBounds(297, 268, 168, 16);
		contentPane.add(lbl_mu_r);
		
		mu_r = new JTextField();
		mu_r.setColumns(10);
		mu_r.setBounds(477, 263, 98, 26);
		contentPane.add(mu_r);
		
		lbl_rho_back = new JLabel("Dopant charge density [C/m^3]");
		lbl_rho_back.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_rho_back.setBounds(578, 236, 211, 16);
		contentPane.add(lbl_rho_back);
		
		rho_back = new JTextField();
		rho_back.setColumns(10);
		rho_back.setBounds(801, 231, 98, 26);
		contentPane.add(rho_back);
		
		lblNewLabel = new JLabel("Name");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(273, 6, 88, 16);
		contentPane.add(lblNewLabel);
		
		name = new JTextField();
		name.setColumns(10);
		name.setBounds(373, 0, 202, 26);
		contentPane.add(name);
		
		type = new JComboBox<>();
		type.setModel(new DefaultComboBoxModel<>(MaterialClass.values()));
		type.setBounds(373, 34, 202, 27);
		contentPane.add(type);
		
		lblType = new JLabel("Type");
		lblType.setHorizontalAlignment(SwingConstants.TRAILING);
		lblType.setBounds(273, 38, 88, 16);
		contentPane.add(lblType);
		
		btn_apply = new JButton("Apply");
		btn_apply.setBounds(651, 301, 124, 29);
		contentPane.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(775, 301, 124, 29);
		contentPane.add(btn_cancel);
	}
	
	public void initialize() {
		this.btn_add.addActionListener(this);
		this.btn_delete.addActionListener(this);
		this.btn_apply.addActionListener(this);
		//this.btn_save.addActionListener(this);
		this.btn_cancel.addActionListener(this);
		this.list.addListSelectionListener(this);
		this.type.addItemListener(this);
		e.opts.gui_material.setModel(new DefaultComboBoxModel<GeneralMaterialType>(makelist()));
		setInputVisibility();
		setLocationRelativeTo(null);
	}
	
	public void resetMaterialList() {
		mat_map.clear();
		id_counter = 0;
	}
	
	public void addmat() {
		Material mat = new Material();
		
		GeneralMaterialType[] arr = makelist();
		
		int ind = Arrays.asList(arr).indexOf(new GeneralMaterialType(MaterialType.VACUUM));

        JList<GeneralMaterialType> tmplist = new JList<>(arr);

        // Optional settings
        tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        tmplist.setVisibleRowCount(5);

        JScrollPane scrollPane = new JScrollPane(tmplist);

        tmplist.setSelectedValue(arr[ind], true);

        int result = JOptionPane.showConfirmDialog(
                null,
                scrollPane,
                "Select template",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE
        );

        if (result == JOptionPane.OK_OPTION) {
        	GeneralMaterialType selected = tmplist.getSelectedValue();

            if (selected != null) {
            	if (selected.cust_id == -1)
            		e.initializeMaterial(mat, selected.type);
            	else
                    mat.copyFrom(mat_map.get(selected.cust_id));
            }
        } else {
        	return;
        }
		
        mat.type = MaterialType.CUSTOM;
		mat.name = ("New material " + id_counter);
		mat.cust_id = id_counter;
		mat_map.put(mat.cust_id, mat);
		id_counter++;
		updateUI();
		e.initializeAllMaterials();
		e.updateAllMaterials(true);
		list.setSelectedValue(mat, true);
	}
	
	public void delete() {
		Material mat = list.getSelectedValue();
		if (mat != null) {
			mat_map.remove(mat.cust_id);
			updateUI();
			e.initializeAllMaterials();
			e.updateAllMaterials(true);
		}
	}

	public void save() {
		Material mat = list.getSelectedValue();
		if (mat != null) {
			loadMaterial(mat);
			updateUI();
			e.initializeAllMaterials();
			e.updateAllMaterials(true);
			list.setSelectedValue(mat, true);
		}
	}
	
	public void updateUI() {
		list.setModel(new DefaultComboBoxModel<Material>(mat_map.values().toArray(new Material[0])));
		e.opts.gui_material.setModel(new DefaultComboBoxModel<GeneralMaterialType>(makelist()));
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == btn_apply) {
			save();
		} else if (e.getSource() == btn_add) {
			addmat();
		} else if (e.getSource() == btn_delete) {
			delete();
		} else if (e.getSource() == btn_cancel) {
			this.setVisible(false);
		}
	}

	@Override
	public void valueChanged(ListSelectionEvent e) {
		Material mat = list.getSelectedValue();
		if (mat != null) {
			storeMaterial(mat);
		}
	}
	
	public void storeMaterial(Material mat) {
		if (mat.semiconducting == 1)
			type.setSelectedItem(MaterialClass.SEMICONDUCTING);
		else if (mat.conducting == 1)
			type.setSelectedItem(MaterialClass.CONDUCTING);
		else
			type.setSelectedItem(MaterialClass.INSULATING);
		
		name.setText(mat.name);
		eps_r.setText(Utils.formatDouble(mat.eps_r));
		mu_r.setText(Utils.formatDouble(mat.mu_r));
		rho_back.setText(Utils.formatDouble(mat.rho_back));
		ni.setText(Utils.formatDouble(mat.ni));
		W.setText(Utils.formatDouble(mat.W/e.eVtoJ));
		Eb.setText(Utils.formatDouble(mat.Eb/e.eVtoJ));
		mu_electron.setText(Utils.formatDouble(mat.D_n*e.beta*e.e_charge));
		mu_hole.setText(Utils.formatDouble(mat.D_p*e.beta*e.e_charge));
		v_sat_n.setText(Utils.formatDouble(mat.v_sat_n));
		v_sat_p.setText(Utils.formatDouble(mat.v_sat_p));
		k_rad.setText(Utils.formatDouble(mat.k_rad));
		k_aug_n.setText(Utils.formatDouble(mat.k_aug_n));
		k_aug_p.setText(Utils.formatDouble(mat.k_aug_p));
		k_SRH_n.setText(Utils.formatDouble(mat.k_SRH_n));
		k_SRH_p.setText(Utils.formatDouble(mat.k_SRH_p));
	}
	
	public void loadMaterial(Material mat) {
		mat.setDefaultParameters();
		
		MaterialClass matclass = (MaterialClass) type.getSelectedItem();
		
		if (matclass == MaterialClass.CONDUCTING || matclass == MaterialClass.SEMICONDUCTING)
			mat.conducting = 1;

		if (matclass == MaterialClass.SEMICONDUCTING)
			mat.semiconducting = 1;
		
		mat.name = name.getText();
		
		mat.eps_r = Double.valueOf(eps_r.getText());
		mat.mu_r = Double.valueOf(mu_r.getText());
		mat.rho_back = Double.valueOf(rho_back.getText());

		if (mat.conducting == 1) {
			mat.ni = Double.valueOf(ni.getText());
			mat.W = Double.valueOf(W.getText())*e.eVtoJ;
			mat.Eb = Double.valueOf(Eb.getText())*e.eVtoJ;
			mat.D_n = Double.valueOf(mu_electron.getText())/(e.beta*e.e_charge);
			mat.D_p = Double.valueOf(mu_hole.getText())/(e.beta*e.e_charge);
			mat.v_sat_n = Double.valueOf(v_sat_n.getText());
			mat.v_sat_p = Double.valueOf(v_sat_p.getText());
			mat.k_rad = Double.valueOf(k_rad.getText());
		}

		if (mat.semiconducting == 1) {
			mat.k_aug_n = Double.valueOf(k_aug_n.getText());
			mat.k_aug_p = Double.valueOf(k_aug_p.getText());
			mat.k_SRH_n = Double.valueOf(k_SRH_n.getText());
			mat.k_SRH_p = Double.valueOf(k_SRH_p.getText());
		}
	}
	
	public void setInputVisibility() {
		MaterialClass matclass = (MaterialClass) type.getSelectedItem();
		
		boolean conducting = (matclass == MaterialClass.CONDUCTING || matclass == MaterialClass.SEMICONDUCTING);
		boolean semiconducting = (matclass == MaterialClass.SEMICONDUCTING);
		
		ni.setVisible(conducting);
		W.setVisible(conducting);
		ni.setVisible(conducting);
		Eb.setVisible(conducting);
		mu_electron.setVisible(conducting);
		mu_hole.setVisible(conducting);
		k_rad.setVisible(conducting);
		
		v_sat_n.setVisible(conducting);
		v_sat_p.setVisible(conducting);
		k_aug_n.setVisible(semiconducting);
		k_aug_p.setVisible(semiconducting);
		k_SRH_n.setVisible(semiconducting);
		k_SRH_p.setVisible(semiconducting);
		
		lbl_ni.setEnabled(conducting);
		lbl_W.setEnabled(conducting);
		lbl_ni.setEnabled(conducting);
		lbl_Eb.setEnabled(conducting);
		lbl_mu_electron.setEnabled(conducting);
		lbl_mu_hole.setEnabled(conducting);
		lbl_k_rad.setEnabled(conducting);
		
		lbl_v_sat_n.setEnabled(conducting);
		lbl_v_sat_p.setEnabled(conducting);
		lbl_k_aug_n.setEnabled(semiconducting);
		lbl_k_aug_p.setEnabled(semiconducting);
		lbl_k_SRH_n.setEnabled(semiconducting);
		lbl_k_SRH_p.setEnabled(semiconducting);
	}
	
	public enum MaterialClass {
		INSULATING("Insulator"),
		CONDUCTING("Conductor"),
		SEMICONDUCTING("Semiconductor");

		String name;
		MaterialClass(String name)
		{
			this.name = name;
		}

		@Override
		public String toString() {
			return name;
		}
	}
	
	public GeneralMaterialType[] makelist() {
		ArrayList<GeneralMaterialType> matlist = new ArrayList<GeneralMaterialType>();
		for (MaterialType t : MaterialType.values()) {
			if (t != MaterialType.CUSTOM)
				matlist.add(new GeneralMaterialType(t));
		}
		
		for (int i : mat_map.keySet()) {
			matlist.add(new GeneralMaterialType(i, mat_map.get(i).name));
		}
		
		return matlist.toArray(new GeneralMaterialType[0]);
	}

	@Override
	public void itemStateChanged(ItemEvent e) {
		setInputVisibility();
	}
}
