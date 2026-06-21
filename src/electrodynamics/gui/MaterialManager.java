package electrodynamics.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileFilter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.stream.JsonReader;

import electrodynamics.GeneralMaterialType;
import electrodynamics.Material;
import electrodynamics.MaterialType;
import electrodynamics.SemiSim;
import electrodynamics.Simulation;
import electrodynamics.Steam;
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
	public JTextField gc;
	public JTextField Ec;
	public JTextField Ev;
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
	private JButton btn_import;
	private JButton btn_export;

	public MaterialManager(Simulation e) {
		this.e = e;
		startingpath = SemiSim.userdir;
		setResizable(false);
		setTitle("Material editor");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 920, 422);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 6, 235, 341);
		contentPane.add(scrollPane);
		
		list = new JList<>();
		scrollPane.setViewportView(list);
		
		btn_add = new JButton("New material");
		btn_add.setBounds(6, 359, 117, 23);
		contentPane.add(btn_add);
		
		btn_delete = new JButton("Delete");
		btn_delete.setBounds(124, 359, 117, 23);
		contentPane.add(btn_delete);

		mu_electron = new JTextField();
		mu_electron.setColumns(10);
		mu_electron.setBounds(476, 265, 98, 23);
		contentPane.add(mu_electron);
		
		lbl_mu_electron = new JLabel("Electron mobility [m^2/(V s)]");
		lbl_mu_electron.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_electron.setBounds(273, 269, 192, 16);
		contentPane.add(lbl_mu_electron);
		
		lbl_mu_hole = new JLabel("Hole mobility [m^2/(V s)]");
		lbl_mu_hole.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_hole.setBounds(297, 302, 168, 16);
		contentPane.add(lbl_mu_hole);
		
		mu_hole = new JTextField();
		mu_hole.setColumns(10);
		mu_hole.setBounds(476, 298, 98, 23);
		contentPane.add(mu_hole);
		
		lbl_ni = new JLabel("Cond. band eff. DOS [1/m^3]");
		lbl_ni.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_ni.setBounds(267, 203, 198, 16);
		contentPane.add(lbl_ni);
		
		gc = new JTextField();
		gc.setColumns(10);
		gc.setBounds(476, 199, 98, 23);
		contentPane.add(gc);
		
		lbl_W = new JLabel("Conduction band energy [eV]");
		lbl_W.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_W.setBounds(253, 137, 212, 16);
		contentPane.add(lbl_W);
		
		Ec = new JTextField();
		Ec.setColumns(10);
		Ec.setBounds(476, 133, 98, 23);
		contentPane.add(Ec);
		
		lbl_Eb = new JLabel("Valence band energy [eV]");
		lbl_Eb.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_Eb.setBounds(283, 170, 182, 16);
		contentPane.add(lbl_Eb);
		
		Ev = new JTextField();
		Ev.setColumns(10);
		Ev.setBounds(476, 166, 98, 23);
		contentPane.add(Ev);
		
		lbl_k_rad = new JLabel("Radiative recomb. rate [m^3/s]");
		lbl_k_rad.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_rad.setBounds(578, 170, 211, 16);
		contentPane.add(lbl_k_rad);
		
		k_rad = new JTextField();
		k_rad.setColumns(10);
		k_rad.setBounds(800, 166, 98, 23);
		contentPane.add(k_rad);
		
		lbl_v_sat_n = new JLabel("Electron sat. velocity [m/s]");
		lbl_v_sat_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_v_sat_n.setBounds(597, 104, 192, 16);
		contentPane.add(lbl_v_sat_n);
		
		v_sat_n = new JTextField();
		v_sat_n.setColumns(10);
		v_sat_n.setBounds(800, 100, 98, 23);
		contentPane.add(v_sat_n);
		
		lbl_v_sat_p = new JLabel("Hole sat. velocity [m/s]");
		lbl_v_sat_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_v_sat_p.setBounds(597, 137, 192, 16);
		contentPane.add(lbl_v_sat_p);
		
		v_sat_p = new JTextField();
		v_sat_p.setColumns(10);
		v_sat_p.setBounds(800, 133, 98, 23);
		contentPane.add(v_sat_p);
		
		lbl_k_SRH_n = new JLabel("SRH recomb. rate n [1/s]");
		lbl_k_SRH_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_SRH_n.setBounds(597, 203, 192, 16);
		contentPane.add(lbl_k_SRH_n);
		
		k_SRH_n = new JTextField();
		k_SRH_n.setColumns(10);
		k_SRH_n.setBounds(800, 199, 98, 23);
		contentPane.add(k_SRH_n);
		
		lbl_k_SRH_p = new JLabel("SRH recomb. rate p [1/s]");
		lbl_k_SRH_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_SRH_p.setBounds(607, 236, 182, 16);
		contentPane.add(lbl_k_SRH_p);
		
		k_SRH_p = new JTextField();
		k_SRH_p.setColumns(10);
		k_SRH_p.setBounds(800, 232, 98, 23);
		contentPane.add(k_SRH_p);
		
		lbl_k_aug_n = new JLabel("Auger recomb. rate n [m^6/s]");
		lbl_k_aug_n.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_aug_n.setBounds(597, 269, 192, 16);
		contentPane.add(lbl_k_aug_n);
		
		k_aug_n = new JTextField();
		k_aug_n.setColumns(10);
		k_aug_n.setBounds(800, 265, 98, 23);
		contentPane.add(k_aug_n);
		
		lbl_k_aug_p = new JLabel("Auger recomb. rate p [m^6/s]");
		lbl_k_aug_p.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_k_aug_p.setBounds(597, 302, 192, 16);
		contentPane.add(lbl_k_aug_p);
		
		k_aug_p = new JTextField();
		k_aug_p.setColumns(10);
		k_aug_p.setBounds(800, 298, 98, 23);
		contentPane.add(k_aug_p);
		
		lbl_eps_r = new JLabel("Dielectric constant");
		lbl_eps_r.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_eps_r.setBounds(297, 71, 168, 16);
		contentPane.add(lbl_eps_r);
		
		eps_r = new JTextField();
		eps_r.setColumns(10);
		eps_r.setBounds(476, 67, 98, 23);
		contentPane.add(eps_r);
		
		lbl_mu_r = new JLabel("Relative permeability");
		lbl_mu_r.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_mu_r.setBounds(297, 104, 168, 16);
		contentPane.add(lbl_mu_r);
		
		mu_r = new JTextField();
		mu_r.setColumns(10);
		mu_r.setBounds(476, 100, 98, 23);
		contentPane.add(mu_r);
		
		lbl_rho_back = new JLabel("Dopant charge density [C/m^3]");
		lbl_rho_back.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_rho_back.setBounds(578, 71, 211, 16);
		contentPane.add(lbl_rho_back);
		
		rho_back = new JTextField();
		rho_back.setColumns(10);
		rho_back.setBounds(800, 67, 98, 23);
		contentPane.add(rho_back);
		
		lblNewLabel = new JLabel("Name");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(273, 6, 88, 16);
		contentPane.add(lblNewLabel);
		
		name = new JTextField();
		name.setColumns(10);
		name.setBounds(372, 3, 202, 23);
		contentPane.add(name);
		
		type = new JComboBox<>();
		type.setModel(new DefaultComboBoxModel<>(MaterialClass.values()));
		type.setBounds(372, 35, 202, 23);
		contentPane.add(type);
		
		lblType = new JLabel("Type");
		lblType.setHorizontalAlignment(SwingConstants.TRAILING);
		lblType.setBounds(273, 38, 88, 16);
		contentPane.add(lblType);
		
		btn_apply = new JButton("Apply");
		btn_apply.setBounds(651, 359, 124, 23);
		contentPane.add(btn_apply);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(775, 359, 124, 23);
		contentPane.add(btn_cancel);
		
		btn_import = new JButton("Import");
		btn_import.setBounds(244, 359, 117, 23);
		contentPane.add(btn_import);
		
		btn_export = new JButton("Export");
		btn_export.setBounds(363, 359, 124, 23);
		contentPane.add(btn_export);
		
		lbl_gv = new JLabel("Valence band eff. DOS [1/m^3]");
		lbl_gv.setHorizontalAlignment(SwingConstants.TRAILING);
		lbl_gv.setBounds(267, 236, 198, 16);
		contentPane.add(lbl_gv);
		
		gv = new JTextField();
		gv.setColumns(10);
		gv.setBounds(476, 232, 98, 23);
		contentPane.add(gv);
		
		btn_calc = new JButton("Band and DOS calculator");
		btn_calc.setBounds(650, 34, 212, 23);
		contentPane.add(btn_calc);
	}
	
	public void initialize() {
		this.btn_add.addActionListener(this);
		this.btn_delete.addActionListener(this);
		this.btn_apply.addActionListener(this);
		this.btn_cancel.addActionListener(this);
		this.btn_import.addActionListener(this);
		this.btn_export.addActionListener(this);
		this.btn_calc.addActionListener(this);
		for (Component c : contentPane.getComponents()) {
			if (c instanceof JTextField)
				((JTextField) c).addActionListener(this);
		}
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
		Material mat = null;
		
		GeneralMaterialType[] arr = makelist();
		
		int ind = Arrays.asList(arr).indexOf(new GeneralMaterialType(MaterialType.VACUUM));

        JList<GeneralMaterialType> tmplist = new JList<>(arr);

        tmplist.setVisibleRowCount(5);

        JScrollPane scrollPane = new JScrollPane(tmplist);
        scrollPane.setPreferredSize(new Dimension(400, 400));

        tmplist.setSelectedValue(arr[ind], true);

        int result = JOptionPane.showConfirmDialog(
                null,
                scrollPane,
                "Select template",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE
        );
        
        String new_name = "New material";

        if (result == JOptionPane.OK_OPTION) {
        	for (GeneralMaterialType selected: tmplist.getSelectedValuesList()) {
        		mat = new Material();
        		
            	if (selected.cust_id == -1)
            		e.initializeMaterial(mat, selected.type);
            	else
                    mat.copyFrom(mat_map.get(selected.cust_id));
            	
            	new_name = selected.toString() + " copy";
            	
                mat.type = MaterialType.CUSTOM;
        		mat.name = new_name;
        		mat.cust_id = id_counter;
        		mat_map.put(mat.cust_id, mat);
        		id_counter++;
        	}
        } else {
        	return;
        }
		
		updateUI();
		e.initializeAllMaterials();
		e.updateAllMaterials(true);
		list.setSelectedValue(mat, true);

		Steam.setAchievement("NEW_MATERIAL");
	}
	
	public void delete() {
		List<Material> mats = list.getSelectedValuesList();
		for (Material mat : mats) {
			mat_map.remove(mat.cust_id);
		}
		updateUI();
		e.initializeAllMaterials();
		e.updateAllMaterials(true);
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
		}else if (e.getSource() == btn_import) {
			this.readFile();
		} else if (e.getSource() == btn_export) {
			this.writeFile();
		} else if (e.getSource() == btn_calc) {
			Calculator calc = new Calculator(this);
			calc.setVisible(true);
		} else if (e.getSource() instanceof JTextField) {
			save();
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
		gc.setText(Utils.formatDouble(mat.gc));
		gv.setText(Utils.formatDouble(mat.gv));
		Ec.setText(Utils.formatDouble(mat.Ec/e.eVtoJ));
		Ev.setText(Utils.formatDouble(mat.Ev/e.eVtoJ));
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
			mat.gc = Double.valueOf(gc.getText());
			mat.gv = Double.valueOf(gv.getText());
			mat.Ec = Double.valueOf(Ec.getText())*e.eVtoJ;
			mat.Ev = Double.valueOf(Ev.getText())*e.eVtoJ;
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
		
		gc.setVisible(conducting);
		Ec.setVisible(conducting);
		gv.setVisible(conducting);
		Ev.setVisible(conducting);
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
		lbl_gv.setEnabled(conducting);
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

	public int current_material_saveversion = 1;
	public String fileextension = ".material";
	public Path startingpath;
	JTextField gv;
	private JLabel lbl_gv;
	private JButton btn_calc;
	
	public void readFile()
	{
		SwingUtilities.invokeLater(() -> {
			File testfile = startingpath.toFile();
			if (!testfile.canRead()) {
				JOptionPane.showMessageDialog(this,
				"Error: Java does not have access to this folder. Please see instructions to fix this issue.");
			}

			JFileChooser fd = new JFileChooser(startingpath.toFile());
			fd.setDialogTitle("Import material(s)");
			fd.setMultiSelectionEnabled(true);
			fd.setFileFilter(new FileFilter(){
				@Override
				public boolean accept(File f) {
					if (f.isDirectory() || f.getName().endsWith(fileextension)) return true;
					return false;
				}
				@Override
				public String getDescription() {
					return fileextension;
				}
			});
			fd.setVisible(true);
			int result = fd.showOpenDialog(this);
			startingpath = fd.getCurrentDirectory().toPath();

			if (result == JFileChooser.APPROVE_OPTION) {
				File[] files = fd.getSelectedFiles();
				for (File infile : files)
					readfile(infile);
			}
		});
	}

	public void readfile(File infile) {
		e.rwLock.writeLock().lock();
		try {
			if (infile == null || !infile.exists()) return;
			
			try {
				JsonReader fstr = new JsonReader(new InputStreamReader(new FileInputStream(infile)));
				Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

				fstr.beginObject();

				assertNextObject(fstr, "version");
				int version = fstr.nextInt();
				
				System.out.println("Loading " + infile.getName() + ", version = " + version);

				if (version > current_material_saveversion) {
					fstr.close();
					throw new IllegalArgumentException("The file was created in a newer version of SemiSim.");
				}

				if (version == current_material_saveversion) {
					assertNextObject(fstr, "materials");
					Material[] mats = (Material[]) gson.fromJson(fstr, Material[].class);

					for (Material mat : mats) {
						mat.cust_id = id_counter;
						mat_map.put(mat.cust_id, mat);
						id_counter++;
						updateUI();
						list.setSelectedValue(mat, true);
					}

					fstr.close();

				}
			} catch (FileNotFoundException ex) {
				return;
			} catch (IOException | IllegalArgumentException ex) {
				JOptionPane.showMessageDialog(this,
				"Unable to load file.\n" + ex.getMessage());
				ex.printStackTrace();
				return;
			}
			return;
		} finally {
			e.rwLock.writeLock().unlock();
		}
	}
	
	public void writeFile()
	{
		SwingUtilities.invokeLater(() -> {
			if (list.getSelectedValue() == null) {
				JOptionPane.showMessageDialog(this,
				"Please select a material to save.");
				return;
			}
			
			File testfile = startingpath.toFile();
			if (!testfile.canWrite()) {
				JOptionPane.showMessageDialog(this,
				"Error: Java does not have access to this folder. Please see instructions to fix this issue.");
				return;
			}
			
			JFileChooser fd = new JFileChooser(startingpath.toFile());
			fd.setSelectedFile(new File(list.getSelectedValue().toString() + ".material"));
			fd.setDialogTitle("Export material(s)");
			fd.setFileFilter(new FileFilter(){
				@Override
				public boolean accept(File f) {
					if (f.isDirectory() || f.getName().endsWith(fileextension)) return true;
					return false;
				}
				@Override
				public String getDescription() {
					return fileextension;
				}
			});
			int result = fd.showSaveDialog(this);
			startingpath = fd.getCurrentDirectory().toPath();

			File outfile = null;
			
			if (result == JFileChooser.APPROVE_OPTION)
				outfile = fd.getSelectedFile();

			if (outfile == null) return;
			if (!outfile.getName().endsWith(fileextension))
				outfile = new File(outfile.getAbsolutePath() + fileextension);

			if (outfile.exists()) {
				String[] options = {"Yes", "No"};

				result = JOptionPane.showOptionDialog(this, "A file named " + outfile.getName() + " already exists. Do you wish to overwrite it?", "Message", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
				if (result != JOptionPane.OK_OPTION)
					return;
			}
			
			writeFile(outfile);
		});
	}

	public void writeFile(File outfile)
	{
		e.rwLock.writeLock().lock();
		try {
			try {
				PrintWriter fstr = new PrintWriter(new FileOutputStream(outfile));

				Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().setPrettyPrinting().create();

				JsonObject save = new JsonObject();
				save.addProperty("version", current_material_saveversion);
				Material[] mats = list.getSelectedValuesList().toArray(new Material[0]);
				save.add("materials", gson.toJsonTree(mats));
				String json = gson.toJson(save);

				fstr.print(json);
				fstr.flush();
				fstr.close();

				JOptionPane.showMessageDialog(this, "Material(s) saved successfully.");
			} catch (FileNotFoundException e) {
				return;
			}
			return;
		} finally {
			e.rwLock.writeLock().unlock();
		}
	}
	

	public void assertNextObject(JsonReader fstr, String name) throws IOException {
		if (!fstr.nextName().equals(name)) {
			fstr.close();
			throw new IllegalArgumentException(name + " not found in file.");
		}
	}
	
	public boolean testNextObject(JsonReader fstr, String name) throws IOException {
		return fstr.nextName().equals(name);
	}
}
