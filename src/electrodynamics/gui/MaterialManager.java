package electrodynamics.gui;

import java.awt.EventQueue;
import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import electrodynamics.Material;
import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.AbstractListModel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.JButton;

public class MaterialManager extends JFrame {

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
	public JTextField k_rad_semi;
	public JTextField recomb_rate_metal;
	public JTextField T;
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
	private JLabel lblDielectricConstant;
	private JTextField eps_r;
	private JLabel lblRelPermeability;
	private JTextField mu_r;
	private JLabel lblNewLabel_5_2;
	private JTextField rho_back;
	private JLabel lblNewLabel;
	private JTextField name;
	private JComboBox type;
	private JLabel lblType;
	private JButton btn_delete;
	private JButton btn_add;
	
	HashMap<Integer, Material> mat_map;
	int id_counter = 0;

	public MaterialManager() {
		setTitle("Custom materials");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 905, 350);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 6, 235, 219);
		contentPane.add(scrollPane);
		
		JList list = new JList();
		list.setModel(new AbstractListModel() {
			String[] values = new String[] {};
			public int getSize() {
				return values.length;
			}
			public Object getElementAt(int index) {
				return values[index];
			}
		});
		scrollPane.setViewportView(list);
		
		btn_add = new JButton("New material");
		btn_add.setBounds(6, 237, 117, 29);
		contentPane.add(btn_add);
		
		btn_delete = new JButton("Delete");
		btn_delete.setBounds(124, 237, 117, 29);
		contentPane.add(btn_delete);

		mu_electron = new JTextField();
		mu_electron.setColumns(10);
		mu_electron.setBounds(477, 66, 98, 26);
		contentPane.add(mu_electron);
		
		JLabel lblNewLabel_3 = new JLabel("Electron mobility [m^2/(V s)]");
		lblNewLabel_3.setToolTipText("Electron mobility");
		lblNewLabel_3.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_3.setBounds(273, 71, 192, 16);
		contentPane.add(lblNewLabel_3);
		
		JLabel lblNewLabel_4 = new JLabel("Hole mobility [m^2/(V s)]");
		lblNewLabel_4.setToolTipText("Hole mobility");
		lblNewLabel_4.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_4.setBounds(297, 104, 168, 16);
		contentPane.add(lblNewLabel_4);
		
		mu_hole = new JTextField();
		mu_hole.setColumns(10);
		mu_hole.setBounds(477, 99, 98, 26);
		contentPane.add(mu_hole);
		
		JLabel lblNewLabel_5 = new JLabel("Intrinsic carrier conc. [1/m^3]");
		lblNewLabel_5.setToolTipText("Semiconductor equilibrium carrier concentration");
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5.setBounds(267, 137, 198, 16);
		contentPane.add(lblNewLabel_5);
		
		ni_semi = new JTextField();
		ni_semi.setColumns(10);
		ni_semi.setBounds(477, 132, 98, 26);
		contentPane.add(ni_semi);
		
		JLabel lblNewLabel_6 = new JLabel("Workfunction [eV]");
		lblNewLabel_6.setToolTipText("Semiconductor work function");
		lblNewLabel_6.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_6.setBounds(297, 170, 168, 16);
		contentPane.add(lblNewLabel_6);
		
		W_semi = new JTextField();
		W_semi.setColumns(10);
		W_semi.setBounds(477, 165, 98, 26);
		contentPane.add(W_semi);
		
		JLabel lblEbsemi = new JLabel("Bandgap [eV]");
		lblEbsemi.setToolTipText("Semiconductor band gap");
		lblEbsemi.setHorizontalAlignment(SwingConstants.TRAILING);
		lblEbsemi.setBounds(297, 203, 168, 16);
		contentPane.add(lblEbsemi);
		
		E_b_semi = new JTextField();
		E_b_semi.setColumns(10);
		E_b_semi.setBounds(477, 198, 98, 26);
		contentPane.add(E_b_semi);
		
		JLabel lblNewLabel_5_1 = new JLabel("Radiative recomb. rate [m^3/s]");
		lblNewLabel_5_1.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_1.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_1.setBounds(578, 5, 211, 16);
		contentPane.add(lblNewLabel_5_1);
		
		k_rad_semi = new JTextField();
		k_rad_semi.setColumns(10);
		k_rad_semi.setBounds(801, 0, 98, 26);
		contentPane.add(k_rad_semi);
		
		lblNewLabel_5_8 = new JLabel("Electron sat. velocity [m/s]");
		lblNewLabel_5_8.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_8.setBounds(597, 170, 192, 16);
		contentPane.add(lblNewLabel_5_8);
		
		v_sat_n = new JTextField();
		v_sat_n.setColumns(10);
		v_sat_n.setBounds(801, 165, 98, 26);
		contentPane.add(v_sat_n);
		
		lblNewLabel_5_9 = new JLabel("Hole sat. velocity [m/s]");
		lblNewLabel_5_9.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_9.setBounds(597, 203, 192, 16);
		contentPane.add(lblNewLabel_5_9);
		
		v_sat_p = new JTextField();
		v_sat_p.setColumns(10);
		v_sat_p.setBounds(801, 198, 98, 26);
		contentPane.add(v_sat_p);
		
		lblNewLabel_5_10 = new JLabel("SRH recomb. rate n [1/s]");
		lblNewLabel_5_10.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_10.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_10.setBounds(597, 38, 192, 16);
		contentPane.add(lblNewLabel_5_10);
		
		k_SRH_n_semi = new JTextField();
		k_SRH_n_semi.setColumns(10);
		k_SRH_n_semi.setBounds(801, 33, 98, 26);
		contentPane.add(k_SRH_n_semi);
		
		lblNewLabel_5_11 = new JLabel("SRH recomb. rate p [1/s]");
		lblNewLabel_5_11.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_11.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_11.setBounds(607, 71, 182, 16);
		contentPane.add(lblNewLabel_5_11);
		
		k_SRH_p_semi = new JTextField();
		k_SRH_p_semi.setColumns(10);
		k_SRH_p_semi.setBounds(801, 66, 98, 26);
		contentPane.add(k_SRH_p_semi);
		
		lblNewLabel_5_12 = new JLabel("Auger recomb. rate n [m^6/s]");
		lblNewLabel_5_12.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_12.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_12.setBounds(597, 104, 192, 16);
		contentPane.add(lblNewLabel_5_12);
		
		k_aug_n_semi = new JTextField();
		k_aug_n_semi.setColumns(10);
		k_aug_n_semi.setBounds(801, 99, 98, 26);
		contentPane.add(k_aug_n_semi);
		
		lblNewLabel_5_13 = new JLabel("Auger recomb. rate p [m^6/s]");
		lblNewLabel_5_13.setToolTipText("Semiconductor carrier recombination rate (in number density units)");
		lblNewLabel_5_13.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_13.setBounds(597, 137, 192, 16);
		contentPane.add(lblNewLabel_5_13);
		
		k_aug_p_semi = new JTextField();
		k_aug_p_semi.setColumns(10);
		k_aug_p_semi.setBounds(801, 132, 98, 26);
		contentPane.add(k_aug_p_semi);
		
		lblDielectricConstant = new JLabel("Dielectric constant");
		lblDielectricConstant.setToolTipText("Semiconductor band gap");
		lblDielectricConstant.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDielectricConstant.setBounds(297, 235, 168, 16);
		contentPane.add(lblDielectricConstant);
		
		eps_r = new JTextField();
		eps_r.setColumns(10);
		eps_r.setBounds(477, 230, 98, 26);
		contentPane.add(eps_r);
		
		lblRelPermeability = new JLabel("Rel. permeability");
		lblRelPermeability.setToolTipText("Semiconductor band gap");
		lblRelPermeability.setHorizontalAlignment(SwingConstants.TRAILING);
		lblRelPermeability.setBounds(297, 268, 168, 16);
		contentPane.add(lblRelPermeability);
		
		mu_r = new JTextField();
		mu_r.setColumns(10);
		mu_r.setBounds(477, 263, 98, 26);
		contentPane.add(mu_r);
		
		lblNewLabel_5_2 = new JLabel("Dopant charge density [C/m^3]");
		lblNewLabel_5_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel_5_2.setBounds(578, 236, 211, 16);
		contentPane.add(lblNewLabel_5_2);
		
		rho_back = new JTextField();
		rho_back.setColumns(10);
		rho_back.setBounds(801, 231, 98, 26);
		contentPane.add(rho_back);
		
		lblNewLabel = new JLabel("Name");
		lblNewLabel.setToolTipText("Electron mobility");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(273, 6, 88, 16);
		contentPane.add(lblNewLabel);
		
		name = new JTextField();
		name.setColumns(10);
		name.setBounds(373, 0, 202, 26);
		contentPane.add(name);
		
		type = new JComboBox();
		type.setBounds(373, 34, 202, 27);
		contentPane.add(type);
		
		lblType = new JLabel("Type");
		lblType.setToolTipText("Electron mobility");
		lblType.setHorizontalAlignment(SwingConstants.TRAILING);
		lblType.setBounds(273, 38, 88, 16);
		contentPane.add(lblType);
	}
}
