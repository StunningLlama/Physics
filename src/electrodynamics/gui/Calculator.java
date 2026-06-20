package electrodynamics.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import electrodynamics.Material;
import electrodynamics.util.Utils;

import javax.swing.JTabbedPane;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

public class Calculator extends JDialog implements ActionListener {

	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	MaterialManager parent;
	private JTextField chi;
	private JTextField Eg;
	private JTextField gc;
	private JTextField gv;
	private JTextField chi_2;
	private JTextField Eg_2;
	private JTextField ni_2;
	private JTextField phi_3;
	private JTextField Eg_3;
	private JTextField gc_3;
	private JTextField gv_3;
	private JTextField phi_4;
	private JTextField Eg_4;
	private JTextField ni_4;
	private JButton okButton;
	private JButton cancelButton;
	private JTabbedPane tabbedPane;

	/**
	 * Create the dialog.
	 */
	public Calculator(MaterialManager parent) {
		setTitle("Band structure calculator");
		this.parent = parent;
		setBounds(100, 100, 505, 349);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new BorderLayout(0, 0));
		{
			tabbedPane = new JTabbedPane(JTabbedPane.TOP);
			contentPanel.add(tabbedPane, BorderLayout.CENTER);
			{
				JPanel panel = new JPanel();
				tabbedPane.addTab("Input χ and g", null, panel, null);
				panel.setLayout(null);
				
				JLabel lblNewLabel = new JLabel("Electron affinity [eV]");
				lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
				lblNewLabel.setBounds(136, 21, 141, 16);
				panel.add(lblNewLabel);
				
				chi = new JTextField();
				chi.setText("0");
				chi.setBounds(289, 18, 130, 23);
				panel.add(chi);
				chi.setColumns(10);
				
				JLabel lblBandGapev = new JLabel("Band gap [eV]");
				lblBandGapev.setHorizontalAlignment(SwingConstants.TRAILING);
				lblBandGapev.setBounds(136, 54, 141, 16);
				panel.add(lblBandGapev);
				
				Eg = new JTextField();
				Eg.setText("0");
				Eg.setColumns(10);
				Eg.setBounds(289, 51, 130, 23);
				panel.add(Eg);
				
				JLabel lblEffectiveCondBand = new JLabel("Effective cond. band DOS [1/m^3]");
				lblEffectiveCondBand.setHorizontalAlignment(SwingConstants.TRAILING);
				lblEffectiveCondBand.setBounds(28, 87, 249, 16);
				panel.add(lblEffectiveCondBand);
				
				gc = new JTextField();
				gc.setText("0");
				gc.setColumns(10);
				gc.setBounds(289, 84, 130, 23);
				panel.add(gc);
				
				JLabel lblEffectiveValenceBand = new JLabel("Effective valence band DOS [1/m^3]");
				lblEffectiveValenceBand.setHorizontalAlignment(SwingConstants.TRAILING);
				lblEffectiveValenceBand.setBounds(28, 120, 249, 16);
				panel.add(lblEffectiveValenceBand);
				
				gv = new JTextField();
				gv.setText("0");
				gv.setColumns(10);
				gv.setBounds(289, 117, 130, 23);
				panel.add(gv);
			}
			
			JPanel panel = new JPanel();
			panel.setLayout(null);
			tabbedPane.addTab("Input χ and ni", null, panel, null);
			
			JLabel lblNewLabel = new JLabel("Electron affinity [eV]");
			lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
			lblNewLabel.setBounds(136, 21, 141, 16);
			panel.add(lblNewLabel);
			
			chi_2 = new JTextField();
			chi_2.setText("0");
			chi_2.setColumns(10);
			chi_2.setBounds(289, 18, 130, 23);
			panel.add(chi_2);
			
			JLabel lblBandGapev = new JLabel("Band gap [eV]");
			lblBandGapev.setHorizontalAlignment(SwingConstants.TRAILING);
			lblBandGapev.setBounds(136, 54, 141, 16);
			panel.add(lblBandGapev);
			
			Eg_2 = new JTextField();
			Eg_2.setText("0");
			Eg_2.setColumns(10);
			Eg_2.setBounds(289, 51, 130, 23);
			panel.add(Eg_2);
			
			JLabel lblEffectiveCondBand = new JLabel("Equilibrium carrier density (ni) [1/m^3]");
			lblEffectiveCondBand.setHorizontalAlignment(SwingConstants.TRAILING);
			lblEffectiveCondBand.setBounds(28, 87, 249, 16);
			panel.add(lblEffectiveCondBand);
			
			ni_2 = new JTextField();
			ni_2.setText("0");
			ni_2.setColumns(10);
			ni_2.setBounds(289, 84, 130, 23);
			panel.add(ni_2);
			
			JLabel lblassumeConductionDos = new JLabel("(Assume conduction DOS = valence DOS)");
			lblassumeConductionDos.setEnabled(false);
			lblassumeConductionDos.setHorizontalAlignment(SwingConstants.TRAILING);
			lblassumeConductionDos.setBounds(85, 115, 334, 16);
			panel.add(lblassumeConductionDos);
			
			JPanel panel_1 = new JPanel();
			panel_1.setLayout(null);
			tabbedPane.addTab("Input Φ and g", null, panel_1, null);
			
			JLabel lblNewLabel_1 = new JLabel("Workfunction [eV]");
			lblNewLabel_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblNewLabel_1.setBounds(136, 21, 141, 16);
			panel_1.add(lblNewLabel_1);
			
			phi_3 = new JTextField();
			phi_3.setText("0");
			phi_3.setColumns(10);
			phi_3.setBounds(289, 18, 130, 23);
			panel_1.add(phi_3);
			
			JLabel lblBandGapev_1 = new JLabel("Band gap [eV]");
			lblBandGapev_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblBandGapev_1.setBounds(136, 54, 141, 16);
			panel_1.add(lblBandGapev_1);
			
			Eg_3 = new JTextField();
			Eg_3.setText("0");
			Eg_3.setColumns(10);
			Eg_3.setBounds(289, 51, 130, 23);
			panel_1.add(Eg_3);
			
			JLabel lblEffectiveCondBand_1 = new JLabel("Effective cond. band DOS [1/m^3]");
			lblEffectiveCondBand_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblEffectiveCondBand_1.setBounds(28, 87, 249, 16);
			panel_1.add(lblEffectiveCondBand_1);
			
			gc_3 = new JTextField();
			gc_3.setText("0");
			gc_3.setColumns(10);
			gc_3.setBounds(289, 84, 130, 23);
			panel_1.add(gc_3);
			
			JLabel lblEffectiveValenceBand_1 = new JLabel("Effective valence band DOS [1/m^3]");
			lblEffectiveValenceBand_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblEffectiveValenceBand_1.setBounds(28, 120, 249, 16);
			panel_1.add(lblEffectiveValenceBand_1);
			
			gv_3 = new JTextField();
			gv_3.setText("0");
			gv_3.setColumns(10);
			gv_3.setBounds(289, 117, 130, 23);
			panel_1.add(gv_3);
			
			JPanel panel_1_1 = new JPanel();
			panel_1_1.setLayout(null);
			tabbedPane.addTab("Input Φ and ni", null, panel_1_1, null);
			
			JLabel lblNewLabel_1_1 = new JLabel("Workfunction [eV]");
			lblNewLabel_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblNewLabel_1_1.setBounds(136, 21, 141, 16);
			panel_1_1.add(lblNewLabel_1_1);
			
			phi_4 = new JTextField();
			phi_4.setText("0");
			phi_4.setColumns(10);
			phi_4.setBounds(289, 18, 130, 23);
			panel_1_1.add(phi_4);
			
			JLabel lblBandGapev_1_1 = new JLabel("Band gap [eV]");
			lblBandGapev_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblBandGapev_1_1.setBounds(136, 54, 141, 16);
			panel_1_1.add(lblBandGapev_1_1);
			
			Eg_4 = new JTextField();
			Eg_4.setText("0");
			Eg_4.setColumns(10);
			Eg_4.setBounds(289, 51, 130, 23);
			panel_1_1.add(Eg_4);
			
			JLabel lblEffectiveCondBand_1_1 = new JLabel("Equilibrium carrier density (ni) [1/m^3]");
			lblEffectiveCondBand_1_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblEffectiveCondBand_1_1.setBounds(28, 87, 249, 16);
			panel_1_1.add(lblEffectiveCondBand_1_1);
			
			ni_4 = new JTextField();
			ni_4.setText("0");
			ni_4.setColumns(10);
			ni_4.setBounds(289, 84, 130, 23);
			panel_1_1.add(ni_4);
			
			JLabel lblassumeConductionDos_1 = new JLabel("(Assume conduction DOS = valence DOS)");
			lblassumeConductionDos_1.setHorizontalAlignment(SwingConstants.TRAILING);
			lblassumeConductionDos_1.setEnabled(false);
			lblassumeConductionDos_1.setBounds(85, 115, 334, 16);
			panel_1_1.add(lblassumeConductionDos_1);
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton = new JButton("Calculate band constants");
				okButton.setActionCommand("OK");
				okButton.addActionListener(this);
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
			{
				cancelButton = new JButton("Cancel");
				cancelButton.setActionCommand("Cancel");
				cancelButton.addActionListener(this);
				buttonPane.add(cancelButton);
			}
		}
	}

	@Override
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == okButton) {
			try {
				Material mat = new Material();
				double eVtoJ = parent.e.eVtoJ;
				double kT = parent.e.k*parent.e.T;
				switch (tabbedPane.getSelectedIndex()) {
				case 0:
					mat.computeBandstructureChiG(getDouble(chi)*eVtoJ, getDouble(Eg)*eVtoJ, getDouble(gc), getDouble(gv), kT);
					break;
				case 1:
					mat.computeBandstructureChi(getDouble(chi_2)*eVtoJ, getDouble(Eg_2)*eVtoJ, getDouble(ni_2), kT);
					break;
				case 2:
					mat.computeBandstructurePhiG(getDouble(phi_3)*eVtoJ, getDouble(Eg_3)*eVtoJ, getDouble(gc_3), getDouble(gv_3), kT);
					break;
				case 3:
					mat.computeBandstructurePhi(getDouble(phi_4)*eVtoJ, getDouble(Eg_4)*eVtoJ, getDouble(ni_4), kT);
					break;
				}

				parent.Ec.setText(Utils.formatDoubleReduced(mat.Ec/eVtoJ));
				parent.Ev.setText(Utils.formatDoubleReduced(mat.Ev/eVtoJ));
				parent.gc.setText(Utils.formatDoubleReduced(mat.gc));
				parent.gv.setText(Utils.formatDoubleReduced(mat.gv));
				this.setVisible(false);
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(this, "Invalid number.", "Error", JOptionPane.OK_OPTION);
			}
		} else if (ev.getSource() == cancelButton) {
			this.dispose();
		}
	}
	
	double getDouble(JTextField f) {
		return Double.valueOf(f.getText());
	}
}
