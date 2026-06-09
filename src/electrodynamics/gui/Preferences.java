// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.stream.JsonReader;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;

import electrodynamics.Simulation;
import electrodynamics.units.Units;

public class Preferences extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
	public JButton btn_apply;
	public JButton btn_cancel;
	public JCheckBox chkbox_undo;
	public JSpinner spinner_imgx;
	public JCheckBox chkbox_potential;
	public JComboBox<Units> gui_units;
	
	Simulation e;
	public int saveversion = 1;
	public File preferences_file = null;
	public JCheckBox chkbox_logscale;
	private JSpinner spinner_imgy;

	public Preferences(Simulation e) {
		setResizable(false);
		this.e = e;
		preferences_file = new File("preferences.json");
		
		setTitle("Preferences");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 315, 281);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JLabel lblNewLabel = new JLabel("Display width [px]");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setToolTipText("Width of simulation domain");
		lblNewLabel.setBounds(17, 18, 138, 16);
		contentPane.add(lblNewLabel);
		
		btn_apply = new JButton("Save");
		btn_apply.setBounds(46, 218, 128, 29);
		contentPane.add(btn_apply);
		
		btn_cancel = new JButton("Reset to defaults");
		btn_cancel.setBounds(173, 218, 136, 29);
		contentPane.add(btn_cancel);
		
		spinner_imgx = new JSpinner();
		spinner_imgx.setBounds(167, 13, 109, 26);
		contentPane.add(spinner_imgx);
		
		chkbox_undo = new JCheckBox("Undo tracks settings");
		chkbox_undo.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_undo.setSelected(true);
		chkbox_undo.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_undo.setBounds(52, 74, 224, 23);
		contentPane.add(chkbox_undo);
		
		chkbox_potential = new JCheckBox("Display potential relative to ground");
		chkbox_potential.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_potential.setSelected(true);
		chkbox_potential.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_potential.setBounds(21, 100, 255, 23);
		contentPane.add(chkbox_potential);
		
		JLabel lblUnitSystem = new JLabel("Unit system");
		lblUnitSystem.setHorizontalAlignment(SwingConstants.TRAILING);
		lblUnitSystem.setBounds(6, 165, 116, 16);
		contentPane.add(lblUnitSystem);
		
		gui_units = new JComboBox<>();
		gui_units.setModel(new DefaultComboBoxModel<>(Units.values()));
		gui_units.setBounds(134, 161, 146, 27);
		contentPane.add(gui_units);
		
		chkbox_logscale = new JCheckBox("Use log scale for density plots");
		chkbox_logscale.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_logscale.setSelected(true);
		chkbox_logscale.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_logscale.setBounds(21, 126, 255, 23);
		contentPane.add(chkbox_logscale);
		
		JLabel lblDisplayHeightpx = new JLabel("Display height [px]");
		lblDisplayHeightpx.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDisplayHeightpx.setToolTipText("Width of simulation domain");
		lblDisplayHeightpx.setBounds(17, 46, 138, 16);
		contentPane.add(lblDisplayHeightpx);
		
		spinner_imgy = new JSpinner();
		spinner_imgy.setBounds(167, 41, 109, 26);
		contentPane.add(spinner_imgy);
		
		resetPrefs();
	}
	
	public void initialize() {
		btn_apply.addActionListener(this);
		btn_cancel.addActionListener(this);
		gui_units.addActionListener(this);
		setLocationRelativeTo(null);
		setVisible(false);
		
		readfile(preferences_file);
		applyPrefs();
	}
	
	public void getPrefs() {
		spinner_imgx.setValue(e.canvas.getWidth());
		spinner_imgy.setValue(e.canvas.getHeight());
	}

	public void applyPrefs() {
		int x = (int)(spinner_imgx.getValue());
		int y = (int)(spinner_imgy.getValue());

		e.canvas.setPreferredSize(new Dimension(x, y));
		e.opts.pack();
	}
	
	public void resetPrefs() {
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int opt_height = 256*(int)Math.floor(0.8*screenSize.getHeight()/256);
		spinner_imgx.setValue(opt_height);
		spinner_imgy.setValue(opt_height);
		chkbox_undo.setSelected(true);
		chkbox_potential.setSelected(true);
		chkbox_logscale.setSelected(true);
		gui_units.setSelectedItem(Units.SI);
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

				if (version > saveversion) {
					fstr.close();
					return;
				}

				if (version == saveversion) {
					assertNextObject(fstr, "preferences");
					fstr.beginObject();
					while (fstr.hasNext()) {
						String name = fstr.nextName();
						switch (name){
						case "imgsize":
							int size = fstr.nextInt();
							spinner_imgx.setValue(size);
							spinner_imgy.setValue(size);
							break;
						case "imgsize_x": spinner_imgx.setValue(fstr.nextInt()); break;
						case "imgsize_y": spinner_imgy.setValue(fstr.nextInt()); break;
						case "undotrackssettings": chkbox_undo.setSelected(fstr.nextBoolean()); break;
						case "potential": chkbox_potential.setSelected(fstr.nextBoolean()); break;
						case "logscale": chkbox_logscale.setSelected(fstr.nextBoolean()); break;
						case "units": gui_units.setSelectedItem(gson.fromJson(fstr, Units.class)); break;
						}
					}
					fstr.endObject();
					fstr.close();
				}
			} catch (FileNotFoundException ex) {
				return;
			} catch (IOException | IllegalArgumentException ex) {
				return;
			}
			return;
		} finally {
			e.rwLock.writeLock().unlock();
		}
	}

	public void writeFile(File outfile)
	{
		e.rwLock.writeLock().lock();
		try {
			try {
				PrintWriter fstr = new PrintWriter(new FileOutputStream(outfile));

				Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

				JsonObject header = new JsonObject();
				header.addProperty("imgsize_x", (int)spinner_imgx.getValue());
				header.addProperty("imgsize_y", (int)spinner_imgy.getValue());
				header.addProperty("undotrackssettings", chkbox_undo.isSelected());
				header.addProperty("potential", chkbox_potential.isSelected());
				header.addProperty("logscale", chkbox_logscale.isSelected());
				header.add("units", gson.toJsonTree((Units) gui_units.getSelectedItem()));


				// Version should always be first
				JsonObject save = new JsonObject();
				save.addProperty("version", saveversion);
				save.add("preferences", header);

				String json = gson.toJson(save);

				fstr.print(json);
				fstr.flush();
				fstr.close();
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

	@Override
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == btn_apply) {
			applyPrefs();
			writeFile(preferences_file);
			setVisible(false);
		} else if (ev.getSource() == btn_cancel) {
			resetPrefs();
		}  else if (ev.getSource() == gui_units) {
			e.units = (Units) gui_units.getSelectedItem();
		} 
	}
}
