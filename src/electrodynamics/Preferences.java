// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

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
import electrodynamics.units.Units;

public class Preferences extends JFrame {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
	public JButton btn_apply;
	public JButton btn_cancel;
	public JCheckBox chkbox_undo;
	public JSpinner spinner_imgsize;
	public JCheckBox chkbox_potential;
	public JComboBox gui_units;
	
	Simulation e;
	int saveversion = 1;
	File preferences_file = null;
	public JCheckBox chkbox_logscale;

	public Preferences(Simulation e) {
		this.e = e;
		preferences_file = new File("preferences.json");
		
		setTitle("Preferences");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 313, 237);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JLabel lblNewLabel = new JLabel("Display size [px]");
		lblNewLabel.setHorizontalAlignment(SwingConstants.CENTER);
		lblNewLabel.setToolTipText("Width of simulation domain");
		lblNewLabel.setBounds(17, 18, 125, 16);
		contentPane.add(lblNewLabel);
		
		btn_apply = new JButton("Save");
		btn_apply.setBounds(72, 174, 100, 29);
		contentPane.add(btn_apply);
		
		btn_cancel = new JButton("Reset to defaults");
		btn_cancel.setBounds(171, 174, 136, 29);
		contentPane.add(btn_cancel);
		
		spinner_imgsize = new JSpinner();
		spinner_imgsize.setBounds(167, 13, 109, 26);
		contentPane.add(spinner_imgsize);
		
		chkbox_undo = new JCheckBox("Undo tracks settings");
		chkbox_undo.setSelected(true);
		chkbox_undo.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_undo.setBounds(21, 42, 224, 23);
		contentPane.add(chkbox_undo);
		
		chkbox_potential = new JCheckBox("Display potential relative to ground");
		chkbox_potential.setSelected(true);
		chkbox_potential.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_potential.setBounds(21, 69, 255, 23);
		contentPane.add(chkbox_potential);
		
		JLabel lblUnitSystem = new JLabel("Unit system");
		lblUnitSystem.setToolTipText("Width of simulation domain");
		lblUnitSystem.setHorizontalAlignment(SwingConstants.CENTER);
		lblUnitSystem.setBounds(17, 128, 91, 16);
		contentPane.add(lblUnitSystem);
		
		gui_units = new JComboBox();
		gui_units.setModel(new DefaultComboBoxModel(Units.values()));
		gui_units.setBounds(130, 124, 145, 27);
		contentPane.add(gui_units);
		
		chkbox_logscale = new JCheckBox("Use log scale for density plots");
		chkbox_logscale.setSelected(true);
		chkbox_logscale.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_logscale.setBounds(21, 95, 255, 23);
		contentPane.add(chkbox_logscale);
	}
	
	public void initialize() {
		btn_apply.addActionListener(e.controls);
		btn_cancel.addActionListener(e.controls);
		gui_units.addActionListener(e.controls);
		setVisible(false);
		
		readfile(preferences_file);
		applyPrefs();
	}
	
	public void getPrefs() {
		spinner_imgsize.setValue(e.renderer.canvas_size);
	}

	public void applyPrefs() {
		int x = (int)(spinner_imgsize.getValue());

		if (e.renderer.canvas_size != x) {
			e.renderer.new_canvas_size = x;
			e.controls.updateimagesize = true;
		}
	}
	
	public void resetPrefs() {
		spinner_imgsize.setValue(768);
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
						case "imgsize": spinner_imgsize.setValue(fstr.nextInt()); break;
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
				header.addProperty("imgsize", (int)spinner_imgsize.getValue());
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
}
