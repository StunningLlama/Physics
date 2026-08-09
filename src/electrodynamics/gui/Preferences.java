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
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.border.EmptyBorder;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.stream.JsonReader;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;

import electrodynamics.Renderer.Text;
import electrodynamics.SemiSim;
import electrodynamics.Simulation;
import electrodynamics.units.Units;
import javax.swing.SpinnerNumberModel;

public class Preferences extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
	private JButton btn_apply;
	private JButton btn_reset;
	private JCheckBox chkbox_undo;
	private JSpinner spinner_imgx;
	private JCheckBox chkbox_potential;
	public JComboBox<Units> gui_units;
	
	Simulation e;
	public int current_pref_saveversion = 1;
	public File preferences_file = null;
	private JSpinner spinner_imgy;
	private JLabel lblUndoHistorySize_2;
	private JSpinner spinner_undosize;
	private JSpinner spinner_fps;
	private JLabel lblFontSize;
	private JSpinner spinner_font;
	private JButton btn_cancel;
	private JCheckBox chkbox_matname;
	private JComboBox<Theme> gui_lookfeel;
	private JCheckBox chkbox_voltage;
	private JSpinner spinner_fps_1;

	public Preferences(Simulation e) {
		setResizable(false);
		this.e = e;
		preferences_file = SemiSim.getUserFile("preferences.json");
		
		setTitle("Preferences");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 581, 307);
		contentPane = new JPanel();
		contentPane.setPreferredSize(new Dimension(581, 250));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JLabel lblNewLabel = new JLabel("Display width [px]");
		lblNewLabel.setHorizontalAlignment(SwingConstants.TRAILING);
		lblNewLabel.setBounds(17, 18, 138, 16);
		contentPane.add(lblNewLabel);
		
		btn_apply = new JButton("Apply");
		btn_apply.setBounds(317, 215, 128, 23);
		contentPane.add(btn_apply);
		
		btn_reset = new JButton("Reset to defaults");
		btn_reset.setBounds(10, 215, 136, 23);
		contentPane.add(btn_reset);
		
		spinner_imgx = new JSpinner();
		spinner_imgx.setModel(new SpinnerNumberModel(768, 1, 10000, 1));
		spinner_imgx.setBounds(167, 15, 109, 23);
		contentPane.add(spinner_imgx);
		
		chkbox_undo = new JCheckBox("Undo tracks settings");
		chkbox_undo.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_undo.setSelected(true);
		chkbox_undo.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_undo.setBounds(338, 13, 224, 23);
		contentPane.add(chkbox_undo);
		
		chkbox_potential = new JCheckBox("Display potential relative to ground");
		chkbox_potential.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_potential.setSelected(true);
		chkbox_potential.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_potential.setBounds(307, 42, 255, 23);
		contentPane.add(chkbox_potential);
		
		JLabel lblUnitSystem = new JLabel("Unit system");
		lblUnitSystem.setHorizontalAlignment(SwingConstants.TRAILING);
		lblUnitSystem.setBounds(288, 134, 116, 16);
		contentPane.add(lblUnitSystem);
		
		gui_units = new JComboBox<>();
		gui_units.setModel(new DefaultComboBoxModel<>(Units.values()));
		gui_units.setBounds(416, 132, 146, 23);
		contentPane.add(gui_units);
		
		JLabel lblDisplayHeightpx = new JLabel("Display height [px]");
		lblDisplayHeightpx.setHorizontalAlignment(SwingConstants.TRAILING);
		lblDisplayHeightpx.setBounds(17, 46, 138, 16);
		contentPane.add(lblDisplayHeightpx);
		
		spinner_imgy = new JSpinner();
		spinner_imgy.setModel(new SpinnerNumberModel(768, 1, 10000, 1));
		spinner_imgy.setBounds(167, 43, 109, 23);
		contentPane.add(spinner_imgy);
		
		JLabel lblUndoHistorySize = new JLabel("Target graphics FPS");
		lblUndoHistorySize.setHorizontalAlignment(SwingConstants.TRAILING);
		lblUndoHistorySize.setBounds(17, 104, 138, 16);
		contentPane.add(lblUndoHistorySize);
		
		spinner_fps = new JSpinner();
		spinner_fps.setModel(new SpinnerNumberModel(60, 1, 1000, 1));
		spinner_fps.setBounds(167, 101, 109, 23);
		contentPane.add(spinner_fps);
		
		lblUndoHistorySize_2 = new JLabel("Undo history size");
		lblUndoHistorySize_2.setToolTipText("Warning: uses memory");
		lblUndoHistorySize_2.setHorizontalAlignment(SwingConstants.TRAILING);
		lblUndoHistorySize_2.setBounds(17, 162, 138, 16);
		contentPane.add(lblUndoHistorySize_2);
		
		spinner_undosize = new JSpinner();
		spinner_undosize.setModel(new SpinnerNumberModel(4, 2, 100, 1));
		spinner_undosize.setBounds(167, 159, 109, 23);
		contentPane.add(spinner_undosize);
		
		lblFontSize = new JLabel("Font size [px]");
		lblFontSize.setHorizontalAlignment(SwingConstants.TRAILING);
		lblFontSize.setBounds(17, 75, 138, 16);
		contentPane.add(lblFontSize);
		
		spinner_font = new JSpinner();
		spinner_font.setModel(new SpinnerNumberModel(12, 1, 100, 1));
		spinner_font.setBounds(167, 72, 109, 23);
		contentPane.add(spinner_font);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(445, 215, 128, 23);
		contentPane.add(btn_cancel);
		
		chkbox_matname = new JCheckBox("Show material name next to cursor");
		chkbox_matname.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_matname.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_matname.setBounds(307, 71, 255, 23);
		contentPane.add(chkbox_matname);
		
		gui_lookfeel = new JComboBox<>();
		gui_lookfeel.setBounds(416, 161, 146, 23);
		contentPane.add(gui_lookfeel);
		
		JLabel lblUiTheme = new JLabel("UI theme");
		lblUiTheme.setHorizontalAlignment(SwingConstants.TRAILING);
		lblUiTheme.setBounds(288, 164, 116, 16);
		contentPane.add(lblUiTheme);

		gui_lookfeel.setModel(new DefaultComboBoxModel<Theme>(Theme.values));
		
		chkbox_voltage = new JCheckBox("Set voltage instead of EMF");
		chkbox_voltage.setSelected(true);
		chkbox_voltage.setActionCommand("");
		chkbox_voltage.setHorizontalTextPosition(SwingConstants.LEADING);
		chkbox_voltage.setHorizontalAlignment(SwingConstants.TRAILING);
		chkbox_voltage.setBounds(307, 99, 255, 23);
		contentPane.add(chkbox_voltage);
		
		JLabel lblTargetSimulationFps = new JLabel("Target simulation FPS");
		lblTargetSimulationFps.setHorizontalAlignment(SwingConstants.TRAILING);
		lblTargetSimulationFps.setBounds(17, 133, 138, 16);
		contentPane.add(lblTargetSimulationFps);
		
		spinner_fps_1 = new JSpinner();
		spinner_fps_1.setModel(new SpinnerNumberModel(60, 1, 1000, 1));
		spinner_fps_1.setBounds(167, 130, 109, 23);
		contentPane.add(spinner_fps_1);
		
		pack();
		
		resetPrefs();
	}
	
	public void initialize() {
		btn_apply.addActionListener(this);
		btn_reset.addActionListener(this);
		btn_cancel.addActionListener(this);
		
		setLocationRelativeTo(null);
		setVisible(false);
		
		readfile(preferences_file);
		
		if (!preferences_file.exists()) {
			writeFile(preferences_file);
		}
		
		applyPrefs();
	}
	
	public void getPrefs() {
		spinner_imgx.setValue(e.canvas.getWidth());
		spinner_imgy.setValue(e.canvas.getHeight());
		chkbox_undo.setSelected(e.controls.undoredo.tracksettings);
		spinner_undosize.setValue(e.controls.undoredo.history_size);
		spinner_fps.setValue((int) e.renderer.targetframerate);
		spinner_fps_1.setValue((int) e.targetframerate);
		chkbox_potential.setSelected(e.renderer.display_relative_voltage);
		spinner_font.setValue(Text.fontsize);
		gui_units.setSelectedItem(e.units);
		chkbox_matname.setSelected(e.renderer.disp_mat_name);
		chkbox_voltage.setSelected(e.controls.setvoltage);
		gui_lookfeel.setSelectedItem(new Theme(UIManager.getLookAndFeel().getClass().getName()));
	}

	public void applyPrefs() {
		
		e.controls.undoredo.tracksettings = chkbox_undo.isSelected();
		e.controls.undoredo.setHistorySize((int) spinner_undosize.getValue());
		
		e.renderer.targetframerate = (double)((int) spinner_fps.getValue());
		e.renderer.frameduration = 1000/e.renderer.targetframerate;

		e.targetframerate = (double)((int) spinner_fps_1.getValue());
		e.frameduration = 1000/e.targetframerate;
		
		e.renderer.display_relative_voltage = chkbox_potential.isSelected();
		Text.setFontSize((int) spinner_font.getValue());

		e.units = (Units) gui_units.getSelectedItem();
		e.renderer.disp_mat_name = chkbox_matname.isSelected();
		
		e.controls.setvoltage = chkbox_voltage.isSelected();
		
		SemiSim.changeLookAndFeel(e, ((Theme)gui_lookfeel.getSelectedItem()).info);
		
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
		chkbox_undo.setSelected(false);
		chkbox_potential.setSelected(true);
		gui_units.setSelectedItem(Units.SI);
		spinner_undosize.setValue(5);
		spinner_fps.setValue(60);
		spinner_fps_1.setValue(60);
		spinner_font.setValue(12);
		gui_lookfeel.setSelectedItem(Theme.system);
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

				if (version > current_pref_saveversion) {
					fstr.close();
					return;
				}

				if (version == current_pref_saveversion) {
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
						case "units": gui_units.setSelectedItem(gson.fromJson(fstr, Units.class)); break;
						case "fps": this.spinner_fps.setValue(fstr.nextInt()); break;
						case "fps_sim": this.spinner_fps_1.setValue(fstr.nextInt()); break;
						case "fontsize": this.spinner_font.setValue(fstr.nextInt()); break;
						case "undosize": this.spinner_undosize.setValue(fstr.nextInt()); break;
						case "matname": this.chkbox_matname.setSelected(fstr.nextBoolean()); break;
						case "voltage": this.chkbox_voltage.setSelected(fstr.nextBoolean()); break;
						case "theme": gui_lookfeel.setSelectedItem(new Theme(fstr.nextString())); break;
						case "windowstate": e.opts.setExtendedState(fstr.nextInt()); break;
						default: fstr.skipValue();
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

				Gson gson = new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create();

				JsonObject header = new JsonObject();
				header.addProperty("imgsize_x", (int)spinner_imgx.getValue());
				header.addProperty("imgsize_y", (int)spinner_imgy.getValue());
				header.addProperty("undotrackssettings", chkbox_undo.isSelected());
				header.addProperty("potential", chkbox_potential.isSelected());
				header.add("units", gson.toJsonTree((Units) gui_units.getSelectedItem()));
				header.addProperty("fps", (int)spinner_fps.getValue());
				header.addProperty("fps_sim", (int)spinner_fps_1.getValue());
				header.addProperty("fontsize", (int)spinner_font.getValue());
				header.addProperty("undosize", (int)spinner_undosize.getValue());
				header.addProperty("matname", chkbox_matname.isSelected());
				header.addProperty("voltage", chkbox_voltage.isSelected());
				header.addProperty("theme", ((Theme) gui_lookfeel.getSelectedItem()).info.getClassName());
				header.addProperty("windowstate", e.opts.getExtendedState());

				// Version should always be first
				JsonObject save = new JsonObject();
				save.addProperty("version", current_pref_saveversion);
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
		} else if (ev.getSource() == btn_reset) {
			resetPrefs();
		} else if (ev.getSource() == btn_cancel) {
			this.setVisible(false);
		}
	}
	
	public static class Theme {
		public static Theme[] values;
		public static Theme system;
		public static void initThemes() {
			ArrayList<Theme> themes = new ArrayList<Theme>(UIManager.getInstalledLookAndFeels().length);
			system = null;
			for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
				themes.add(new Theme(info));
				if (info.getClassName() == UIManager.getSystemLookAndFeelClassName()) system = new Theme(info);
			}
			values = themes.toArray(new Theme[0]);
		}
		
		LookAndFeelInfo info;
		
		public Theme(LookAndFeelInfo info) {
			this.info = info;
		}
		
		public Theme(String classname) {
			for (Theme theme : values) {
				if (theme.info.getClassName().equals(classname)) {
					info = theme.info;
				}
			}
		}
		
		public String toString() {
			return info.getName();
		}
		
		@Override
		public boolean equals(Object other) {
			if (!(other instanceof Theme)) return false;
			return info.getClassName().equals(((Theme) other).info.getClassName());
		}
	}
}