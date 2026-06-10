// Copyright (c) Brandon Li 2025-2026
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.function.Predicate;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.filechooser.FileFilter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.reflect.TypeToken;
import com.google.gson.stream.JsonReader;

import electrodynamics.Renderer.ScalarMode;
import electrodynamics.Renderer.VectorMode;
import electrodynamics.Simulation.BoundaryCondition;
import electrodynamics.probe.ChargeProbe;
import electrodynamics.probe.CurrentProbe;
import electrodynamics.probe.FluxProbe;
import electrodynamics.probe.Ground;
import electrodynamics.probe.Probe;
import electrodynamics.probe.VoltageProbe;

public class SaveManager {
	Simulation e;
	
	/* Saving and loading */

	public static File infile;
	public static File outfile;
	public static File currentfile;
	public int current_saveversion = 4;
	public String fileextension = ".semisim";
	public String startingpath = ".";

	public String defaultsettings;
	
	public SaveManager(Simulation e) {
		this.e = e;
	}

	public void readFile()
	{
		SwingUtilities.invokeLater(() -> {
			File testfile = new File(startingpath);
			if (!testfile.canRead()) {
				JOptionPane.showMessageDialog(e.opts,
				"Error: Java does not have access to this folder. Please see instructions to fix this issue.");
			}

			JFileChooser fd = new JFileChooser(startingpath);
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
			int result = fd.showOpenDialog(e.opts);
			startingpath = fd.getCurrentDirectory().getPath();

			if (result == JFileChooser.APPROVE_OPTION)
				infile = fd.getSelectedFile();
			else
				infile = null;
			
			readfile(infile);
		});
	}

	@SuppressWarnings("unchecked")
	public void readfile(File infile) {
		e.rwLock.writeLock().lock();
		try {
			if (infile == null || !infile.exists()) return;
			
			if (e.controls.changesmade) {
				String[] options = {"Yes", "No"};
				int result = JOptionPane.showOptionDialog(e.opts, "There are unsaved changes. Do you still wish to open this file?", "Message", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
				if (result != JOptionPane.OK_OPTION)
					return;
			}
			
			try {
				JsonReader fstr = new JsonReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile))));
				Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

				fstr.beginObject();

				assertNextObject(fstr, "version");
				int version = fstr.nextInt();
				
				System.out.println("Loading " + infile.getName() + ", version = " + version);

				if (version > current_saveversion) {
					fstr.close();
					throw new IllegalArgumentException("The file was created in a newer version of SemiSim.");
				}

				JOptionPane optionPane = new JOptionPane("Loading file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
				JDialog dialog = optionPane.createDialog("Loading");

				dialog.setModal(false);
				dialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
				dialog.setVisible(true);

				if (version == 2 || version == 3 || version == 4) {
					setDefaults();
					VectorMode_old vecmode_old = null;
					assertNextObject(fstr, "header");
					fstr.beginObject();
					while (fstr.hasNext()) {
						String name = fstr.nextName();
						switch (name){
						case "resolution": e.default_resolution = fstr.nextInt(); break;
						case "width": e.default_width = fstr.nextDouble(); break;
						case "time": e.time = fstr.nextDouble(); break;
						case "phase": e.AC_phase = fstr.nextDouble(); break;
						case "description": e.description = fstr.nextString(); break;
						
						// Version 2
						case "gui_view": e.controls.scalarview.setOption(gson.fromJson(fstr, Renderer.ScalarView.class)); break;
						case "gui_view_vec": e.controls.vectorview.setOption(gson.fromJson(fstr, Renderer.VectorView.class)); break;
						case "gui_view_vec_mode": vecmode_old = gson.fromJson(fstr, VectorMode_old.class); break;

						// Version 3
						case "scalarview": e.controls.scalarview.setOption(gson.fromJson(fstr, Renderer.ScalarView.class)); break;
						case "vectorview": e.controls.vectorview.setOption(gson.fromJson(fstr, Renderer.VectorView.class)); break;
						case "scalarmode": e.controls.scalarmode.setOption(gson.fromJson(fstr, Renderer.ScalarMode.class)); break;
						case "vectormode": e.controls.vectormode.setOption(gson.fromJson(fstr, Renderer.VectorMode.class)); break;
						case "gui_bc": e.opts.gui_bc.setSelectedItem(gson.fromJson(fstr, BoundaryCondition.class)); break;
						
						default:
							if (e.opts.boolean_names.containsKey(name)) e.opts.boolean_names.get(name).setSelected(fstr.nextBoolean());
							else if (e.opts.integer_names.containsKey(name)) e.opts.integer_names.get(name).setValue(fstr.nextInt());
							break;
						}
					}
					if (vecmode_old != null)
						vecmode_old.applySetting(e);
					fstr.endObject();
					e.opts.setRedundantOptions();
					e.reset(true, null);

					assertNextObject(fstr, "data");
					fstr.beginObject();
					while (fstr.hasNext()) {
						String name = fstr.nextName();
						switch (name){
						case "ex": e.Ex = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "ey": e.Ey = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "hz": e.Hz = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "rho_c": e.rho_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_n": e.rho_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_p": e.rho_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_back": e.rho_back = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_free": e.rho_free = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "jx_c": e.Jx_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_c": e.Jy_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jx_n": e.Jx_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_n": e.Jy_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jx_p": e.Jx_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_p": e.Jy_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "materials": e.materials = validateArraySize((Material[][]) gson.fromJson(fstr, Material[][].class)); break;
						case "materialmap": e.materialmanager.mat_map = (HashMap<Integer, Material>) gson.fromJson(fstr, new TypeToken<HashMap<Integer, Material>>(){}.getType()); break;
						case "last_material_id": e.materialmanager.id_counter = fstr.nextInt(); break;
						
						case "voltageprobes": e.probes.addAll(Arrays.asList(
						(VoltageProbe[]) gson.fromJson(fstr, VoltageProbe[].class))); break;
						case "currentprobes": e.probes.addAll(Arrays.asList(
						(CurrentProbe[]) gson.fromJson(fstr, CurrentProbe[].class))); break;
						case "chargeprobes": e.probes.addAll(Arrays.asList(
						(ChargeProbe[]) gson.fromJson(fstr, ChargeProbe[].class))); break;
						case "fluxprobes": e.probes.addAll(Arrays.asList(
						(FluxProbe[]) gson.fromJson(fstr, FluxProbe[].class))); break;
						case "ground": e.probes.add((Ground) gson.fromJson(fstr, Ground.class)); break;

						default: fstr.skipValue(); break; // skip others
						}
					}
					fstr.endObject();

					for (Probe p : e.probes) p.data.fixWeirdIssue();

					if (testNextObject(fstr, "advsettings")) {
						fstr.beginObject();
						e.adv_opts.readAdvancedSettings(gson, fstr, version < 4? Preset.VERSION_1 : Preset.DEFAULT);
						fstr.endObject();
						e.advsettings_tweaked = true;
					} else {
						e.advsettings_tweaked = false;
					}

					fstr.close();

					e.opts.textPane.setText(e.description);
					e.opts.textPane.setEditable(false);
					e.opts.textPane.setCaretPosition(0);
					e.materialmanager.updateUI();
					if (version < 4) {
						e.initializeAllMaterials();
					}
					e.updateAllMaterials(false);
					e.calcMiscFields(true);
					updateLabels();
					e.controls.undoredo.captureState(e);
				} else if (version == 1) {
					e.reset(true, Preset.VERSION_1);
					
					int view_vec_mode = -1;

					while (fstr.hasNext()) {
						String name = fstr.nextName();
						switch (name){
						case "time": e.time = fstr.nextDouble(); break;
						case "gui_paused": e.opts.gui_paused.setSelected(fstr.nextBoolean()); break;
						case "gui_tooltip": e.opts.menu_tooltip.setSelected(fstr.nextBoolean()); break;
						case "gui_text_bg": e.opts.menu_text_bg.setSelected(fstr.nextBoolean()); break;
						case "gui_view": e.controls.scalarview.setOption(fstr.nextInt()); break;
						case "gui_view_vec": e.controls.vectorview.setOption(fstr.nextInt()); break;
						case "gui_view_vec_mode": view_vec_mode = fstr.nextInt(); break;
						case "gui_simspeed": e.opts.gui_simspeed.setValue(fstr.nextInt()); break;
						case "gui_simspeed_2": e.opts.gui_simspeed_2.setValue(fstr.nextInt()); break;
						case "gui_brightness": e.opts.gui_brightness.setValue(fstr.nextInt()); break;
						case "gui_brightness_vec": e.opts.gui_brightness_vec.setValue(fstr.nextInt()); break;
						case "gui_elem_colors": e.opts.menu_elem_colors.setSelected(fstr.nextBoolean()); break;
						case "gui_bc": e.opts.gui_bc.setSelectedIndex(fstr.nextInt()); break;
						case "description": e.description = fstr.nextString(); break;

						case "ex": e.Ex = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "ey": e.Ey = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "hz": e.Hz = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "rho_c": e.rho_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_n": e.rho_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_p": e.rho_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_back": e.rho_back = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "rho_free": e.rho_free = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "jx_c": e.Jx_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_c": e.Jy_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jx_n": e.Jx_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_n": e.Jy_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jx_p": e.Jx_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
						case "jy_p": e.Jy_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

						case "materials": e.materials = validateArraySize((Material[][]) gson.fromJson(fstr, Material[][].class)); break;

						case "voltageprobes": e.probes.addAll(Arrays.asList(
						(VoltageProbe[]) gson.fromJson(fstr, VoltageProbe[].class))); break;
						case "currentprobes": e.probes.addAll(Arrays.asList(
						(CurrentProbe[]) gson.fromJson(fstr, CurrentProbe[].class))); break;

						case "ground": e.probes.add((Ground) gson.fromJson(fstr, Ground.class)); break;

						default: fstr.skipValue(); break; // skip others
						}
					}

					if (view_vec_mode != -1) {
						VectorMode_old vecmode_old = VectorMode_old.values()[view_vec_mode];
						vecmode_old.applySetting(e);
					}

					e.opts.menu_interface.setSelected(true);
					fstr.endObject();
					fstr.close();

					e.opts.textPane.setText(e.description);
					e.opts.textPane.setEditable(false);
					e.opts.textPane.setCaretPosition(0);
					e.initializeAllMaterials();
					e.updateAllMaterials(false);
					e.calcMiscFields(true);
					updateLabels();
					e.controls.undoredo.captureState(e);
				}

				dialog.dispose();
				e.opts.setTitle(SemiSim.name + " - " + infile.getName());
				e.controls.changesmade = false;
				currentfile = infile;
			} catch (FileNotFoundException ex) {
				return;
			} catch (IOException | IllegalArgumentException ex) {
				JOptionPane.showMessageDialog(e.opts,
				"Unable to load file.\n" + ex.getMessage());
				ex.printStackTrace();
				return;
			}
			return;
		} finally {
			e.rwLock.writeLock().unlock();
		}
	}
	
	public void updateLabels() {
		for (Probe p : e.probes)
			if (p.labelcoord.x == -1) p.calculateDefaultLabelCoords();
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

	public Material[][] validateArraySize(Material[][] array) throws RuntimeException {
		Material[][] field = new Material[e.nx][e.ny];

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				if (i < e.nx && j < e.ny) {
					field[i][j] = array[i][j];
				}
			}
		}

		for (int i = 0; i < e.nx; i++) {
			for (int j = 0; j < e.ny; j++) {
				if (field[i][j] == null)
					field[i][j] = new Material();
			}
		}

		return field;
	}

	public double[][] validateArraySize(double[][] array) throws RuntimeException {
		double[][] field = new double[e.nx][e.ny];

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				if (i < e.nx && j < e.ny) {
					field[i][j] = array[i][j];
				}
			}
		}

		return field;
	}

	public void writeFile(boolean saveas)
	{
		SwingUtilities.invokeLater(() -> {
			if (!saveas && currentfile != null && currentfile.exists()) {
				writeFile(currentfile);
				return;
			}
			
			File testfile = new File(startingpath);
			if (!testfile.canWrite()) {
				JOptionPane.showMessageDialog(e.opts,
				"Error: Java does not have access to this folder. Please see instructions to fix this issue.");
			}

			JFileChooser fd = new JFileChooser(startingpath);
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
			int result = fd.showSaveDialog(e.opts);
			startingpath = fd.getCurrentDirectory().getPath();

			if (result == JFileChooser.APPROVE_OPTION)
				outfile = fd.getSelectedFile();
			else
				outfile = null;

			if (outfile == null) return;
			if (!outfile.getName().endsWith(fileextension))
				outfile = new File(outfile.getAbsolutePath() + fileextension);

			if (outfile.exists()) {
				String[] options = {"Yes", "No"};

				result = JOptionPane.showOptionDialog(e.opts, "A file named " + outfile.getName() + " already exists. Do you wish to overwrite it?", "Message", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
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
				PrintWriter fstr = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outfile)));

				Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

				JOptionPane optionPane = new JOptionPane("Saving file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
				JDialog dialog = optionPane.createDialog("Saving");

				dialog.setModal(false);
				dialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
				dialog.setVisible(true);

				JsonObject header = new JsonObject();
				header.addProperty("resolution", e.default_resolution);
				header.addProperty("width", e.default_width);
				header.addProperty("time", e.time);
				header.addProperty("phase", e.AC_phase);
				
				for (String s : e.opts.boolean_names.keySet()) header.addProperty(s, e.opts.boolean_names.get(s).isSelected());
				for (String s : e.opts.integer_names.keySet()) header.addProperty(s, e.opts.integer_names.get(s).getValue());
				
				header.addProperty("description", e.description);
				header.add("scalarview", gson.toJsonTree(e.controls.scalarview.getOption()));
				header.add("vectorview", gson.toJsonTree(e.controls.vectorview.getOption()));
				header.add("scalarmode", gson.toJsonTree(e.controls.scalarmode.getOption()));
				header.add("vectormode", gson.toJsonTree(e.controls.vectormode.getOption()));
				header.add("gui_bc", gson.toJsonTree(e.opts.gui_bc.getSelectedItem()));

				JsonObject data = new JsonObject();
				data.add("ex", gson.toJsonTree(e.Ex));
				data.add("ey", gson.toJsonTree(e.Ey));
				data.add("hz", gson.toJsonTree(e.Hz));
				data.add("rho_c", gson.toJsonTree(e.rho_abs));
				data.add("rho_n", gson.toJsonTree(e.rho_n));
				data.add("rho_p", gson.toJsonTree(e.rho_p));
				data.add("rho_back", gson.toJsonTree(e.rho_back));
				data.add("rho_free", gson.toJsonTree(e.rho_free));
				data.add("jx_c", gson.toJsonTree(e.Jx_abs));
				data.add("jy_c", gson.toJsonTree(e.Jy_abs));
				data.add("jx_n", gson.toJsonTree(e.Jx_n));
				data.add("jy_n", gson.toJsonTree(e.Jy_n));
				data.add("jx_p", gson.toJsonTree(e.Jx_p));
				data.add("jy_p", gson.toJsonTree(e.Jy_p));
				data.add("materials", gson.toJsonTree(e.materials));
				data.add("voltageprobes", gson.toJsonTree(filterByType(e.probes, (p)->!(p instanceof VoltageProbe && !(p instanceof Ground))).toArray()));
				data.add("currentprobes", gson.toJsonTree(filterByType(e.probes, (p)->!(p instanceof CurrentProbe)).toArray()));
				data.add("chargeprobes", gson.toJsonTree(filterByType(e.probes, (p)->!(p instanceof ChargeProbe)).toArray()));
				data.add("fluxprobes", gson.toJsonTree(filterByType(e.probes, (p)->!(p instanceof FluxProbe)).toArray()));
				data.add("ground", gson.toJsonTree(e.getGround()));
				data.add("materialmap", gson.toJsonTree(e.materialmanager.mat_map));
				data.addProperty("last_material_id", e.materialmanager.id_counter);

				JsonObject advsettings = new JsonObject();
				e.adv_opts.writeAdvancedSettings(gson, advsettings);

				// Version should always be first
				JsonObject save = new JsonObject();
				save.addProperty("version", current_saveversion);
				save.add("header", header);
				save.add("data", data);
				save.add("advsettings", advsettings);

				String json = gson.toJson(save);

				fstr.print(json);
				fstr.flush();
				fstr.close();

				dialog.dispose();
				e.opts.setTitle(SemiSim.name + " - " + outfile.getName());
				e.controls.changesmade = false;
				currentfile = outfile;
			} catch (FileNotFoundException e) {
				return;
			} catch (IOException e) {
				e.printStackTrace();
				return;
			}
			return;
		} finally {
			e.rwLock.writeLock().unlock();
		}
	}
	
	public ArrayList<Probe> filterByType(List<Probe> list, Predicate<? super Probe> filter) {
		ArrayList<Probe> listcopy = new ArrayList<Probe>(list);
		listcopy.removeIf(filter);
		return listcopy;
	}
	
	public void setDefaults() {
		e.setDefaultParameters();
		e.opts.setDefaults(e);
	}
	
	// Old version of VectorMode for backwards file compatibility
	enum VectorMode_old {

		ARROWS,
		LINES,
		DOTS,
		CONTOUR,
		SPECIES;
		
		public void applySetting(Simulation e) {
			switch(this) {
			case ARROWS:
				e.controls.vectormode.setOption(VectorMode.ARROWS);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case CONTOUR:
				e.controls.vectormode.setOption(VectorMode.NONE);
				e.controls.scalarmode.setOption(ScalarMode.CONTOUR_COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case DOTS:
				e.controls.vectormode.setOption(VectorMode.DOTS);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case LINES:
				e.controls.vectormode.setOption(VectorMode.LINES);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(false);
				break;
			case SPECIES:
				e.controls.vectormode.setOption(VectorMode.NONE);
				e.controls.scalarmode.setOption(ScalarMode.COLORS);
				e.opts.gui_carriers.setSelected(true);
				e.opts.gui_carrier_density.setValue(e.opts.gui_brightness_vec.getValue());
				break;
			default:
				break;
			}
		}
	}
}
