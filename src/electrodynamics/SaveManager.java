// Copyright (c) Brandon Li 2025
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
import java.util.Arrays;
import java.util.concurrent.CopyOnWriteArrayList;
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
import com.google.gson.JsonIOException;
import com.google.gson.JsonObject;
import com.google.gson.JsonSyntaxException;
import com.google.gson.stream.JsonReader;

public class SaveManager {
	Simulation e;
	
	/* Saving and loading */

	public static File infile;
	public static File outfile;
	int saveversion = 2;
	String fileextension = ".semisim";
	String startingpath = ".";
	
	public SaveManager(Simulation e) {
		this.e = e;
	}

	public void readFile()
	{
		SwingUtilities.invokeLater(() -> {
			e.rwLock.writeLock().lock();
			try {
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

				if (infile == null || !infile.exists()) return;
				try {
					JsonReader fstr = new JsonReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile))));
					Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

					fstr.beginObject();

					assertNextObject(fstr, "version");
					int version = fstr.nextInt();

					if (version > saveversion) {
						fstr.close();
						throw new IllegalArgumentException("The file was created in a newer version of SemiSim.");
					}

					JOptionPane optionPane = new JOptionPane("Loading file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
					JDialog dialog = optionPane.createDialog("Loading");

					dialog.setModal(false);
					dialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
					dialog.setVisible(true);

					if (version == saveversion) {

						int resolution_tmp = e.default_resolution;
						double width_tmp = e.default_width;
						assertNextObject(fstr, "header");
						fstr.beginObject();
						while (fstr.hasNext()) {
							String name = fstr.nextName();
							switch (name){
							case "resolution": resolution_tmp = fstr.nextInt(); break;
							case "width": width_tmp = fstr.nextDouble(); break;
							case "time": e.time = fstr.nextDouble(); break;
							case "gui_paused": e.opts.gui_paused.setSelected(fstr.nextBoolean()); break;
							case "gui_tooltip": e.opts.gui_tooltip.setSelected(fstr.nextBoolean()); break;
							case "gui_text_bg": e.opts.gui_text_bg.setSelected(fstr.nextBoolean()); break;
							case "gui_elem_colors": e.opts.gui_elem_colors.setSelected(fstr.nextBoolean()); break;
							case "gui_interface": e.opts.gui_interface.setSelected(fstr.nextBoolean()); break;
							case "gui_simspeed": e.opts.gui_simspeed.setValue(fstr.nextInt()); break;
							case "gui_simspeed_2": e.opts.gui_simspeed_2.setValue(fstr.nextInt()); break;
							case "gui_brightness": e.opts.gui_brightness.setValue(fstr.nextInt()); break;
							case "gui_brightness_vec": e.opts.gui_brightness_vec.setValue(fstr.nextInt()); break;
							case "description": e.opts.textPane.setText(fstr.nextString()); break;
							case "gui_view": e.opts.gui_view.setSelectedItem(gson.fromJson(fstr, Renderer.ScalarView.class)); break;
							case "gui_view_vec": e.opts.gui_view_vec.setSelectedItem(gson.fromJson(fstr, Renderer.VectorView.class)); break;
							case "gui_view_vec_mode": e.opts.gui_view_vec_mode.setSelectedItem(gson.fromJson(fstr, Renderer.VectorMode.class)); break;
							case "gui_bc": e.opts.gui_bc.setSelectedItem(gson.fromJson(fstr, BoundaryCondition.class)); break;
							default: fstr.skipValue(); break; // skip others
							}
						}
						fstr.endObject();

						e.setSize(resolution_tmp, width_tmp);
						e.resetFields(true);

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

							case "voltageprobes": e.voltageprobes = new CopyOnWriteArrayList<>(Arrays.asList(
							(VoltageProbe[]) gson.fromJson(fstr, VoltageProbe[].class))); break;
							case "currentprobes": e.currentprobes = new CopyOnWriteArrayList<>(Arrays.asList(
							(CurrentProbe[]) gson.fromJson(fstr, CurrentProbe[].class))); break;

							case "ground": e.ground = (VoltageProbe) gson.fromJson(fstr, VoltageProbe.class); break;

							default: fstr.skipValue(); break; // skip others
							}
						}
						fstr.endObject();

						if (testNextObject(fstr, "advsettings")) {
							fstr.beginObject();
							readAdvancedSettings(gson, fstr);
							fstr.endObject();
							e.advsettings_tweaked = true;
						} else {
							e.advsettings_tweaked = false;
						}

						fstr.close();

						e.opts.textPane.setEditable(false);
						e.opts.textPane.setCaretPosition(0);
						e.constructBoundary();
						e.updateAllMaterials(false);
						e.calcMiscFields(true);
					} else if (version == 1) {

						e.setSize(e.default_resolution, e.default_width);
						e.resetFields(true);

						while (fstr.hasNext()) {
							String name = fstr.nextName();
							switch (name){
							case "time": e.time = fstr.nextDouble(); break;
							case "gui_paused": e.opts.gui_paused.setSelected(fstr.nextBoolean()); break;
							case "gui_tooltip": e.opts.gui_tooltip.setSelected(fstr.nextBoolean()); break;
							case "gui_text_bg": e.opts.gui_text_bg.setSelected(fstr.nextBoolean()); break;
							case "gui_view": e.opts.gui_view.setSelectedIndex(fstr.nextInt()); break;
							case "gui_view_vec": e.opts.gui_view_vec.setSelectedIndex(fstr.nextInt()); break;
							case "gui_view_vec_mode": e.opts.gui_view_vec_mode.setSelectedIndex(fstr.nextInt()); break;
							case "gui_simspeed": e.opts.gui_simspeed.setValue(fstr.nextInt()); break;
							case "gui_simspeed_2": e.opts.gui_simspeed_2.setValue(fstr.nextInt()); break;
							case "gui_brightness": e.opts.gui_brightness.setValue(fstr.nextInt()); break;
							case "gui_brightness_vec": e.opts.gui_brightness_vec.setValue(fstr.nextInt()); break;
							case "gui_elem_colors": e.opts.gui_elem_colors.setSelected(fstr.nextBoolean()); break;
							case "gui_bc": e.opts.gui_bc.setSelectedIndex(fstr.nextInt()); break;
							case "description": e.opts.textPane.setText(fstr.nextString()); break;

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

							case "voltageprobes": e.voltageprobes = new CopyOnWriteArrayList<>(Arrays.asList(
							(VoltageProbe[]) gson.fromJson(fstr, VoltageProbe[].class))); break;
							case "currentprobes": e.currentprobes = new CopyOnWriteArrayList<>(Arrays.asList(
							(CurrentProbe[]) gson.fromJson(fstr, CurrentProbe[].class))); break;

							case "ground": e.ground = (VoltageProbe) gson.fromJson(fstr, VoltageProbe.class); break;

							default: fstr.skipValue(); break; // skip others
							}
						}
						e.opts.gui_interface.setSelected(true);
						fstr.endObject();
						fstr.close();

						e.opts.textPane.setEditable(false);
						e.opts.textPane.setCaretPosition(0);
						e.constructBoundary();
						e.updateAllMaterials(false);
						e.calcMiscFields(true);
					}

					dialog.dispose();
					e.opts.setTitle("Brandon's semiconductor simulator - " + infile.getName());
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
		});
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

	public void writeFile()
	{
		SwingUtilities.invokeLater(() -> {
			e.rwLock.writeLock().lock();
			try {
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
					result = JOptionPane.showConfirmDialog(e.opts, "A file with that name already exists. Do you wish to overwrite it?", "Message", JOptionPane.YES_NO_OPTION);
					if (result != JOptionPane.OK_OPTION)
						return;
				}

				try {
					PrintWriter fstr = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outfile)));

					Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();

					JOptionPane optionPane = new JOptionPane("Saving file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
					JDialog dialog = optionPane.createDialog("Saving");

					dialog.setModal(false);
					dialog.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
					dialog.setVisible(true);

					JsonObject header = new JsonObject();
					header.addProperty("resolution", e.resolution);
					header.addProperty("width", e.width);
					header.addProperty("time", e.time);
					header.addProperty("gui_paused", e.opts.gui_paused.isSelected());
					header.addProperty("gui_tooltip", e.opts.gui_tooltip.isSelected());
					header.addProperty("gui_text_bg", e.opts.gui_text_bg.isSelected());
					header.addProperty("gui_elem_colors", e.opts.gui_elem_colors.isSelected());
					header.addProperty("gui_interface", e.opts.gui_interface.isSelected());
					header.addProperty("gui_simspeed", e.opts.gui_simspeed.getValue());
					header.addProperty("gui_simspeed_2", e.opts.gui_simspeed_2.getValue());
					header.addProperty("gui_brightness", e.opts.gui_brightness.getValue());
					header.addProperty("gui_brightness_vec", e.opts.gui_brightness_vec.getValue());
					header.addProperty("description", e.opts.textPane.getText());
					header.add("gui_view", gson.toJsonTree(e.opts.gui_view.getSelectedItem()));
					header.add("gui_view_vec", gson.toJsonTree(e.opts.gui_view_vec.getSelectedItem()));
					header.add("gui_view_vec_mode", gson.toJsonTree(e.opts.gui_view_vec_mode.getSelectedItem()));
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
					data.add("voltageprobes", gson.toJsonTree(e.voltageprobes.toArray(new VoltageProbe[e.voltageprobes.size()])));
					data.add("currentprobes", gson.toJsonTree(e.currentprobes.toArray(new CurrentProbe[e.currentprobes.size()])));
					data.add("ground", gson.toJsonTree(e.ground));

					JsonObject advsettings = new JsonObject();
					writeAdvancedSettings(gson, advsettings);

					// Version should always be first
					JsonObject save = new JsonObject();
					save.addProperty("version", saveversion);
					save.add("header", header);
					save.add("data", data);
					save.add("advsettings", advsettings);

					String json = gson.toJson(save);

					fstr.print(json);
					fstr.flush();
					fstr.close();

					dialog.dispose();
					e.opts.setTitle("Brandon's semiconductor simulator - " + outfile.getName());
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
		});
	}

	public void writeAdvancedSettings() {
		e.adv_opts.width				.setText(formatDouble(e.width						));
		e.adv_opts.resolution			.setText(Integer.toString(e.resolution				));
		e.adv_opts.depth				.setText(formatDouble(e.depth						));
		e.adv_opts.mu_electron			.setText(formatDouble(e.mu_electron					));
		e.adv_opts.mu_hole				.setText(formatDouble(e.mu_hole						));
		e.adv_opts.ni_semi				.setText(formatDouble(e.ni_semi						));
		e.adv_opts.W_semi				.setText(formatDouble(e.W_semi/e.eVtoJ				));
		e.adv_opts.E_b_semi				.setText(formatDouble(e.E_b_semi/e.eVtoJ			));
		e.adv_opts.ni_metal				.setText(formatDouble(e.ni_metal					));
		e.adv_opts.W_metal				.setText(formatDouble(e.W_metal_default/e.eVtoJ		));
		e.adv_opts.E_b_metal			.setText(formatDouble(e.E_b_metal/e.eVtoJ			));
		e.adv_opts.W_metal_high			.setText(formatDouble(e.W_metal_high/e.eVtoJ		));
		e.adv_opts.W_metal_low			.setText(formatDouble(e.W_metal_low/e.eVtoJ			));
		e.adv_opts.recomb_rate_semi		.setText(formatDouble(e.recomb_rate_semi			));
		e.adv_opts.recomb_rate_metal	.setText(formatDouble(e.recomb_rate_metal			));
		e.adv_opts.T					.setText(formatDouble(e.T							));
	}
	
	public void readAdvancedSettings() {
		
		
		double width_tmp			= Double.valueOf(e.adv_opts.width				.getText());	
		int resolution_tmp			= Integer.valueOf(e.adv_opts.resolution			.getText());
		e.depth						= Double.valueOf(e.adv_opts.depth				.getText());
		e.mu_electron				= Double.valueOf(e.adv_opts.mu_electron			.getText());
		e.mu_hole					= Double.valueOf(e.adv_opts.mu_hole				.getText());
		e.ni_semi					= Double.valueOf(e.adv_opts.ni_semi				.getText());
		e.W_semi					= Double.valueOf(e.adv_opts.W_semi				.getText())*e.eVtoJ;
		e.E_b_semi					= Double.valueOf(e.adv_opts.E_b_semi			.getText())*e.eVtoJ;
		e.ni_metal					= Double.valueOf(e.adv_opts.ni_metal			.getText());
		e.W_metal_default			= Double.valueOf(e.adv_opts.W_metal				.getText())*e.eVtoJ;
		e.E_b_metal					= Double.valueOf(e.adv_opts.E_b_metal			.getText())*e.eVtoJ;
		e.W_metal_high				= Double.valueOf(e.adv_opts.W_metal_high		.getText())*e.eVtoJ;
		e.W_metal_low				= Double.valueOf(e.adv_opts.W_metal_low			.getText())*e.eVtoJ;
		e.recomb_rate_semi			= Double.valueOf(e.adv_opts.recomb_rate_semi	.getText());
		e.recomb_rate_metal			= Double.valueOf(e.adv_opts.recomb_rate_metal	.getText());
		e.T							= Double.valueOf(e.adv_opts.T					.getText());

		e.updateConstants();
		e.setSize(resolution_tmp, width_tmp);
		e.lastsimspeed = -1;
		e.advsettings_tweaked = true;
	}

	public void writeAdvancedSettings (Gson gson, JsonObject advsettings) {
		advsettings.addProperty("depth", e.depth 						);
		advsettings.addProperty("mu_electron", e.mu_electron			);
		advsettings.addProperty("mu_hole", e.mu_hole					);
		advsettings.addProperty("ni_semi", e.ni_semi					);
		advsettings.addProperty("W_semi", e.W_semi						);
		advsettings.addProperty("E_b_semi", e.E_b_semi					);
		advsettings.addProperty("ni_metal", e.ni_metal					);
		advsettings.addProperty("W_metal_default", e.W_metal_default		);
		advsettings.addProperty("E_b_metal", e.E_b_metal					);
		advsettings.addProperty("W_metal_high", e.W_metal_high				);
		advsettings.addProperty("W_metal_low", e.W_metal_low				);
		advsettings.addProperty("recomb_rate_semi", e.recomb_rate_semi		);
		advsettings.addProperty("recomb_rate_metal", e.recomb_rate_metal	);
		advsettings.addProperty("T", e.T									);
	}
	

	public void readAdvancedSettings (Gson gson, JsonReader fstr) throws JsonIOException, JsonSyntaxException, IOException, RuntimeException {

		while (fstr.hasNext()) {
			String name = fstr.nextName();
			switch (name){
			case "depth": e.depth 						= fstr.nextDouble(); break;
			case "mu_electron": e.mu_electron			= fstr.nextDouble(); break;
			case "mu_hole": e.mu_hole					= fstr.nextDouble(); break;
			case "ni_semi": e.ni_semi					= fstr.nextDouble(); break;
			case "W_semi": e.W_semi						= fstr.nextDouble(); break;
			case "E_b_semi": e.E_b_semi					= fstr.nextDouble(); break;
			case "ni_metal": e.ni_metal					= fstr.nextDouble(); break;
			case "W_metal_default": e.W_metal_default		= fstr.nextDouble(); break;
			case "E_b_metal": e.E_b_metal					= fstr.nextDouble(); break;
			case "W_metal_high": e.W_metal_high				= fstr.nextDouble(); break;
			case "W_metal_low": e.W_metal_low				= fstr.nextDouble(); break;
			case "recomb_rate_semi": e.recomb_rate_semi		= fstr.nextDouble(); break;
			case "recomb_rate_metal": e.recomb_rate_metal	= fstr.nextDouble(); break;
			case "T": e.T									= fstr.nextDouble(); break;

			default: fstr.skipValue(); break; // skip others
			}
		}
	}

	public String formatDouble(double d) {
		return Double.toString(d);
	}
}
