// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Timer;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileSystemView;

import com.formdev.flatlaf.FlatDarkLaf;
import com.formdev.flatlaf.FlatLaf;
import com.formdev.flatlaf.FlatLightLaf;
import com.formdev.flatlaf.intellijthemes.FlatHighContrastIJTheme;
import com.formdev.flatlaf.intellijthemes.FlatSolarizedLightIJTheme;
import com.formdev.flatlaf.intellijthemes.materialthemeuilite.FlatMTMaterialDarkerIJTheme;
import electrodynamics.Simulation.SimulationThread;
import electrodynamics.gui.Preferences.Theme;
import electrodynamics.plot.Plot;

public class SemiSim {
	
	public static SemiSim instance;
	public static int n_threads = Runtime.getRuntime().availableProcessors();
	public static String name = "Brandon's semiconductor simulator";
	public static String about = "<html><body><p style='width: 250px;'>Brandon's Semiconductor Simulator / SemiSim.<br>"
									+ "Version 2.0.2<br>"
									+ "(c) 2026 Brandon Li<br><br>"
									+ "Thanks to Paul Falstad, Ariel Baksh, and retconaway for providing help, feedback, and suggestions.<br><br>"
									+ "Data taken from:<br>"
									+ "Sitlisky, Vadim. &ldquo;New Semiconductor Materials. Characteristics and Properties&rdquo;. <a href=\"http://www.ioffe.ru\"><em>www.ioffe.ru</em></a>.<br> Retrieved June 2026.<br>"
									+ "Schroder, D. K. (2006). <em>Semiconductor material<br> and device characterization</em>. John Wiley &amp; Sons.</p></body></html>";

	public static Path rootdir = Paths.get(".");
	public static Path userdir = Paths.get(".");
	public static OS os;
	ArrayList<SimulationThread> sim_threads = new ArrayList<>();
	ArrayList<Renderer.GraphicsThread> graphics_threads = new ArrayList<>();
	Timer master_timer = new Timer();
	Timer graphics_timer = new Timer();
	Timer misc_timer = new Timer();
	ScheduledThreadPoolExecutor threadPool = new ScheduledThreadPoolExecutor(3);
	Simulation sim;
	
	public SemiSim() {
		sim = new Simulation();
		
		for (int i = 0; i < n_threads; i++) {
			sim_threads.add(sim.new SimulationThread(i, n_threads));
		}

		for (int i = 0; i < n_threads; i++) {
			graphics_threads.add(sim.renderer.new GraphicsThread(i, n_threads));
		}

		for (int i = 0; i < n_threads; i++) {
			sim_threads.get(i).start();
			graphics_threads.get(i).start();
		}

		threadPool.schedule(sim, 0, TimeUnit.MILLISECONDS);
		threadPool.schedule(sim.renderer, 0, TimeUnit.MILLISECONDS);
		threadPool.schedule(sim.potentialSolver, 0, TimeUnit.MILLISECONDS);
	}
	
	public static File getRootFile(String path) {
		return rootdir.resolve(Paths.get(path)).toFile();
	}

	public static File getUserFile(String path) {
		return userdir.resolve(Paths.get(path)).toFile();
	}

	public static void displayErrorMessage(Exception e) {
		if (instance == null) return;
		try {
			Runnable r = () -> {
				JOptionPane.showMessageDialog(instance.sim.opts, e.toString(), "Fatal error!", JOptionPane.OK_OPTION);
				e.printStackTrace();
				try {
					PrintWriter pw = new PrintWriter(new FileOutputStream(getUserFile("error_log.txt")));
					e.printStackTrace(pw);
					pw.flush();
					pw.close();
				} catch (FileNotFoundException e1) {
					System.exit(-1);
				}     
				System.exit(-1);
			};

			if (SwingUtilities.isEventDispatchThread()) {
				r.run();
			} else {
				SwingUtilities.invokeAndWait(r);
			}
		} catch (InvocationTargetException | InterruptedException e1) {
			e1.printStackTrace();
		}
	}
	
	public static void displayWarningMessage(String title, String description) {
		if (instance == null) return;
		try {
			Runnable r = () -> {
				JOptionPane.showMessageDialog(instance.sim.opts, description, title, JOptionPane.OK_OPTION);
			};

			if (SwingUtilities.isEventDispatchThread()) {
				r.run();
			} else {
				SwingUtilities.invokeAndWait(r);
			}
		} catch (InvocationTargetException | InterruptedException e1) {
			e1.printStackTrace();
		}
	}

	public static void detect64Bit() {
		if (instance == null) return;
		if (!System.getProperty("sun.arch.data.model").equals("64"))
		{
			int result = JOptionPane.showConfirmDialog(instance.sim.opts, "Running this application on a 32-bit platform may cause some issues. Do you still wish to proceed?", "Message", JOptionPane.YES_NO_OPTION);
			if (result != JOptionPane.OK_OPTION)
			{
				System.exit(0);
			}
		}
	}
	
	public static void detectOS() {
		String osname = System.getProperty("os.name").toLowerCase();
		
		if (osname.contains("windows")) {
			os = OS.WINDOWS;
		} else if (osname.contains("mac")) {
			os = OS.MAC;
		} else if (osname.contains("linux")) {
			os = OS.LINUX;
		} else {
			os = OS.UNKNOWN;
		}
	}

	public static void setDirectory() {
		try {
			rootdir = Paths.get(SemiSim.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getParent();
		} catch (URISyntaxException e1) {
			e1.printStackTrace();
		}

		userdir = new JFileChooser().getFileSystemView().getDefaultDirectory().toPath();
		System.out.println(FileSystemView.getFileSystemView().getDefaultDirectory().toPath().toString());
		
		if (os == OS.WINDOWS) {
			userdir = userdir.resolve(Paths.get("SemiSim"));
		} else if (os == OS.MAC) {
			userdir = userdir.resolve(Paths.get("Documents/SemiSim"));
		} else if (os == OS.LINUX) {
			userdir = userdir.resolve(Paths.get("SemiSim"));
		} else {
			userdir = userdir.resolve(Paths.get("SemiSim"));
		}
		
		if (!userdir.toFile().exists()) {
			userdir.toFile().mkdir();
		}

		try {
			File mat = getUserFile("materials");
			if (!mat.exists()) {
				mat.mkdir();
				for (File f : getRootFile("materials").listFiles()) {
					File dest = getUserFile("materials/" + f.getName());
					if (!dest.exists()) {
						Files.copy(f.toPath(), dest.toPath());
					}
				}
			}

			File screenshots = getUserFile("screenshots");
			if (!screenshots.exists()) {
				screenshots.mkdir();
			}
			
			File simulations = getUserFile("simulations");
			if (!simulations.exists()) {
				simulations.mkdir();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void initializeLookAndFeel() {
		try {
			UIManager.setLookAndFeel(
					UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}
		
		installLaf(new FlatLightLaf());
		installLaf(new FlatDarkLaf());
		installLaf(new FlatSolarizedLightIJTheme());
		//installLaf(new FlatSolarizedDarkIJTheme());
		//installLaf(new FlatMTMaterialLighterIJTheme());
		installLaf(new FlatMTMaterialDarkerIJTheme());
		installLaf(new FlatHighContrastIJTheme());
		
		Theme.initThemes();
	}
	
	public static <T extends LookAndFeel> void installLaf(T t) {
		FlatLaf.installLafInfo(t.getName(), t.getClass());
	}
	

	public static void changeLookAndFeel(Simulation sim, LookAndFeelInfo info) {
		if (UIManager.getLookAndFeel().getClass().getName().equals(info.getClassName())) return;
		try {
			UIManager.setLookAndFeel(info.getClassName());
			SwingUtilities.updateComponentTreeUI(sim.opts);
			SwingUtilities.updateComponentTreeUI(sim.adv_opts);
			SwingUtilities.updateComponentTreeUI(sim.materialmanager);
			SwingUtilities.updateComponentTreeUI(sim.materialviewer);
			SwingUtilities.updateComponentTreeUI(sim.prefs);

			for (Plot p : sim.plots) SwingUtilities.updateComponentTreeUI(p.frame);
			
			if (Steam.downloadui != null) SwingUtilities.updateComponentTreeUI(Steam.downloadui);
			if (Steam.uploadui != null) SwingUtilities.updateComponentTreeUI(Steam.uploadui);
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args)
	{
		detectOS();
		
		setDirectory();
		
		initializeLookAndFeel();
		
		Steam.initialize();
		
		SwingUtilities.invokeLater(() -> {
			instance = new SemiSim();
		});
	}
	
	public enum OS {
		WINDOWS, MAC, LINUX, UNKNOWN;
	}
}
