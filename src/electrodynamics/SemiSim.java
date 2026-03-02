// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Timer;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

import electrodynamics.Simulation.SimulationThread;

public class SemiSim {
	
	public static SemiSim instance;
	public static int n_threads = Runtime.getRuntime().availableProcessors();
	public static String name = "Brandon's semiconductor simulator";
	public static String about = "<html><body><p style='width: 250px;'>Brandon's Semiconductor Simulator / SemiSim.<br>"
									+ "Version 1.12<br>"
									+ "(c) 2026 Brandon Li<br><br>"
									+ "Thanks to Paul Falstad, Ariel Baksh, and retconaway for providing help, feedback, and suggestions.</p></body></html>";

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
			sim_threads.add(sim.new SimulationThread(i, n_threads, sim.nx));
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
		//master_timer.schedule(sim, 0, sim.renderer.frameduration);
		//graphics_timer.schedule(sim.renderer, 0, sim.renderer.frameduration);
		//misc_timer.schedule(sim.potentialSolver, 0, sim.renderer.frameduration);
	}
	
	public static void displayErrorMessage(Exception e) {
		if (instance == null) return;
		try {
			SwingUtilities.invokeAndWait(() -> {
				JOptionPane.showMessageDialog(instance.sim.opts, e.toString(), "Error", JOptionPane.OK_OPTION);
				e.printStackTrace();
				try {
					PrintWriter pw = new PrintWriter(new FileOutputStream("error_log.txt"));
				    e.printStackTrace(pw);
				    pw.flush();
				    pw.close();
				} catch (FileNotFoundException e1) {
					System.exit(-1);
				}     
				System.exit(-1);
			});
		} catch (InvocationTargetException | InterruptedException e1) {
			System.exit(-1);
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
	
	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(
					UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}
		
		SwingUtilities.invokeLater(() -> {
			instance = new SemiSim();
		});
	}
}
