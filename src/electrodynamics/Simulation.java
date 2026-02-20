// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
import java.awt.BorderLayout;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.TimerTask;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import javax.swing.ButtonGroup;
import javax.swing.InputMap;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import electrodynamics.Controls.Brush;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.plot.BandPlot;
import electrodynamics.plot.CarrierPlot;
import electrodynamics.plot.Plot;
import electrodynamics.plot.ScalarPlot;
import electrodynamics.util.MenuBuilder;
import electrodynamics.util.PeriodicTask;
import electrodynamics.util.Timer;
import electrodynamics.util.Utils;

public class Simulation extends PeriodicTask {
	/*
	 * Demos:
	 * 	- PN diode (LED)
	 * 	- Schottky diode
	 * 	- MOSFET (p-channel, n-channel), depletion & enhancement
	 * 	- BJT (NPN, PNP)
	 * 	- JFET (p-channel, n-channel)
	 *  - IGBT - doesn't work well
	 *  - MESFET
	 *  - SCR, maybe
	 *  - UJT
	 *  - Darlington pair
	 */
	
	/* Parts */
	
	public RenderCanvas canvas;
	public MainWindow opts;
	public AdvancedOptions adv_opts;
	public Renderer renderer;
	public Controls controls;
	public SaveManager savemanager;
	
	public ArrayList<Plot> plots = new ArrayList<>();
	public BandPlot bandplot;
	public ScalarPlot scalarplot;
	public CarrierPlot carrierplot;
	
	
	/* Multithreading */
	
	ReentrantReadWriteLock rwLock = new ReentrantReadWriteLock();
	ReentrantLock poissonLock = new ReentrantLock(true);

	CyclicBarrier start_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	CyclicBarrier stop_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	CyclicBarrier mid_barrier = new CyclicBarrier(SemiSim.n_threads);

	
	/* Domain parameters */

	public int default_resolution = 256;
	public double default_width = 2.56e-5;
	
	public int resolution = 0;
	public int nx;						// Number of grid popublic ints in x-dimension
	public int ny;						// Number of grid popublic ints in y-dimension
	public double width;				// Width, in SI
	public double ds;					// Spatial discretization
	public double dt;					// Timestep
	public double depth = 1e-3;		// Extent of circuit in z-dimension, only used to get reasonable values for current probe
	
	public int absorber_width;
	public double absorbing_coeff;
	public double Hz_dissipation;
	public double dt_maximum;
	public int parity = -1;			// Negative sign resulting from flipped y-axis in graphics coordinate system

	public double time = 0.0;
	public long stepnumber = 0;
	public long frame = 0;
	public int lastsimspeed = 0;
	public int iteration_multiplier = 0;

	public double error_detection_threshold = 1e-5;
	public boolean sign_violation = false;
	public boolean numerical_overflow = false;
	
	public boolean advsettings_tweaked = false;
	
	public double AC_phase;
	public double AC_freq;
	public double AC_amplitude;
	public boolean AC_source_exists;

	
	/* Physical constants */

	public double eps0 = 8.85e-12;
	public double mu0 = 1.257e-6;
	public double c = 1/Math.sqrt(eps0*mu0);
	public double k = 1.381e-23;
	public double e_charge = 1.6e-19;
	public double e_mass = 9e-31;
	public double eVtoJ = e_charge;			// eV to J conversion factor
	public double q_n = -e_charge;				// Charge of single electron
	public double q_p = e_charge;				// Charge of single hole

	
	/* Material properties */

	public double T = 297.61;									// System temperature
	public double beta = 1/(k*T);
	
	public double mu_electron = 0.1400*1500;					// Electron mobility
	public double D_electron = mu_electron/(beta*e_charge);	// Diffusion constant, determined by Einstein relation

	public double mu_hole = 0.5*mu_electron;					// Hole mobility, slightly lower than electron
	public double D_hole = mu_hole/(beta*e_charge);

	public double ni_semi = 1e16;					// Semiconductor equilibrium concentration
	public double W_semi = 4.7*eVtoJ;				// Semiconductor work function
	public double E_b_semi = 1.12*eVtoJ;			// Semiconductor band gap
	public double recomb_rate_semi = 1e7/ni_semi;

	public double ni_metal = 5e20;					// Metal charge carrier concentratio
	public double ni_metal_high = 2.5*ni_metal;
	public double ni_metal_low = 0.25*ni_metal;
	public double W_metal_default = 4.7*eVtoJ;		// Metal work function, same as semiconductor
	public double W_metal_high = W_semi + 0.3*eVtoJ;
	public double W_metal_low = W_semi - 0.3*eVtoJ;
	public double E_b_metal = 1.12*eVtoJ;			// Same as semiconductor
	public double recomb_rate_metal = 1e8/ni_semi;

	// Only reason activation energy needed is to smoothly public interpolate rate constant
	public double recomb_cross_section = 1e-10;
	public double arrhenius_prefactor = recomb_cross_section*Math.sqrt(8*(k*T)/(Math.PI*e_mass/2));
	public double E_a_semi = -Math.log(recomb_rate_semi/arrhenius_prefactor);
	public double E_a_metal = -Math.log(recomb_rate_metal/arrhenius_prefactor);

	public double n_default_doping_concentration = 5e19;
	public double p_default_doping_concentration = 5e19;
	public double n_light_doping_concentration = 1e19;
	public double p_light_doping_concentration = 1e19;
	public double n_heavy_doping_concentration = 2.5e20;
	public double p_heavy_doping_concentration = 2.5e20;

	public double dielectric_eps_r = 25.0;
	public double ferromagnet_mu_r = 250.0;
	public double staticcharge_density = 10.0;

	public double E_sat = 5e5;						// Electric field at which carrier velocity saturates

	public int junction_size = 3;					// Free energy smoothing distance: junction_size < 3 causes instability
	
	
	/* Dynamical simulation variables */

	public double[][] Ex;			// x component of E field
	public double[][] Ey;			// y component of E field
	public double[][] Hz;			// z component of B field
	public double[][] Hz_laplacian;

	public double[][] rho_n;		// Electrons
	public double[][] rho_p;		// Holes
	public double[][] rho_back;	// Background
	public double[][] rho_abs;		// Absorber charge
	public double[][] rho_free;	// Total free charges
	public double[][] mobility_factor;

	public double[][] Jx_n;		// Electron current
	public double[][] Jy_n;
	public double[][] Jx_p;		// Hole current
	public double[][] Jy_p;
	public double[][] Jx_abs;		// Absorber current
	public double[][] Jy_abs;
	public double[][] Jx_free;		// Total free current
	public double[][] Jy_free;

	/* Miscellaneous fields used for display */

	public double[][] Dx;				// Electric displacement field
	public double[][] Dy;
	public double[][] Bz;				// B field
	public double[][] Sx;				// Poynting vector
	public double[][] Sy;
	public double[][] u;				// EM energy density
	public double[][] phi;				// Electric scalar potential

	public double[][] G;				// Generation rate
	public double[][] R;				// Recombination rate
	public double[][] F_n;				// Total chemical potential of electrons (quasi-Fermi level minus the electrostatic potential)
	public double[][] F_p;				// Total chemical potential of holes
	public double[][] grad_E0x_n;		// Gradient of chemical energy
	public double[][] grad_E0y_n;
	public double[][] grad_E0x_p;
	public double[][] grad_E0y_p;
	public double[][] grad_Fx_n;		// Gradient of total chemical potential
	public double[][] grad_Fy_n;
	public double[][] grad_Fx_p;
	public double[][] grad_Fy_p;
	public double[][] Q;				// Heat dissipation rate
	public double[][] S;				// Free energy dissipation rate
	public double[][] F;				// Average total electrochemical potential

	public boolean updateMiscFields = false;
	public double[][] debug;

	
	/* Material parameters */

	public Material[][] materials;

	public double[][] F0_n;		// Standard chemical potential of electrons
	public double[][] F0_p;		// Standard chemical potential of holes
	public double[][] E0_n;		// Standard chemical energy of electrons
	public double[][] E0_p;		// Standard chemical energy of holes
	public double[][] K;			// Charge carrier equilibrium constant
	public double[][] E_a;			// Activation energy of recombination
	public double[][] r;			// Recombination rate constant
	public double[][] L;			// Carrier generation due to incoming light

	public double[][] cmfx_n;		// Chemical-motive force for electrons
	public double[][] cmfy_n;

	public double[][] cmfx_p;		// Chemical-motive force for holes
	public double[][] cmfy_p;

	public double[][] emfx;		// External electromotive force
	public double[][] emfy;
	public double[][] epsx;		// Dielectric constant
	public double[][] epsy;
	public double[][] mu_z;		// Relative permeability

	public int[][] conducting;		// Does the material have partially filled bands
	public int[][] conducting_x;
	public int[][] conducting_y;

	public int[][] ac_x;
	public int[][] ac_y;

	public double[][] absorptivity;
	public double[][] absorptivity_x;
	public double[][] absorptivity_y;

	
	/* Multigrid Poisson eq solver */

	public int log2_resolution;
	public double[][] MG_rho0;
	public double[][][] MG_rho;
	public double[][][] MG_epsx;
	public double[][][] MG_epsy;
	public double[][] MG_eps_avg;
	public double[][] MG_phi1;
	public double[][] MG_phi2;

	
	/* Pathfinding */

	public int[][] distance;
	public boolean[][] visited;
	public int[][] abs_depth;
	
	
	/* Probes */
	public List<VoltageProbe> voltageprobes = new CopyOnWriteArrayList<>();
	public List<CurrentProbe> currentprobes = new CopyOnWriteArrayList<>();
	public List<ChargeProbe> chargeprobes = new CopyOnWriteArrayList<>();
	public VoltageProbe ground = null;
	public String datafilename = "probedata.txt";
	public File datafile;
	public PrintWriter datastream;


	/* Performance profiling */
	
	Timer t4 = new Timer("Poisson constraint solver", 1, true);
	Timer t7 = new Timer("Poisson potential solver", 10, true);
	Timer t6 = new Timer("Iterate simulation", 40, true);
	Timer t8 = new Timer("Calc misc fields", 20, true);
	Timer t9 = new Timer("Debug", 20, false);
	Timer simFPStimer = new Timer("Simulation FPS", 10, true);

	public Simulation() {
		controls = new Controls(this);
		renderer = new Renderer(this);
		savemanager = new SaveManager(this);
		opts = new MainWindow();
		canvas = new RenderCanvas(this);
		canvas.setFocusable(true);

		bandplot = new BandPlot(); plots.add(bandplot);
		scalarplot = new ScalarPlot(); plots.add(scalarplot);
		carrierplot = new CarrierPlot(); plots.add(carrierplot);

		SemiSim.detect64Bit();

		setSize(default_resolution, default_width);
		resetFields(true);
		
		opts.add(canvas, BorderLayout.CENTER);
		canvas.addMouseListener(controls);
		canvas.addMouseMotionListener(controls);
		canvas.addMouseWheelListener(controls);
		canvas.addKeyListener(controls);
		opts.gui_reset.addActionListener(controls);
		opts.gui_view.addActionListener(controls);
		opts.gui_view_vec.addActionListener(controls);
		opts.gui_brush.addActionListener(controls);
		opts.gui_material.addActionListener(controls);
		opts.gui_adv_settings.addActionListener(controls);
		
		opts.menu_open.addActionListener(controls);
		opts.menu_save.addActionListener(controls);
		opts.menu_about.addActionListener(controls);
		opts.menu_help.addActionListener(controls);
		opts.menu_undo.addActionListener(controls);
		opts.menu_redo.addActionListener(controls);
		opts.menu_save.addActionListener(controls);
		opts.menu_cut.addActionListener(controls);
		opts.menu_copy.addActionListener(controls);
		opts.menu_paste.addActionListener(controls);
		opts.menu_editdesc.addActionListener(controls);
		opts.menu_new.addActionListener(controls);
		
		opts.gui_brush.addItemListener(controls);
		
		controls.brushbuttonmap = new HashMap<Brush, JRadioButtonMenuItem>();
		controls.buttongroup = new ButtonGroup();
		for (Brush b : Brush.values()) {
			if (b == Brush.DRAW || b == Brush.VOLTAGE || b == Brush.BANDS)
				opts.menu_tools.add(new JSeparator());
			
			JRadioButtonMenuItem button = new JRadioButtonMenuItem(b.name);
			controls.brushbuttonmap.put(b, button);
			button.addActionListener(controls);
			controls.buttongroup.add(button);
			opts.menu_tools.add(button);
		}
		controls.buttongroup.setSelected(controls.brushbuttonmap.get(Brush.INTERACT).getModel(), true);
		
		MenuBuilder.addDirectoryToMenu(opts.menu_examples, new File("examples"), savemanager.fileextension, (File f) -> savemanager.readfile(f));

		//opts.gui_material.removeItem(MaterialType.ABSORBER);
		opts.gui_view.removeItem(ScalarView.DEBUG);

		opts.gui_parameter1.setEnabled(true);
		opts.gui_parameter1.setVisible(true);
		opts.gui_parameter1_text.setEnabled(true);
		opts.gui_parameter1_text.setVisible(true);

		opts.pack();

		controls.addKeyBinds(canvas);
		controls.addKeyBinds(opts.panel);

		adv_opts = new AdvancedOptions();
		adv_opts.btn_apply.addActionListener(controls);
		adv_opts.btn_cancel.addActionListener(controls);
		adv_opts.setVisible(false);
		
		InputMap im = (InputMap)UIManager.get("Button.focusInputMap");
		im.put(KeyStroke.getKeyStroke("pressed SPACE"), "none");
		im.put(KeyStroke.getKeyStroke("released SPACE"), "none");


		datafile = new File(datafilename);

		try {
			datastream = new PrintWriter(new FileOutputStream(datafile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		rwLock.readLock().lock();
        try {
    		try {
    			iteration_multiplier = opts.gui_simspeed_2.getValue();

    			if (controls.clear) {
    				SwingUtilities.invokeLater(() -> {
    					resetFields(false);
    					multigridSolve(true, false);
    				});
    				time = 0.0;
    				controls.clear = false;
    			}

    			if (controls.reset) {
    				SwingUtilities.invokeLater(() -> {
        				int result = JOptionPane.showConfirmDialog(opts, "Do you wish to reset the entire simulation?", "Message", JOptionPane.YES_NO_OPTION);
        				if (result == JOptionPane.OK_OPTION)
        				{
        					resetFields(true);
        					time = 0.0;
        				}
    				});
    				controls.reset = false;
    			}

    			if (opts.gui_simspeed.getValue() != lastsimspeed) {
    				lastsimspeed = opts.gui_simspeed.getValue();
    				dt = dt_maximum*(lastsimspeed/20.0);
    			}

    			
    			controls.handleUndoRedo();
    			
    			controls.handleMouseInput();

    			if (!opts.gui_paused.isSelected() || controls.advanceframe) {
    				for (int i = 0; i < iteration_multiplier ; i++) {
    					start_barrier.await();
    					stop_barrier.await();
    				}
    				calcMiscFields(false);
    				frame++;
    			} else if (updateMiscFields) {
    				calcMiscFields(false);
    			}
    			
    			if (numerical_overflow)
    				opts.gui_paused.setSelected(true);


    			if (controls.save) {
    				savemanager.writeFile();
    				controls.save = false;
    			}

    			if (controls.load) {
    				savemanager.readFile();
    				controls.load = false;
    			}

    			simFPStimer.stop();
    			simFPStimer.start();

    		} catch (Exception e) {
    			SemiSim.displayErrorMessage(e);
    		}
        } finally {
            rwLock.readLock().unlock();
        }

        SemiSim.instance.threadPool.schedule(this, nextDelay(renderer.frameduration), TimeUnit.MILLISECONDS);
	}

	TimerTask potentialSolver = new PeriodicTask() {
		@Override
		public void run() {
	        rwLock.readLock().lock();
	        try {
				if (!opts.gui_paused.isSelected() || controls.advanceframe || updateMiscFields) {
					multigridSolve(false, true);
					//calcMiscFields(true);
				}
	        } catch( Exception e) {
	        	SemiSim.displayErrorMessage(e);
	        } finally {
	            rwLock.readLock().unlock();
	        }
	        SemiSim.instance.threadPool.schedule(this, nextDelay(renderer.frameduration), TimeUnit.MILLISECONDS);
		}
	};
	
	public void updateConstants() {
		c = 1/Math.sqrt(eps0*mu0);
		beta = 1/(k*T);
		D_electron = mu_electron/(beta*e_charge);
		D_hole = mu_hole/(beta*e_charge);
		arrhenius_prefactor = recomb_cross_section*Math.sqrt(8*(k*T)/(Math.PI*e_mass/2));
		E_a_semi = -Math.log(recomb_rate_semi/arrhenius_prefactor);
		E_a_metal = -Math.log(recomb_rate_metal/arrhenius_prefactor);
	}

	public boolean setSize(int resolution, double width) {
		
		this.width = width;
		ds = width/resolution;
		width = nx*ds;
		dt_maximum = 0.9*ds/(Math.sqrt(2)*c);
		Hz_dissipation = 0.01*ds*ds/dt_maximum;

		if (resolution == this.resolution) return false;

		// Simulation and graphics threads should not be active when variables are initialized
		rwLock.writeLock().lock();
		try {
			if(resolution < 4) {
				JOptionPane.showMessageDialog(opts, "Resolution too small.", "Error", JOptionPane.OK_OPTION);
				resolution = 4;
			}
			
			log2_resolution = (int) Math.round(Math.log(resolution)/Math.log(2));
			if(1 << log2_resolution != resolution) {
				JOptionPane.showMessageDialog(opts, "Resolution must be a power of 2.", "Error", JOptionPane.OK_OPTION);
				resolution = 1 << log2_resolution;
			}
			
			this.resolution = resolution;
			nx = resolution;
			ny = resolution;

			Ex = new double[nx][ny];
			Ey = new double[nx][ny];
			Bz = new double[nx][ny];
			Hz_laplacian = new double[nx][ny];
			rho_abs = new double[nx][ny];
			rho_n = new double[nx][ny];
			rho_p = new double[nx][ny];
			rho_back = new double[nx][ny];
			rho_free = new double[nx][ny];
			mobility_factor = new double[nx][ny];

			Jx_abs = new double[nx][ny];
			Jy_abs = new double[nx][ny];
			Jx_n = new double[nx][ny];
			Jy_n = new double[nx][ny];
			Jx_p = new double[nx][ny];
			Jy_p = new double[nx][ny];
			Jx_free = new double[nx][ny];
			Jy_free = new double[nx][ny];

			materials = new Material[nx][ny];
			K = new double[nx][ny];
			F0_n = new double[nx][ny];
			F0_p = new double[nx][ny];
			F_n = new double[nx][ny];
			F_p = new double[nx][ny];
			E0_n = new double[nx][ny];
			E0_p = new double[nx][ny];
			E_a = new double[nx][ny];
			r = new double[nx][ny];
			L = new double[nx][ny];
			cmfx_n = new double[nx][ny];
			cmfy_n = new double[nx][ny];
			cmfx_p = new double[nx][ny];
			cmfy_p = new double[nx][ny];
			conducting = new int[nx][ny];
			conducting_x = new int[nx][ny];
			conducting_y = new int[nx][ny];
			ac_x = new int[nx][ny];
			ac_y = new int[nx][ny];
			absorptivity = new double[nx][ny];
			absorptivity_x = new double[nx][ny];
			absorptivity_y = new double[nx][ny];
			emfx = new double[nx][ny];
			emfy = new double[nx][ny];
			epsx = new double[nx][ny];
			epsy = new double[nx][ny];
			mu_z = new double[nx][ny];

			Dx = new double[nx][ny];
			Dy = new double[nx][ny];
			Hz = new double[nx][ny];
			Sx = new double[nx][ny];
			Sy = new double[nx][ny];
			u = new double[nx][ny];
			phi = new double[nx][ny];

			G = new double[nx][ny];
			R = new double[nx][ny];
			F_n = new double[nx][ny];
			F_p = new double[nx][ny];
			grad_E0x_n = new double[nx][ny];
			grad_E0y_n = new double[nx][ny];
			grad_E0x_p = new double[nx][ny];
			grad_E0y_p = new double[nx][ny];
			grad_Fx_n = new double[nx][ny];
			grad_Fy_n = new double[nx][ny];
			grad_Fx_p = new double[nx][ny];
			grad_Fy_p = new double[nx][ny];
			Q = new double[nx][ny];
			S = new double[nx][ny];
			F = new double[nx][ny];
			debug = new double[nx][ny];

			MG_rho0 = new double[nx][ny];
			MG_rho = new double[log2_resolution+1][nx][ny];
			MG_epsx = new double[log2_resolution+1][nx][ny];
			MG_epsy = new double[log2_resolution+1][nx][ny];
			MG_eps_avg = new double[nx][ny];
			MG_phi1 = new double[nx][ny];
			MG_phi2 = new double[nx][ny];

			distance = new int[nx][ny];
			visited = new boolean[nx][ny];
			abs_depth = new int[nx][ny];

			opts.setVisible(true);

			controls.setResolution(resolution);
			renderer.setResolution(resolution);

			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					materials[i][j] = new Material();
					controls.selection[i][j] = new ClipboardMaterial();
					controls.clipboard[i][j] = new ClipboardMaterial();
				}
			}

			resetFields(true);
		} finally {
			rwLock.writeLock().unlock();
		}
		
		return true;
	}

	public void resetFields(boolean resetall) {

		rwLock.writeLock().lock();
		try {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{

					if (resetall) {

						materials[i][j].erase();
						controls.selection[i][j].erase();
						controls.clipboard[i][j].erase();

						F0_n[i][j] = 0;
						F0_p[i][j] = 0;
						F_n[i][j] = 0;
						F_p[i][j] = 0;
						E0_n[i][j] = 0;
						E0_p[i][j] = 0;
						K[i][j] = 0;
						E_a[i][j] = 0;
						r[i][j] = 0;
						L[i][j] = 0;

						cmfx_n[i][j] = 0;
						cmfy_n[i][j] = 0;
						cmfx_p[i][j] = 0;
						cmfy_p[i][j] = 0;

						emfx[i][j] = 0.0;
						emfy[i][j] = 0.0;

						epsx[i][j] = eps0;
						epsy[i][j] = eps0;
						mu_z[i][j] = mu0;

						conducting[i][j] = 0;
						conducting_x[i][j] = 0;
						conducting_y[i][j] = 0;

						ac_x[i][j] = 0;
						ac_y[i][j] = 0;

						absorptivity[i][j] = 0;
						absorptivity_x[i][j] = 0;
						absorptivity_y[i][j] = 0;
					}

					Ex [i][j] = 0.0;
					Ey [i][j] = 0.0;
					Bz [i][j] = 0.0;
					Hz_laplacian [i][j] = 0.0;

					rho_abs[i][j] = 0.0;
					rho_n[i][j] = 0.0;
					rho_p[i][j] = 0.0;
					rho_back[i][j] = 0.0;
					rho_free[i][j] = 0.0;
					mobility_factor[i][j] = 0.0;

					Jx_abs[i][j] = 0.0;
					Jy_abs[i][j] = 0.0;
					Jx_n[i][j] = 0.0;
					Jy_n[i][j] = 0.0;
					Jx_p[i][j] = 0.0;
					Jy_p[i][j] = 0.0;
					Jx_free[i][j] = 0.0;
					Jy_free[i][j] = 0.0;

					MG_rho0 [i][j] = 0.0;
					MG_eps_avg [i][j] = 0.0;
					MG_phi1 [i][j] = 0.0;
					MG_phi2 [i][j] = 0.0;
					for (int n = 0; n < log2_resolution+1; n++) {
						MG_rho[n][i][j] = 0;
						MG_epsx[n][i][j] = 0;
						MG_epsy[n][i][j] = 0;
					}

					distance[i][j] = Integer.MAX_VALUE;
					visited[i][j] = false;

					Dx[i][j] = 0.0;
					Dy[i][j] = 0.0;
					Hz[i][j] = 0.0;
					Sx[i][j] = 0.0;
					Sy[i][j] = 0.0;
					u[i][j] = 0.0;
					phi[i][j] = 0.0;

					G[i][j] = 0.0;
					R[i][j] = 0.0;
					F_n[i][j] = 0.0;
					F_p[i][j] = 0.0;
					grad_E0x_n[i][j] = 0.0;
					grad_E0y_n[i][j] = 0.0;
					grad_E0x_p[i][j] = 0.0;
					grad_E0y_p[i][j] = 0.0;
					grad_Fx_n[i][j] = 0.0;
					grad_Fy_n[i][j] = 0.0;
					grad_Fx_p[i][j] = 0.0;
					grad_Fy_p[i][j] = 0.0;
					Q[i][j] = 0.0;
					S[i][j] = 0.0;
					F[i][j] = 0.0;
					debug[i][j] = 0.0;

					controls.selected[i][j] = false;
					controls.selected_EMF[i][j] = false;
				}
			}

			numerical_overflow = false;
			
			if (resetall) {
				voltageprobes.clear();
				currentprobes.clear();
				chargeprobes.clear();
				ground = null;

				for (Plot p: plots) {
					p.frame.setVisible(false);
				}

				opts.setTitle("Brandon's semiconductor simulator");

				controls.resetUndoHistory();
			}

			renderer.reset();

			initializeAllMaterials();
			updateAllMaterials(true);
			controls.captureState();
			checkCFL();
		}
		finally {
			rwLock.writeLock().unlock();
		}

	}

	public void constructBoundary() {
		absorber_width = (int)Math.ceil(0.045*nx);
		absorbing_coeff = 50*c/(ds*nx);
		double max_stretch = 10;
		double coeff = Math.log(max_stretch);

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				double dx = (double)Math.max(0, absorber_width-Math.min(i, nx-1-i))/absorber_width;
				double dy = (double)Math.max(0, absorber_width-Math.min(j, ny-1-j))/absorber_width;
				double depth = Math.sqrt(dx*dx+dy*dy);
				if (depth > 0) {
					if ((BoundaryCondition)opts.gui_bc.getSelectedItem() == BoundaryCondition.DISSIPATIVE) {
						if (materials[i][j].type == MaterialType.VACUUM) {
							initializeMaterial(i, j, MaterialType.ABSORBER);
							materials[i][j].auto_placed = true;
						}
					} else if (materials[i][j].type == MaterialType.ABSORBER) {
						if (materials[i][j].auto_placed) {
							eraseMaterial(i, j);
						}
					}
				}
			}
		}

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				abs_depth[i][j] = 0;
			}
		}
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (materials[i][j].type != MaterialType.ABSORBER) {
					set(i+1, j, 1);
					set(i-1, j, 1);
					set(i, j+1, 1);
					set(i, j-1, 1);
				}
			}
		}


		for (int d = 1; d < nx; d++) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					if (abs_depth[i][j] == d) {
						set(i+1, j, d+1);
						set(i-1, j, d+1);
						set(i, j+1, d+1);
						set(i, j-1, d+1);
					}
				}
			}
		}
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (materials[i][j].type == MaterialType.ABSORBER) {
					double depth = (double)abs_depth[i][j]/absorber_width;
					double stretchfactor = Math.exp(coeff*depth);
					materials[i][j].eps_r = stretchfactor;
					materials[i][j].mu_r = stretchfactor;
					materials[i][j].absorptivity = 1;
				}
			}
		}
	}
	
	public void set(int i, int j, int d) {
		if (i < 0 || j < 0 || i >= nx || j >= ny) return;
		if (materials[i][j].type == MaterialType.ABSORBER && abs_depth[i][j] == 0)
			abs_depth[i][j] = d;
	}

	class SimulationThread extends Thread {

		int i_min;
		int i_max;
		int n_thread;

		public SimulationThread(int n, int n_threads, int nx) {
			n_thread = n;
			System.out.println("Thread " + n_thread + ": " + i_min + " < i <= " + i_max);
		}

		@Override
		public void run() {
			try {
				while (true) {
					start_barrier.await();

					i_min = (n_thread*nx)/SemiSim.n_threads;
					i_max = (n_thread+1)*nx/SemiSim.n_threads-1;

					if (n_thread == 0) {
						t6.start();

						stepnumber++;
					}

					/* Interior B field & charge */
					for (int i = 1; i < nx-2; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-2; j++)
							{
								Hz_laplacian[i][j] = (Hz[i+1][j]+Hz[i][j+1]+Hz[i-1][j]+Hz[i][j-1]-4*Hz[i][j])/(ds*ds);
							}
						}
					}

					mid_barrier.await();

					for (int i = 0; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 0; j < ny-1; j++)
							{
								double sigma = absorptivity[i][j]*mu_z[i][j]*absorbing_coeff;

								Hz[i][j] = (Hz[i][j]*(1-0.5*dt*sigma/mu_z[i][j])
										+ (-(Ey[i+1][j] - Ey[i][j]) + (Ex[i][j+1] - Ex[i][j]))*dt/(ds*mu_z[i][j])
										+ Hz_dissipation*dt*Hz_laplacian[i][j])
										/(1+0.5*dt*sigma/mu_z[i][j]);

								// Includes unphysical magnetic monopole current term to make absorption very effective
							}
						}
					}

					/* Update charge carriers */
					for (int i = 1; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-1; j++)
							{
								double generation_rate = conducting[i][j]*(r[i][j]*(K[i][j] - rho_n[i][j]*rho_p[i][j]/(q_n*q_p)) + L[i][j]);

								rho_n[i][j] = rho_n[i][j] - (Jx_n[i][j]-Jx_n[i-1][j] + Jy_n[i][j]-Jy_n[i][j-1])*dt/ds + dt*q_n*generation_rate;
								rho_p[i][j] = rho_p[i][j] - (Jx_p[i][j]-Jx_p[i-1][j] + Jy_p[i][j]-Jy_p[i][j-1])*dt/ds + dt*q_p*generation_rate;
								rho_abs[i][j] = rho_abs[i][j] - (Jx_abs[i][j]-Jx_abs[i-1][j] + Jy_abs[i][j]-Jy_abs[i][j-1])*dt/ds;

								rho_free[i][j] = rho_abs[i][j]+rho_n[i][j]+rho_p[i][j]+rho_back[i][j];

								double E_avg = Math.sqrt(0.5*(Ex[i][j]*Ex[i][j] + Ex[i-1][j]*Ex[i-1][j] + Ey[i][j]*Ey[i][j] + Ey[i][j-1]*Ey[i][j-1]));
								mobility_factor[i][j] = Math.min(1, E_sat/E_avg);
							}
						}
					}

					mid_barrier.await();

					/* Update E field and currents */
					for (int i = 0; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-1; j++)
							{
								double ex_prev = Ex[i][j];

								double mf = Math.min(mobility_factor[i+1][j], mobility_factor[i][j]);
								//double mobility_factor = Math.min(1, E_sat/Math.abs(Ex[i][j]));

								double sigma_n = conducting_x[i][j]*mf*mu_electron*Utils.logmean(-rho_n[i+1][j],-rho_n[i][j]);
								double sigma_p = conducting_x[i][j]*mf*mu_hole*Utils.logmean(rho_p[i+1][j], rho_p[i][j]);
								
								double emf_phase = ac_x[i][j]*AC_amplitude+(1-ac_x[i][j]);

								Jx_abs[i][j] = 0;

								Jx_n[i][j] = conducting_x[i][j]*(-mf*D_electron*(rho_n[i+1][j] - rho_n[i][j])/ds
										+ sigma_n*(emf_phase*emfx[i][j] + cmfx_n[i][j]/q_n));

								Jx_p[i][j] = conducting_x[i][j]*(-mf*D_hole*(rho_p[i+1][j] - rho_p[i][j])/ds
										+ sigma_p*(emf_phase*emfx[i][j] + cmfx_p[i][j]/q_p));

								double sigma = sigma_n + sigma_p + absorptivity_x[i][j]*epsx[i][j]*absorbing_coeff;

								double jx = Jx_abs[i][j] + Jx_n[i][j] + Jx_p[i][j];

								Ex[i][j] = (Ex[i][j]*(1-0.5*dt*sigma/epsx[i][j]) + ((Hz[i][j]-Hz[i][j-1])/ds - jx)*dt/epsx[i][j])
										/(1+0.5*dt*sigma/epsx[i][j]);

								Jx_abs[i][j] += 0.5*(absorptivity_x[i][j]*epsx[i][j]*absorbing_coeff)*(ex_prev + Ex[i][j]);
								Jx_n[i][j] += 0.5*sigma_n*(ex_prev + Ex[i][j]);
								Jx_p[i][j] += 0.5*sigma_p*(ex_prev + Ex[i][j]);
							}
						}
					}

					for (int i = 1; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 0; j < ny-1; j++)
							{
								double ey_prev = Ey[i][j];

								double mf = Math.min(mobility_factor[i][j+1], mobility_factor[i][j]);
								//double mobility_factor = Math.min(1, E_sat/Math.abs(Ey[i][j]));

								double sigma_n = conducting_y[i][j]*mf*mu_electron*Utils.logmean(-rho_n[i][j+1],-rho_n[i][j]);
								double sigma_p = conducting_y[i][j]*mf*mu_hole*Utils.logmean(rho_p[i][j+1], rho_p[i][j]);

								double emf_phase = ac_y[i][j]*AC_amplitude+(1-ac_y[i][j]);

								Jy_abs[i][j] = 0;

								Jy_n[i][j] = conducting_y[i][j]*(-mf*D_electron*(rho_n[i][j+1] - rho_n[i][j])/ds
										+ sigma_n*(emf_phase*emfy[i][j] + cmfy_n[i][j]/q_n));

								Jy_p[i][j] = conducting_y[i][j]*(-mf*D_hole*(rho_p[i][j+1] - rho_p[i][j])/ds
										+ sigma_p*(emf_phase*emfy[i][j] + cmfy_p[i][j]/q_p));

								double sigma = sigma_n + sigma_p + absorptivity_y[i][j]*epsy[i][j]*absorbing_coeff;

								double jy = Jy_abs[i][j] + Jy_n[i][j] + Jy_p[i][j];

								Ey[i][j] = (Ey[i][j]*(1-0.5*dt*sigma/epsy[i][j]) + (-(Hz[i][j]-Hz[i-1][j])/ds - jy)*dt/epsy[i][j])
										/(1+0.5*dt*sigma/epsy[i][j]);

								Jy_abs[i][j] += 0.5*(absorptivity_y[i][j]*epsy[i][j]*absorbing_coeff)*(ey_prev + Ey[i][j]);
								Jy_n[i][j] += 0.5*sigma_n*(ey_prev + Ey[i][j]);
								Jy_p[i][j] += 0.5*sigma_p*(ey_prev + Ey[i][j]);
							}
						}
					}

					if (n_thread == 0) {
						/* Apply boundary condition */
						for (int j = 0; j < ny-1; j++)
						{
							Ey[0][j] = 0;
							Ey[nx-1][j] = 0;
						}

						for (int i = 0; i < nx-1; i++)
						{
							Ex[i][0] = 0;
							Ex[i][ny-1] = 0;
						}

						t6.stop();
					}

					stop_barrier.await();
					
					if (n_thread == 0) {
						/* Enforce Gauss law constraint */
						if (stepnumber%500 == 0) {
							multigridSolve(true, false);
						}

						time += dt;
						AC_phase += 2*Math.PI*AC_freq*dt;
						AC_amplitude = Math.cos(AC_phase);

						controls.advanceframe = false;
					}
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void calcMiscFields(boolean updatePhi) {
		ScalarView view_scalar = (ScalarView) opts.gui_view.getSelectedItem();

		if (updatePhi)
			multigridSolve(false, true);

		t8.start();

		sign_violation = false;
		double rho_n_max_tmp = 0;
		double rho_p_max_tmp = 0;
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				// Prevent errors when taking log of 0
				if (conducting[i][j] == 0) {
					F_n[i][j] = Double.NaN;
					F_p[i][j] = Double.NaN;
					F[i][j] = Double.NaN;
				} else {
					F_n[i][j] = F0_n[i][j] + k*T*Math.log(rho_n[i][j]/q_n);
					F_p[i][j] = F0_p[i][j] + k*T*Math.log(rho_p[i][j]/q_p);
					// Think about why we should take average


					//double probe_n = -e_charge*ni_metal;
					//double probe_p = e_charge*ni_metal;
					//double sigma_n = mu_electron*logmean(-rho_n[i][j], -probe_n);
					//double sigma_p = mu_hole*logmean(rho_p[i][j], probe_p);
					double sigma_n = mu_electron*(-rho_n[i][j]);
					double sigma_p = mu_hole*rho_p[i][j];

					F[i][j] = (sigma_n*(F_n[i][j]/q_n+phi[i][j]) + sigma_p*(F_p[i][j]/q_p+phi[i][j]))
							/(sigma_n + sigma_p) - W_semi/eVtoJ;
				}

				G[i][j] = conducting[i][j]*(r[i][j]*K[i][j] + L[i][j]);
				R[i][j] = conducting[i][j]*r[i][j]*rho_n[i][j]*rho_p[i][j]/(q_n*q_p);

				if (rho_n[i][j] > error_detection_threshold || rho_p[i][j] < -error_detection_threshold)
					sign_violation = true;

				if (-rho_n[i][j] > rho_n_max_tmp)
					rho_n_max_tmp = -rho_n[i][j];

				if (rho_p[i][j] > rho_p_max_tmp)
					rho_p_max_tmp = rho_p[i][j];
			}
		}
		renderer.rho_n_max = rho_n_max_tmp;
		renderer.rho_p_max = rho_p_max_tmp;

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				Jx_free[i][j] = Jx_n[i][j] + Jx_p[i][j];
				Dx[i][j] = epsx[i][j]*Ex[i][j];
				Sy[i][j] = -Ex[i][j]*0.5*(Hz[i][j] + Hz[i][j-1]);
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				Jy_free[i][j] = Jy_n[i][j] + Jy_p[i][j];
				Dy[i][j] = epsy[i][j]*Ey[i][j];
				Sx[i][j] = Ey[i][j]*0.5*(Hz[i][j] + Hz[i-1][j]);
			}
		}

		for (int i = 1; i < nx-2; i++)
		{
			for (int j = 1; j < ny-2; j++)
			{
				Bz[i][j] = Hz[i][j]*mu_z[i][j];
				u[i][j] = 0.25*(Ex[i][j]*Dx[i][j] + Ex[i][j+1]*Dx[i][j+1])
						+ 0.25*(Ey[i][j]*Dy[i][j] + Ey[i+1][j]*Dy[i+1][j])
						+ 0.5*Hz[i][j]*Bz[i][j];
				
				if (!Double.isFinite(Hz[i][j]))
					numerical_overflow = true;
			}
		}

		if (view_scalar == ScalarView.HEAT) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					grad_E0x_n[i][j] = -conducting_x[i][j]*(E0_n[i+1][j] - E0_n[i][j])/ds;
					grad_E0x_p[i][j] = -conducting_x[i][j]*(E0_p[i+1][j] - E0_p[i][j])/ds;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					grad_E0y_n[i][j] = -conducting_y[i][j]*(E0_n[i][j+1] - E0_n[i][j])/ds;
					grad_E0y_p[i][j] = -conducting_y[i][j]*(E0_p[i][j+1] - E0_p[i][j])/ds;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					double n_contrib = -(G[i][j]-R[i][j])*E0_n[i][j];
					double jn_contrib = 0.5*(grad_E0x_n[i-1][j]*Jx_n[i-1][j]+grad_E0x_n[i][j]*Jx_n[i][j]
							+grad_E0y_n[i][j-1]*Jy_n[i][j-1]+grad_E0y_n[i][j]*Jy_n[i][j])/q_n;

					double p_contrib = -(G[i][j]-R[i][j])*E0_p[i][j];
					double jp_contrib = 0.5*(grad_E0x_p[i-1][j]*Jx_p[i-1][j]+grad_E0x_p[i][j]*Jx_p[i][j]
							+grad_E0y_p[i][j-1]*Jy_p[i][j-1]+grad_E0y_p[i][j]*Jy_p[i][j])/q_p;

					double ohm_contrib = 0.5*(Ex[i-1][j]*Jx_free[i-1][j]+Ex[i][j]*Jx_free[i][j]+Ey[i][j-1]*Jy_free[i][j-1]+Ey[i][j]*Jy_free[i][j]);

					Q[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
				}
			}
		}

		if (view_scalar == ScalarView.ENTROPY) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					grad_Fx_n[i][j] = conducting_x[i][j] == 1? -(F_n[i+1][j] - F_n[i][j])/ds : 0;
					grad_Fx_p[i][j] = conducting_x[i][j] == 1? -(F_p[i+1][j] - F_p[i][j])/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					grad_Fy_n[i][j] = conducting_y[i][j] == 1? -(F_n[i][j+1] - F_n[i][j])/ds : 0;
					grad_Fy_p[i][j] = conducting_y[i][j] == 1? -(F_p[i][j+1] - F_p[i][j])/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					double n_contrib = -(G[i][j]-R[i][j])*F_n[i][j];
					double jn_contrib = 0.5*(grad_Fx_n[i-1][j]*Jx_n[i-1][j]+grad_Fx_n[i][j]*Jx_n[i][j]
							+grad_Fy_n[i][j-1]*Jy_n[i][j-1]+grad_Fy_n[i][j]*Jy_n[i][j])/q_n;

					double p_contrib = -(G[i][j]-R[i][j])*F_p[i][j];
					double jp_contrib = 0.5*(grad_Fx_p[i-1][j]*Jx_p[i-1][j]+grad_Fx_p[i][j]*Jx_p[i][j]
							+grad_Fy_p[i][j-1]*Jy_p[i][j-1]+grad_Fy_p[i][j]*Jy_p[i][j])/q_p;

					double ohm_contrib = 0.5*(Ex[i-1][j]*Jx_free[i-1][j]+Ex[i][j]*Jx_free[i][j]+Ey[i][j-1]*Jy_free[i][j-1]+Ey[i][j]*Jy_free[i][j]);

					S[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
				}
			}
		}

		if (ground != null)
			ground.potential = F[ground.x][ground.y];

		for (VoltageProbe p: voltageprobes) {
			p.potential = F[p.x][p.y];
		}

		for (CurrentProbe p: currentprobes) {
			p.current = calcCurrent(p.x1, p.y1, p.x2, p.y2);
		}
		
		for (ChargeProbe p: chargeprobes) {
			p.charge = calcCharge(p);
		}

		if (controls.logdata) {
			String data = "";
			data += ("t = " + Utils.getSI(time, "s"));

			for (VoltageProbe p: voltageprobes) {
				if (ground != null)
					data += (", V(" + Utils.getSI(p.x*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y+1)*ds, "m") + ") = " + Utils.getSI(p.potential-ground.potential, "V"));
				else
					data += (", V(" + Utils.getSI(p.x*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y+1)*ds, "m") + ") = " + Utils.getSI(p.potential, "V"));
			}

			for (CurrentProbe p: currentprobes) {
				data += (", I(" + Utils.getSI(p.x1*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y1+1)*ds, "m") + " - " + Utils.getSI(p.x2*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y2+1)*ds, "m") + ") = " + Utils.getSI(p.current*depth, "A"));
			}
			for (ChargeProbe p: chargeprobes) {
				data += (", Q(" + Utils.getSI(p.x1*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y1+1)*ds, "m") + " - " + Utils.getSI(p.x2*ds, "m") + ", " + Utils.getSI(ds*ny-(p.y2+1)*ds, "m") + ") = " + Utils.getSI(p.charge*depth, "C"));
			}

			data += "\n";

			datastream.print(data);
			datastream.flush();

			renderer.probetexttimer = 30;
			controls.logdata = false;
		}

		updateMiscFields = false;

		t8.stop();
	}
	
	public void computeChemicalForces()
	{

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 1) {
					computeFreeEnergy(i, j, materials[i][j].ni, materials[i][j].W, materials[i][j].Eb);
				} else {
					F0_n[i][j] = 0;
					F0_p[i][j] = 0;
					E0_n[i][j] = 0;
					E0_p[i][j] = 0;
				}

				E_a[i][j] = materials[i][j].Ea;
			}
		}

		// Smoothing out free energy in space makes simulation more stable
		double[][] F_n_tmp = copyArray(F0_n);
		double[][] F_p_tmp = copyArray(F0_p);
		double[][] E_n_tmp = copyArray(E0_n);
		double[][] E_p_tmp = copyArray(E0_p);
		double[][] Ea_tmp = copyArray(E_a);

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 1) {
					double F_n_sum = 0;
					double F_p_sum = 0;
					double E_n_sum = 0;
					double E_p_sum = 0;
					double Ea_sum = 0;
					double neighbors = 0;

					markNeighborhood(i, j, junction_size);

					for (int di = -junction_size; di <= junction_size; di++) {
						for (int dj = -junction_size; dj <= junction_size; dj++) {
							if (i+di >= 0 && j+dj >= 0 && i+di < nx && j+dj < ny && conducting[i+di][j+dj] == 1 && visited[i+di][j+dj]) {
								F_p_sum += F_p_tmp[i+di][j+dj];
								F_n_sum += F_n_tmp[i+di][j+dj];
								E_p_sum += E_p_tmp[i+di][j+dj];
								E_n_sum += E_n_tmp[i+di][j+dj];
								Ea_sum += Ea_tmp[i+di][j+dj];
								neighbors += 1;
								distance[i+di][j+dj] = Integer.MAX_VALUE;
								visited[i+di][j+dj] = false;
							}
						}
					}
					F0_n[i][j] = F_n_sum/neighbors;
					F0_p[i][j] = F_p_sum/neighbors;
					E0_n[i][j] = E_n_sum/neighbors;
					E0_p[i][j] = E_p_sum/neighbors;
					E_a[i][j] = Ea_sum/neighbors;
					K[i][j] = Math.exp(-beta*(F0_n[i][j]+F0_p[i][j]));
					r[i][j] = arrhenius_prefactor*Math.exp(-E_a[i][j]);
				} else {
					E_a[i][j] = 0;
					K[i][j] = 0;
					r[i][j] = 0;
				}
			}
		}


		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				cmfx_n[i][j] = -conducting_x[i][j]*(F0_n[i+1][j] - F0_n[i][j])/ds;
				cmfx_p[i][j] = -conducting_x[i][j]*(F0_p[i+1][j] - F0_p[i][j])/ds;
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				cmfy_n[i][j] = -conducting_y[i][j]*(F0_n[i][j+1] - F0_n[i][j])/ds;
				cmfy_p[i][j] = -conducting_y[i][j]*(F0_p[i][j+1] - F0_p[i][j])/ds;
			}
		}
	}

	public void updateAllMaterials(boolean updateRho) {
		constructBoundary();
		
		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				mu_z[i][j] = 0.25*mu0*(materials[i][j].mu_r+materials[i+1][j].mu_r+materials[i][j+1].mu_r+materials[i+1][j+1].mu_r);
				absorptivity[i][j] = 0.25*(materials[i][j].absorptivity+materials[i+1][j].absorptivity+materials[i][j+1].absorptivity+materials[i+1][j+1].absorptivity);
			}
		}

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				rho_back[i][j] = materials[i][j].rho_back;
				conducting[i][j] = materials[i][j].conducting*materials[i][j].activated;
			}
		}
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				conducting_x[i][j] = Math.min(materials[i+1][j].conducting*materials[i+1][j].activated, materials[i][j].conducting*materials[i][j].activated);
				emfx[i][j] = 0.5*(materials[i+1][j].emf*Math.cos(materials[i+1][j].emf_direction)*materials[i+1][j].activated
						+ materials[i][j].emf*Math.cos(materials[i][j].emf_direction)*materials[i][j].activated);
				epsx[i][j] = eps0*0.5*(materials[i+1][j].eps_r + materials[i][j].eps_r);
				absorptivity_x[i][j] = Math.min(materials[i+1][j].absorptivity, materials[i][j].absorptivity);
				ac_x[i][j] = (materials[i][j].type == MaterialType.AC_EMF || materials[i+1][j].type == MaterialType.AC_EMF) ? 1:0;
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				conducting_y[i][j] = Math.min(materials[i][j+1].conducting*materials[i][j+1].activated, materials[i][j].conducting*materials[i][j].activated);
				emfy[i][j] = 0.5*(materials[i][j+1].emf*Math.sin(materials[i][j+1].emf_direction)*materials[i][j+1].activated
						+ materials[i][j].emf*Math.sin(materials[i][j].emf_direction)*materials[i][j].activated);
				epsy[i][j] = eps0*0.5*(materials[i][j+1].eps_r + materials[i][j].eps_r);
				absorptivity_y[i][j] = Math.min(materials[i][j+1].absorptivity, materials[i][j].absorptivity);
				ac_y[i][j] = (materials[i][j].type == MaterialType.AC_EMF || materials[i][j+1].type == MaterialType.AC_EMF) ? 1:0;
			}
		}

		computeChemicalForces();

		AC_source_exists = false;
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (K[i][j] > 0 || materials[i][j].type == MaterialType.SWITCH) {
					if (updateRho || (rho_n[i][j] == 0 && rho_p[i][j] == 0)) {
						rho_n[i][j] = calcEquilibriumElectronCharge(rho_back[i][j], K[i][j]);
						rho_p[i][j] = calcEquilibriumHoleCharge(rho_back[i][j], K[i][j]);
					}
				} else {
					rho_n[i][j] = 0;
					rho_p[i][j] = 0;
				}

				if (materials[i][j].absorptivity == 0) {
					rho_abs[i][j] = 0;
				}

				rho_free[i][j] = rho_abs[i][j]+rho_n[i][j]+rho_p[i][j]+rho_back[i][j];
				
				if (materials[i][j].type == MaterialType.AC_EMF)
					AC_source_exists = true;
			}
		}

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				Jx_n[i][j] *= conducting_x[i][j];
				Jx_p[i][j] *= conducting_x[i][j];
				Jx_abs[i][j] *= (absorptivity_x[i][j] > 0)? 1:0;
			}
		}

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				Jy_n[i][j] *= conducting_y[i][j];
				Jy_p[i][j] *= conducting_y[i][j];
				Jy_abs[i][j] *= (absorptivity_y[i][j] > 0)? 1:0;
			}
		}

		prescaleDielectric();
	}


	public void updateJustEMFs() {
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				emfx[i][j] = 0.5*(materials[i+1][j].emf*Math.cos(materials[i+1][j].emf_direction)*materials[i+1][j].activated
						+ materials[i][j].emf*Math.cos(materials[i][j].emf_direction)*materials[i][j].activated);
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				emfy[i][j] = 0.5*(materials[i][j+1].emf*Math.sin(materials[i][j+1].emf_direction)*materials[i][j+1].activated
						+ materials[i][j].emf*Math.sin(materials[i][j].emf_direction)*materials[i][j].activated);
			}
		}
	}

	public void computeFreeEnergy(int i, int j, double ni, double W, double B)
	{
		double K = ni*ni;

		E0_n[i][j] = -W + 0.5*B;
		E0_p[i][j] = W + 0.5*B;

		double TS_n = (Math.log(K)/beta + E0_n[i][j] + E0_p[i][j])/2.0;
		double TS_p = TS_n;

		F0_n[i][j] = E0_n[i][j] - TS_n;
		F0_p[i][j] = E0_p[i][j] - TS_p;
	}

	public double calcEquilibriumElectronCharge(double rho_back, double K) {
		double B = rho_back/e_charge;
		return -e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}

	public double calcEquilibriumHoleCharge(double rho_back, double K) {
		double B = -rho_back/e_charge;
		return e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}

	// Mark all points on grid that are connected to given point and within a specified distance using Dijkstra's algorithm
	public void markNeighborhood(int i1, int j1, int max_distance) {

		for (int i = i1 - max_distance; i <= i1+max_distance; i++) {
			for (int j = j1 - max_distance; j <= j1+max_distance; j++) {
				if (i >= 0 && j >= 0 && i < nx && j < ny) {
					distance[i][j] = Integer.MAX_VALUE;
					visited[i][j] = false;
				}
			}
		}
		distance[i1][j1] = 0;

		while (true) {
			int smallest_length = Integer.MAX_VALUE;
			int smallest_i = 0;
			int smallest_j = 0;
			for (int i = i1 - max_distance; i <= i1+max_distance; i++) {
				for (int j = j1 - max_distance; j <= j1+max_distance; j++) {
					if (i >= 0 && j >= 0 && i < nx && j < ny) {
						if (distance[i][j] < smallest_length && !visited[i][j] && conducting[i][j] == 1) {
							smallest_length = distance[i][j];
							smallest_i = i;
							smallest_j = j;
						}
					}
				}
			}

			if (smallest_length == Integer.MAX_VALUE)
				break;

			for (int di = -1; di <= 1; di ++) {
				for (int dj = -1; dj <= 1; dj ++) {
					if (Math.abs(di)+Math.abs(dj) > 0 /* == 1 for vN neighborhood */
							&& smallest_i+di >= 0 && smallest_j+dj >= 0 && smallest_i+di < nx && smallest_j+dj < ny
							&& distance[smallest_i][smallest_j] < max_distance
							&& conducting[smallest_i+di][smallest_j+dj] == 1) {
						distance[smallest_i+di][smallest_j+dj] = Math.min(distance[smallest_i+di][smallest_j+dj], distance[smallest_i][smallest_j]+1);
					}
				}
			}

			visited[smallest_i][smallest_j] = true;
		}
	}

	public double[][] copyArray(double[][] array) {
		double [][] newarray = new double[array.length][];
		for(int i = 0; i < array.length; i++)
			newarray[i] = array[i].clone();
		return newarray;
	}
	
	public void eraseMaterial(int i, int j) {
		materials[i][j].erase();
		rho_p[i][j] = 0;
		rho_n[i][j] = 0;
	}

	public void initializeAllMaterials() {
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				initializeMaterial(i, j);
			}
		}
	}

	public void initializeMaterial(int i, int j) {
		initializeMaterial(i, j, materials[i][j].type);
	}

	public void initializeMaterial(int i, int j, MaterialType material) {
		if (i < 0 || j < 0 || i >= nx || j >= ny || materials[i][j].modified)
			return;

		materials[i][j].type = material;

		if (material == MaterialType.DIELECTRIC) materials[i][j].eps_r = dielectric_eps_r;
		else if (material == MaterialType.FERROMAGNET) materials[i][j].mu_r = ferromagnet_mu_r;
		else if (material == MaterialType.POS_CHARGE) materials[i][j].rho_back = staticcharge_density;
		else if (material == MaterialType.NEG_CHARGE) materials[i][j].rho_back = -staticcharge_density;

		if (MaterialType.isConducting(material))
		{
			materials[i][j].conducting = 1;
			materials[i][j].ni = ni_metal;
			materials[i][j].W = W_metal_default;
			materials[i][j].Eb = E_b_metal;
			materials[i][j].Ea = E_a_metal;

			if (material == MaterialType.METAL_HIGH_W) materials[i][j].W = W_metal_high;
			else if (material == MaterialType.METAL_LOW_W) materials[i][j].W = W_metal_low;
			else if (material == MaterialType.METAL_HIGH_C) materials[i][j].ni = ni_metal_high;
			else if (material == MaterialType.METAL_LOW_C) materials[i][j].ni = ni_metal_low;
		}

		if (MaterialType.isSemiconducting(material))
		{
			materials[i][j].conducting = 1;
			materials[i][j].semiconducting = 1;
			materials[i][j].ni = ni_semi;
			materials[i][j].W = W_semi;
			materials[i][j].Eb = E_b_semi;
			materials[i][j].Ea = E_a_semi;

			if (material == MaterialType.SEMI_P_TYPE) materials[i][j].rho_back = -p_default_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_N_TYPE) materials[i][j].rho_back = n_default_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_HEAVY_P_TYPE) materials[i][j].rho_back = -p_heavy_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_HEAVY_N_TYPE) materials[i][j].rho_back = n_heavy_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_LIGHT_P_TYPE) materials[i][j].rho_back = -p_light_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_LIGHT_N_TYPE) materials[i][j].rho_back = n_light_doping_concentration*e_charge;
		}

		materials[i][j].auto_placed = false;
	}

	public void checkCFL() {
		System.out.println("Wave equation CFL ratio = " + (Math.sqrt(2)*c)/(ds/dt_maximum));
		System.out.println("Electron diffusion CFL ratio = " + dt_maximum/(ds*ds/(4*D_electron)));
		System.out.println("Hole diffusion CFL ratio = " + dt_maximum/(ds*ds/(4*D_hole)));
	}


	public double calcCharge(ChargeProbe p) {
		double Q = 0;

		int n_min = 0;
		int n_max = 0;
		int m_min = 0;
		int m_max = 0;

		n_min = Math.min(p.x1, p.x2);
		n_max = Math.max(p.x1, p.x2);
		m_min = Math.min(p.y1, p.y2);
		m_max = Math.max(p.y1, p.y2);

		for (int n = n_min; n <= n_max; n++) {
			for (int m = m_min; m <= m_max; m++) {
				Q += rho_free[n][m]*(ds*ds);
			}
		}

		return Q;
	}

	public double calcCurrent(int x0, int y0, int x1, int y1) {
		int dy = y1 - y0;
		int dx = x1 - x0;
		float t = (float) 0.5;
		float J = 0;

		if (Math.abs(dx) > Math.abs(dy)) {
			float m = (float) dy / (float) dx;
			t += y0;
			dx = (dx < 0) ? -1 : 1;
			m *= dx;
			while (x0 != x1) {
				int x0_prev = x0;
				float t_prev = t;

				x0 += dx;
				t += m;

				J += accumCurrent(x0_prev, (int)t_prev, x0, (int)t, Jx_free, Jy_free);

			}
		} else {
			float m = (float) dx / (float) dy;
			t += x0;
			dy = (dy < 0) ? -1 : 1;
			m *= dy;
			while (y0 != y1) {
				int y0_prev = y0;
				float t_prev = t;

				y0 += dy;
				t += m;

				J += accumCurrent((int)t_prev, y0_prev, (int)t, y0, Jx_free, Jy_free);
			}
		}

		return J;
	}

	public double accumCurrent(int x0, int y0, int x1, int y1, double[][] Jx, double[][] Jy) {
		int dx = x1-x0;
		int dy = y1-y0;
		if (dx == 1 && dy == 0) {
			return -Jy[x0+1][y0]*ds;
		} else if (dx == -1 && dy == 0) {
			return Jy[x0][y0]*ds;
		} else if (dx == 0 && dy == 1) {
			return Jx[x0][y0+1]*ds;
		} else if (dx == 0 && dy == -1) {
			return -Jx[x0][y0]*ds;
		} else if (dx == 1 && dy == 1) {
			return (Jx[x0][y0+1] - Jy[x0+1][y0+1])*ds;
		} else if (dx == 1 && dy == -1) {
			return (-Jx[x0][y0] - Jy[x0+1][y0-1])*ds;
		} else if (dx == -1 && dy == 1) {
			return (+Jx[x0][y0+1] + Jy[x0][y0+1])*ds;
		} else if (dx == -1 && dy == -1) {
			return (-Jx[x0][y0] + Jy[x0][y0-1])*ds;
		}
		return 0;
	}

	public void prescaleDielectric() {
		downscale_x_vector(epsx, MG_epsx, log2_resolution);
		downscale_y_vector(epsy, MG_epsy, log2_resolution);
	}

	public void multigridSolve(boolean correctEfield, boolean computePhi) {
		assert(!(correctEfield && computePhi));

        poissonLock.lock();
        try {
			if (correctEfield)
				t4.start();
			if (computePhi)
				t7.start();

			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					MG_phi1[i][j] = 0;
					MG_phi2[i][j] = 0;
					for (int k = 0; k <= log2_resolution; k++) {
						MG_rho[k][i][j] = 0;
					}
				}
			}

			if (correctEfield)
			{
				for (int i = 1; i < nx-1; i++) {
					for (int j = 1; j < ny-1; j++) {
						MG_rho0[i][j] = ((Ex[i][j]*epsx[i][j]-Ex[i-1][j]*epsx[i-1][j] + Ey[i][j]*epsy[i][j]-Ey[i][j-1]*epsy[i][j-1])/ds) - rho_free[i][j];
						debug[i][j] = MG_rho0[i][j];
					}
				}

				double num = 0;
				double denom = 0;
				for (int i = 1; i < nx-1; i++) {
					for (int j = 1; j < ny-1; j++) {
						num += MG_rho0[i][j]*MG_rho0[i][j];
						denom += rho_free[i][j]*rho_free[i][j];
					}
				}
				System.out.println("Starting poisson residual: " + Math.sqrt(num/denom));
			}
			if (computePhi) {
				for (int i = 1; i < nx-1; i++) {
					for (int j = 1; j < ny-1; j++) {
						MG_rho0[i][j] = (Ex[i][j]-Ex[i-1][j] + Ey[i][j]-Ey[i][j-1])/(ds) + (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]-4*phi[i][j])/(ds*ds);
					}
				}
			}

			//int[] stepsarray = {0, 0, 200, 200, 200, 200, 200, 50, 20};
			int[] stepsarray = {0, 0, 100, 100, 100, 100, 50, 25, 20};

			if (log2_resolution <= 1)
				throw new RuntimeException("Grid size too small!");

			downscale(MG_rho0, MG_rho, log2_resolution);

			for (int fineness = 2; fineness <= log2_resolution; fineness++) {
				int nx_tmp = (1 << fineness);
				int ny_tmp = (1 << fineness);
				double gridsize = width/(1 << fineness);

				int poissonsteps = stepsarray[(fineness > 8)? 8 : fineness];
				double alpha = (gridsize*gridsize);

				JacobiIteration(poissonsteps, nx_tmp, ny_tmp, alpha, fineness, computePhi);

				if (fineness == log2_resolution)
					break;

				for (int i = 0; i < nx_tmp; i++) {
					for (int j = 0; j < ny_tmp; j++) {
						MG_phi2[i][j] = MG_phi1[i][j];
					}
				}

				for (int i = 0; i < nx_tmp; i++)
				{
					for (int j = 0; j < nx_tmp; j++) {
						MG_phi1[2*i][2*j] = MG_phi2[i][j];
						MG_phi1[2*i+1][2*j] = MG_phi2[i][j];
						MG_phi1[2*i][2*j+1] = MG_phi2[i][j];
						MG_phi1[2*i+1][2*j+1] = MG_phi2[i][j];
					}
				}

				for (int i = 0; i < 2*nx_tmp; i++) {
					for (int j = 0; j < 2*ny_tmp; j++) {
						MG_phi2[i][j] = MG_phi1[i][j];
					}
				}
			}

			if (correctEfield) {
				for (int i = 0; i < nx-1; i++)
				{
					for (int j = 0; j < ny-1; j++)
					{
						Ex[i][j] = Ex[i][j] + (MG_phi1[i+1][j]-MG_phi1[i][j])/ds;
						Ey[i][j] = Ey[i][j] + (MG_phi1[i][j+1]-MG_phi1[i][j])/ds;
					}
				}
			}

			if (computePhi)
			{
				for (int i = 0; i < nx; i++)
				{
					for (int j = 0; j < ny; j++)
					{
						phi[i][j] = phi[i][j] + MG_phi1[i][j];
					}
				}
			}

			if (correctEfield) {

				double num = 0;
				double denom = 0;
				for (int i = 1; i < nx-1; i++) {
					for (int j = 1; j < ny-1; j++) {
						double drho = ((Ex[i][j]*epsx[i][j]-Ex[i-1][j]*epsx[i-1][j] + Ey[i][j]*epsy[i][j]-Ey[i][j-1]*epsy[i][j-1])/ds) - rho_free[i][j];
						num += drho*drho;
						denom += rho_free[i][j]*rho_free[i][j];
					}
				}

				System.out.println("Poisson residual: " + Math.sqrt(num/denom));
			}

			if (correctEfield)
				t4.stop();

			if (computePhi)
				t7.stop();
		} finally {
            poissonLock.unlock();
        }
	}

	public void JacobiIteration(int steps, int xmax, int ymax, double alpha, int fineness, boolean calcPhi) {

		if (calcPhi) {
			for (int poissonit = 0; poissonit < steps; poissonit++) {
				for (int i = 1; i < xmax-1; i++) {
					for (int j = 1; j < ymax-1; j++) {
						MG_phi2[i][j] = (0.1*MG_phi1[i][j] + ((MG_phi1[i-1][j] + MG_phi1[i+1][j] + MG_phi1[i][j-1] + MG_phi1[i][j+1]) + MG_rho[fineness][i][j]*alpha)/4.0)/1.1;
					}
				}
				for (int i = 1; i < xmax-1; i++) {
					for (int j = 1; j < ymax-1; j++) {
						MG_phi1[i][j] = (0.1*MG_phi2[i][j] + ((MG_phi2[i-1][j] + MG_phi2[i+1][j] + MG_phi2[i][j-1] + MG_phi2[i][j+1]) + MG_rho[fineness][i][j]*alpha)/4.0)/1.1;
					}
				}
			}
		} else {
			for (int i = 1; i < xmax-1; i++) {
				for (int j = 1; j < ymax-1; j++) {
					MG_eps_avg[i][j] = (MG_epsx[fineness][i-1][j]+MG_epsx[fineness][i][j]+MG_epsy[fineness][i][j-1]+MG_epsy[fineness][i][j]);
				}
			}

			for (int poissonit = 0; poissonit < steps; poissonit++) {
				for (int i = 1; i < xmax-1; i++) {
					for (int j = 1; j < ymax-1; j++) {
						MG_phi2[i][j] = (0.1*MG_phi1[i][j] + ((MG_phi1[i-1][j]*MG_epsx[fineness][i-1][j]
								+ MG_phi1[i+1][j]*MG_epsx[fineness][i][j]
								+ MG_phi1[i][j-1]*MG_epsy[fineness][i][j-1]
								+ MG_phi1[i][j+1]*MG_epsy[fineness][i][j])
								+ MG_rho[fineness][i][j]*alpha)/MG_eps_avg[i][j])/1.1;
					}
				}
				for (int i = 1; i < xmax-1; i++) {
					for (int j = 1; j < ymax-1; j++) {
						MG_phi1[i][j] = (0.1*MG_phi2[i][j] + ((MG_phi2[i-1][j]*MG_epsx[fineness][i-1][j]
								+ MG_phi2[i+1][j]*MG_epsx[fineness][i][j]
								+ MG_phi2[i][j-1]*MG_epsy[fineness][i][j-1]
								+ MG_phi2[i][j+1]*MG_epsy[fineness][i][j])
								+ MG_rho[fineness][i][j]*alpha)/MG_eps_avg[i][j])/1.1;
					}
				}
			}
		}
	}

	public void downscale(double[][] source, double[][][] dest, int steps) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				dest[steps][i][j] = source[i][j];
			}
		}
		int nx_d = nx;
		int ny_d = ny;
		for (int k = steps-1; k >= 0; k--) {
			nx_d = nx_d/2;
			ny_d = ny_d/2;

			for (int i = 0; i < nx_d; i++) {
				for (int j = 0; j < ny_d; j++) {
					dest[k][i][j] = 0.25*(dest[k+1][2*i][2*j]+dest[k+1][2*i+1][2*j]+dest[k+1][2*i][2*j+1]+dest[k+1][2*i+1][2*j+1]);
				}
			}
		}
	}

	public void downscale_x_vector(double[][] source, double[][][] dest, int steps) {
		for (int i = 0; i < nx-1; i++) {
			for (int j = 0; j < ny; j++) {
				dest[steps][i][j] = source[i][j];
			}
		}
		int nx_d = nx;
		int ny_d = ny;
		for (int k = steps-1; k >= 0; k--) {
			nx_d = nx_d/2;
			ny_d = ny_d/2;

			for (int i = 0; i < nx_d-1; i++) {
				for (int j = 0; j < ny_d; j++) {
					dest[k][i][j] = 0.125*(dest[k+1][2*i][2*j]+dest[k+1][2*i][2*j+1]
							+2*dest[k+1][2*i+1][2*j]+2*dest[k+1][2*i+1][2*j+1]
							+dest[k+1][2*i+2][2*j]+dest[k+1][2*i+2][2*j+1]);
				}
			}
		}
	}

	public void downscale_y_vector(double[][] source, double[][][] dest, int steps) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny-1; j++) {
				dest[steps][i][j] = source[i][j];
			}
		}
		int nx_d = nx;
		int ny_d = ny;
		for (int k = steps-1; k >= 0; k--) {
			nx_d = nx_d/2;
			ny_d = ny_d/2;

			for (int i = 0; i < nx_d; i++) {
				for (int j = 0; j < ny_d-1; j++) {
					dest[k][i][j] = 0.125*(dest[k+1][2*i][2*j]+dest[k+1][2*i+1][2*j]
							+2*dest[k+1][2*i][2*j+1]+2*dest[k+1][2*i+1][2*j+1]
							+dest[k+1][2*i][2*j+2]+dest[k+1][2*i+1][2*j+2]);
				}
			}
		}
	}
}

enum BoundaryCondition {
	DISSIPATIVE("Absorbing boundary"),
	CONDUCTING("Conducting boundary");

	String name;
	BoundaryCondition(String name)
	{
		this.name = name;
	}

	@Override
	public String toString() {
		return name;
	}
}

class CurrentProbe {
	int x1 = 0;
	int y1 = 0;
	int x2 = 0;
	int y2 = 0;

	double current = 0;
}

class VoltageProbe {
	int x = 0;
	int y = 0;

	double potential = 0;
}

class ChargeProbe {
	int x1;
	int y1;

	int x2;
	int y2;

	double charge = 0;
}

enum Species {
	ELECTRON, HOLE;
}

enum MaterialType
{

	EMF					("Voltage source (Adjustable)",			230, 216, 46, 230),
	AC_EMF				("AC voltage source (Adjustable)",		230, 150, 216, 230),
	SWITCH				("Switch",								194, 194, 194, 120),
	METAL				("Metal",								153, 153, 153, 130),
	METAL_HIGH_C		("Conductive metal",					191, 191, 191, 130),
	METAL_LOW_C			("Resistive metal",						94, 94, 94, 130),
	METAL_HIGH_W		("High workfunction metal",				163, 116, 116, 130),
	METAL_LOW_W			("Low workfunction metal",				116, 121, 163, 130),
	SEMI				("Intrinsic semiconductor",				207,  161, 212, 110),
	SEMI_P_TYPE			("P-type semiconductor",				191,  74,  34, 110),
	SEMI_N_TYPE			("N-type semiconductor",				 84, 123, 191, 110),
	SEMI_HEAVY_P_TYPE	("Heavily doped P-type semiconductor",	204,  41,  41, 110),
	SEMI_HEAVY_N_TYPE	("Heavily doped N-type semiconductor",	 39,  52, 194, 110),
	SEMI_LIGHT_P_TYPE	("Lightly doped P-type semiconductor",	201, 131,  73, 110),
	SEMI_LIGHT_N_TYPE	("Lightly doped N-type semiconductor",	137, 188, 204, 110),
	DIELECTRIC			("Dielectric",							 81, 171,  51, 80),
	FERROMAGNET			("Ferromagnet",							116, 50, 117, 80),
	POS_CHARGE			("Positive static charge",				116, 50, 50, 80),
	NEG_CHARGE			("Negative static charge",				50, 50, 117, 80),
	DECO				("Decoration",							255, 255, 255, 255),
	ABSORBER			("Absorber",							 50,  50,  50),
	VACUUM				("Vacuum",								 20,  20,  20);

	String name;
	int color_r;
	int color_g;
	int color_b;
	int color_grayscale;

	MaterialType(String name, int r, int g, int b) {
		this.name = name;
		color_r = r;
		color_g = g;
		color_b = b;
		color_grayscale = (int)(0.7*Math.max(Math.max(r, g), b));
	}

	MaterialType(String name, int r, int g, int b, int grayscale_brightness) {
		this.name = name;
		color_r = r;
		color_g = g;
		color_b = b;
		color_grayscale = grayscale_brightness;
	}

	public static boolean isConducting(MaterialType material) {
		return (material == MaterialType.EMF
				|| material == MaterialType.AC_EMF
				|| material == MaterialType.SWITCH
				|| material == MaterialType.METAL
				|| material == MaterialType.METAL_HIGH_W
				|| material == MaterialType.METAL_LOW_W
				|| material == MaterialType.METAL_HIGH_C
				|| material == MaterialType.METAL_LOW_C);
	}

	public static boolean isSemiconducting(MaterialType material) {
		return (material == MaterialType.SEMI_P_TYPE
				|| material == MaterialType.SEMI_N_TYPE
				|| material == MaterialType.SEMI
				|| material == MaterialType.SEMI_HEAVY_P_TYPE
				|| material == MaterialType.SEMI_HEAVY_N_TYPE
				|| material == MaterialType.SEMI_LIGHT_P_TYPE
				|| material == MaterialType.SEMI_LIGHT_N_TYPE);
	}

	@Override
	public String toString() {
		return "Material: " + name;
	}
}

class Material implements Cloneable {
	MaterialType type = MaterialType.VACUUM;

	boolean modified = false;
	boolean auto_placed = true;
	int activated = 1;
	int conducting = 0;
	int semiconducting = 0;
	double emf = 0.0;			// EMF strength
	double emf_direction = 0.0;	// EMF direction
	double eps_r = 1.0;			// Permittivity
	double mu_r = 1.0;			// Permeability
	double rho_back = 0.0;		// Background charge density
	double ni = 0;				// Equilibrium carrier density
	double W = 0;				// Work function
	double Eb = 0;				// Bandgap
	double Ea = 0;				// Recombination activation energy
	double absorptivity = 0.0;

	public void erase() {
		type = MaterialType.VACUUM;
		modified = false;
		auto_placed = true;
		activated = 1;
		conducting = 0;
		semiconducting = 0;
		emf = 0.0;
		emf_direction = 0.0;
		eps_r = 1.0;
		mu_r = 1.0;
		rho_back = 0.0;
		ni = 0;
		W = 0;
		Eb = 0;
		Ea = 0;
		absorptivity = 0;
	}

    @Override
    public Material clone() {
        try {
			return (Material) super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
    }
}

class ClipboardMaterial implements Cloneable {
	double rho_n = 0;
	double rho_p = 0;
	Material m = new Material();
	
	public ClipboardMaterial() {};

	public ClipboardMaterial(Material m) {
		this.m = m.clone();
	}
	
    public ClipboardMaterial(Simulation e, int i, int j) {
    	this(e.materials[i][j]);
    	rho_n = e.rho_n[i][j];
    	rho_p = e.rho_p[i][j];
    }
    
    public void paste(Simulation e, int i, int j) {
    	e.materials[i][j] = m.clone();
    	e.rho_n[i][j] = rho_n;
    	e.rho_p[i][j] = rho_p;
    }
	
	public void erase() {
		m.erase();
		rho_n = 0;
		rho_p = 0;
	}
	
    @Override
    public ClipboardMaterial clone() {
        try {
			return (ClipboardMaterial) super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
    }
}