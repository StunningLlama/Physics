// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.TimerTask;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.plot.BandPlot;
import electrodynamics.plot.CarrierPlot;
import electrodynamics.plot.Plot;
import electrodynamics.plot.ProbePlot;
import electrodynamics.plot.ScalarPlot;
import electrodynamics.probe.ChargeProbe;
import electrodynamics.probe.CurrentProbe;
import electrodynamics.probe.VoltageProbe;
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
	
	public Renderer.RenderCanvas canvas;
	public MainWindow opts;
	public AdvancedOptions adv_opts;
	public Renderer renderer;
	public Controls controls;
	public SaveManager savemanager;
	
	public ArrayList<Plot> plots = new ArrayList<>();
	public BandPlot bandplot;
	public ScalarPlot scalarplot;
	public CarrierPlot carrierplot;
	public ProbePlot voltageprobeplot;
	public ProbePlot currentprobeplot;
	public ProbePlot chargeprobeplot;
	
	
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

	public long stepnumber = 0;
	public long frame = 0;
	public int lastsimspeed = 0;
	public int iteration_multiplier = 0;

	public double error_detection_threshold = 1e-5;
	public boolean sign_violation = false;
	public boolean numerical_overflow = false;
	
	public boolean advsettings_tweaked = false;
	
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

	public double ni_currentsource = ni_metal;
	public double currentsource_mobility = 0.0002;
	public double currentsource_sigma = ni_currentsource*e_charge*(mu_electron + mu_hole)*currentsource_mobility;

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
	
	public double time = 0.0;
	public double AC_phase;
	
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

	public double[][] emfx;					// External electromotive force
	public double[][] emfy;
	public double[][] epsx;					// Dielectric constant
	public double[][] epsy;
	public double[][] mu_z;					// Relative permeability

	public double[][] relative_mobility;	// Mobility relative to metal
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
		canvas = renderer.new RenderCanvas(this);
		canvas.setFocusable(true);
		adv_opts = new AdvancedOptions();
		datafile = new File(datafilename);

		bandplot = new BandPlot(); plots.add(bandplot);
		scalarplot = new ScalarPlot(); plots.add(scalarplot);
		carrierplot = new CarrierPlot(); plots.add(carrierplot);
		voltageprobeplot = new ProbePlot("Voltage probe plot", "Voltage [mV]", "V", 1e3); plots.add(voltageprobeplot);
		currentprobeplot = new ProbePlot("Current probe plot", "Current [mA]", "I", 1e3); plots.add(currentprobeplot);
		chargeprobeplot = new ProbePlot("Charge probe plot", "Charge [fC]", "Q", 1e15); plots.add(chargeprobeplot);
		
		for (Plot p : plots)
			p.initalize();

		SemiSim.detect64Bit();

		setSize(default_resolution, default_width);
		resetFields(true);
		
		opts.initialize(this);
		adv_opts.initialize(this);

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
    				opts.setTitle(SemiSim.name);
    				opts.textPane.setText("Description of simulation");
    				SaveManager.currentfile = null;
    				controls.reset = false;
    			}

    			lastsimspeed = opts.gui_simspeed.getValue();
    			dt = dt_maximum*(lastsimspeed/20.0);

    			if (controls.undo) {
    				controls.undoredo.undo(this);
    				controls.undo = false;
    			}
    			
    			if (controls.redo) {
    				controls.undoredo.redo(this);
    				controls.redo = false;
    			}
    			
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

    			setDebugInfo();
    			
    			if (numerical_overflow)
    				opts.gui_paused.setSelected(true);


    			if (controls.save) {
    				savemanager.writeFile(false);
    				controls.save = false;
    			}
    			
    			if (controls.saveas) {
    				savemanager.writeFile(true);
    				controls.saveas = false;
    			}

    			if (controls.load) {
    				savemanager.readFile();
    				controls.load = false;
    			}

    			if (controls.exit) {
    				SwingUtilities.invokeLater(() -> {
    					controls.windowClosing(null);
    				});
    				controls.exit = false;
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
			relative_mobility = new double[nx][ny];
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

						relative_mobility[i][j] = 1.0;

						conducting[i][j] = 0;
						conducting_x[i][j] = 0;
						conducting_y[i][j] = 0;

						ac_x[i][j] = 0;
						ac_y[i][j] = 0;

						absorptivity[i][j] = 0;
						absorptivity_x[i][j] = 0;
						absorptivity_y[i][j] = 0;

						controls.selected[i][j] = false;
						controls.selected_EMF[i][j] = false;
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
				}
			}
			
			if (resetall) {
				controls.EMF_selected = false;
				controls.changesmade = false;
				opts.gui_bc.setSelectedItem(BoundaryCondition.DISSIPATIVE);
				controls.prev_boundary = BoundaryCondition.DISSIPATIVE;
			}

			numerical_overflow = false;
			sign_violation = false;
			
			if (resetall) {
				voltageprobes.clear();
				currentprobes.clear();
				chargeprobes.clear();
				ground = null;

				for (Plot p: plots) {
					p.frame.setVisible(false);
				}

				controls.undoredo.resetUndoHistory(this);
				controls.resetZoom();
			}
			
			for (VoltageProbe p : voltageprobes) p.data.resetData();
			for (CurrentProbe p : currentprobes) p.data.resetData();
			for (ChargeProbe p : chargeprobes) p.data.resetData();
			if (ground != null) ground.data.resetData();

			renderer.resetChargeDots();

			initializeAllMaterials();
			updateAllMaterials(true);
			controls.undoredo.captureState(this);
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
								mobility_factor[i][j] = Math.min(relative_mobility[i][j], E_sat/E_avg);
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
		ScalarView view_scalar = controls.scalarview.getOption();

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
					grad_E0x_n[i][j] = conducting_x[i][j] == 1? -conducting_x[i][j]*(E0_n[i+1][j] - E0_n[i][j])/ds : 0;
					grad_E0x_p[i][j] = conducting_x[i][j] == 1? -conducting_x[i][j]*(E0_p[i+1][j] - E0_p[i][j])/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					grad_E0y_n[i][j] = conducting_x[i][j] == 1? -conducting_y[i][j]*(E0_n[i][j+1] - E0_n[i][j])/ds : 0;
					grad_E0y_p[i][j] = conducting_x[i][j] == 1? -conducting_y[i][j]*(E0_p[i][j+1] - E0_p[i][j])/ds : 0;
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
			ground.measure(this, false);

		for (VoltageProbe p: voltageprobes) {
			p.measure(this, frame%controls.plotinterval == 0);
		}

		for (CurrentProbe p: currentprobes) {
			p.measure(this, frame%controls.plotinterval == 0);
		}
		
		for (ChargeProbe p: chargeprobes) {
			p.measure(this, frame%controls.plotinterval == 0);
		}

		if (controls.logdata) {
			logProbeData();
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
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 0) {
					F0_n[i][j] = Double.NaN;
					F0_p[i][j] = Double.NaN;
					E0_n[i][j] = Double.NaN;
					E0_p[i][j] = Double.NaN;
				}
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
				relative_mobility[i][j] = (materials[i][j].type == MaterialType.CURRENT)? currentsource_mobility : 1;
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

		if (material.isConducting())
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
			else if (material == MaterialType.CURRENT) materials[i][j].ni = ni_currentsource;
		}

		if (material.isSemiconducting())
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
	
	public void logProbeData() {
		String data = "";
		data += ("t = " + Utils.getSI(time, "s"));

		int index = 0;
		for (VoltageProbe p: voltageprobes) {
			data += (", V" + getProbeName(index) + " = " + Utils.getSI(p.potential, "V"));
			index++;
		}
		
		index = 0;
		for (CurrentProbe p: currentprobes) {
			data += (", I" + getProbeName(index) + " = " + Utils.getSI(p.current, "A"));
			index++;
		}
		
		index = 0;
		for (ChargeProbe p: chargeprobes) {
			data += (", Q" + getProbeName(index) + " = " + Utils.getSI(p.charge, "C"));
			index++;
		}

		data += "\n";

		datastream.print(data);
		datastream.flush();

		renderer.probetexttimer = 30;
	}
	
	public String prevDesc = "";
	
	public void setDebugInfo() {
		if (!controls.debugging) {
			if (prevDesc != "") {
				opts.textPane.setText(prevDesc);
				prevDesc = "";
			}

			return;
		}
		
		if (prevDesc == "")
			prevDesc = opts.textPane.getText();
		
		int mx = controls.mx;
		int my = controls.my;
		Material mat = materials[mx][my];
		
		String str = "";
		str += "Mouse\n";
		str += ("x\t"  						+	Utils.getSI(mx*ds, "m") + "\n");
		str += ("y\t"  						+	Utils.getSI(ds*ny-(my+1)*ds, "m") + "\n");
		str += "\nFields\n";
		str += ("E\t"  						+	Utils.getSI(Utils.bilinearinterp_length(Ex, Ey, mx, my, nx, ny), "V/m") + "\n");
		str += ("D\t"  						+	Utils.getSI(Utils.bilinearinterp_length(Dx, Dy, mx, my, nx, ny), "C/m^2") + "\n");
		str += ("B\t"  						+	Utils.getSI(parity*Utils.bilinearinterp(Bz, mx-0.5, my-0.5, nx, ny), "T") + "\n");
		str += ("H\t"  						+	Utils.getSI(parity*Utils.bilinearinterp(Hz, mx-0.5, my-0.5, nx, ny), "A/m") + "\n");
		str += "\nCharge/current\n";
		str += ("\u03c1\u2099\t"  			+	Utils.getSI(rho_n[mx][my], "C/m^3") + "\n");
		str += ("\u03c1\u209A\t"  			+	Utils.getSI(rho_p[mx][my], "C/m^3") + "\n");
		str += ("\u03c1\u2080\t"			+	Utils.getSI(rho_back[mx][my], "C/m^3") + "\n");
		str += ("\u03c1 abs\t"				+	Utils.getSI(rho_abs[mx][my], "C/m^3") + "\n");
		str += ("\u03c1\t"  				+	Utils.getSI(rho_free[mx][my], "C/m^3") + "\n");
		str += ("J\u2099\t" 				+	Utils.getSI(Utils.bilinearinterp_length(Jx_n, Jy_n, mx, my, nx, ny), "A/m^2") + "\n");
		str += ("J\u209A\t"  				+	Utils.getSI(Utils.bilinearinterp_length(Jx_p, Jy_p, mx, my, nx, ny), "A/m^2") + "\n");
		str += ("J abs\t"  					+	Utils.getSI(Utils.bilinearinterp_length(Jx_abs, Jy_abs, mx, my, nx, ny), "A/m^2") + "\n");
		str += ("J\t"  						+	Utils.getSI(Utils.bilinearinterp_length(Jx_free, Jy_free, mx, my, nx, ny), "A/m^2") + "\n");
		str += "\nThermodynamic quantities\n";
		str += ("S\t"  						+	Utils.getSI(Utils.bilinearinterp_length(Sx, Sy, mx, my, nx, ny), "W/m^2") + "\n");
		str += ("u\t"  						+	Utils.getSI(u[mx][my], "J/m^3") + "\n");
		str += ("\u03d5\t"  				+	Utils.getSI(phi[mx][my], "V") + "\n");
		str += ("F\u2099\t"  				+	Utils.getSI(F_n[mx][my]/q_n + phi[mx][my] - W_semi/eVtoJ, "V") + "\n");
		str += ("F\u209a\t"  				+	Utils.getSI(F_p[mx][my]/q_p + phi[mx][my] - W_semi/eVtoJ, "V") + "\n");
		str += ("F\t"  						+	Utils.getSI(F[mx][my], "V") + "\n");
		str += ("G\t"  						+	Utils.getSI(G[mx][my], "/m^3 s") + "\n");
		str += ("R\t"  						+	Utils.getSI(R[mx][my], "/m^3 s") + "\n");
		str += ("Electron CMF\t"  			+	Utils.getSI(Utils.bilinearinterp_length(cmfy_n, cmfy_n, mx, my, nx, ny), "J/m") + "\n");
		str += ("Hole CMF\t"  				+	Utils.getSI(Utils.bilinearinterp_length(cmfy_p, cmfy_p, mx, my, nx, ny), "J/m") + "\n");
		str += ("EMF\t"  					+	Utils.getSI(Utils.bilinearinterp_length(emfx, emfy, mx, my, nx, ny), "V/m") + "\n");
		str += "\nMaterial properties\n";
		str += ("Type\t"  					+	mat.type.name + "\n");
		str += ("\u2130\t"  				+	Utils.getSI(mat.emf, "V/m") + "\n");
		str += ("\u03b5/\u03b5\u2080\t" 	+	Utils.getSI(mat.eps_r, "") + "\n");
		str += ("\u03bc/\u03bc\u2080\t"  	+	Utils.getSI(mat.mu_r, "") + "\n");
		str += ("ni\t"  					+	Utils.getSI(mat.ni, "/m^3") + "\n");
		str += ("W\t"  						+	Utils.getSI(mat.W, "eV") + "\n");
		str += ("Eb\t"  					+	Utils.getSI(mat.Eb, "eV") + "\n");
		str += ("Ea\t"  					+	Utils.getSI(mat.Ea, "J") + "\n");
		str += ("Absorptivity\t"  			+	Utils.getSI(mat.absorptivity, "") + "\n");
		str += ("Conducting\t"  			+	mat.conducting + "\n");
		str += ("Semicond.\t"  				+	mat.semiconducting + "\n");
		
		opts.textPane.setText(str);
	}
	
	public String getProbeName(int index) {
		if (index < 26) return String.valueOf((char)('a'+index));
		return String.valueOf(index - 26);
	}
	
	public enum BoundaryCondition {
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

	public enum Species {
		ELECTRON, HOLE;
	}
}