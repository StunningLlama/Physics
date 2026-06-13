// Copyright (c) Brandon Li 2025-2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
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

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import electrodynamics.Renderer.ScalarView;
import electrodynamics.Renderer.VectorView;
import electrodynamics.gui.AdvancedOptions;
import electrodynamics.gui.MainWindow;
import electrodynamics.gui.MaterialManager;
import electrodynamics.gui.MaterialViewer;
import electrodynamics.gui.Preferences;
import electrodynamics.plot.BandPlot;
import electrodynamics.plot.CarrierPlot;
import electrodynamics.plot.Plot;
import electrodynamics.plot.ProbePlot;
import electrodynamics.plot.ScalarPlot;
import electrodynamics.plot.XYPlot;
import electrodynamics.probe.ChargeProbe;
import electrodynamics.probe.CurrentProbe;
import electrodynamics.probe.FluxProbe;
import electrodynamics.probe.Ground;
import electrodynamics.probe.Probe;
import electrodynamics.probe.Ruler;
import electrodynamics.probe.VoltageProbe;
import electrodynamics.units.Quantity;
import electrodynamics.units.Units;
import electrodynamics.util.FastExp;
import electrodynamics.util.PeriodicTask;
import electrodynamics.util.Timer;
import electrodynamics.util.Utils;

public class Simulation extends PeriodicTask {
	/*
	 * Demos:
	 * 	- PN diode (LED) [done]
	 * 	- Schottky diode [done]
	 * 	- MOSFET (p-channel, n-channel), depletion & enhancement [done]
	 * 	- BJT (NPN, PNP) [done]
	 * 	- JFET (p-channel, n-channel) [done]
	 *  - IGBT - [done, okay]
	 *  - MESFET [done]
	 *  - SCR  [done]
	 *  - Darlington pair [done]
	 *  - LEDs
	 *  - Saturation velocity
	 *  - Recombination
	 */
	
	/* Parts */
	
	public Renderer.RenderCanvas canvas;
	public MainWindow opts;
	public AdvancedOptions adv_opts;
	public Renderer renderer;
	public Controls controls;
	public SaveManager savemanager;
	public Preferences prefs;
	public MaterialManager materialmanager;
	public MaterialViewer materialviewer;
	
	public ArrayList<Plot> plots = new ArrayList<>();
	public BandPlot bandplot;
	public ScalarPlot scalarplot;
	public CarrierPlot carrierplot;
	public XYPlot xyplot;
	
	public Units units = Units.SI;
	public String description = "Description of simulation";
	
	/* Multithreading */
	
	public ReentrantReadWriteLock rwLock = new ReentrantReadWriteLock();
	public ReentrantLock poissonLock = new ReentrantLock(true);

	private CyclicBarrier start_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	private CyclicBarrier stop_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	private CyclicBarrier mid_barrier = new CyclicBarrier(SemiSim.n_threads);

	
	/* Domain parameters */

	public int default_resolution;
	public double default_width;
	
	public int resolution;
	public int nx;						// Number of grid points in x-dimension
	public int ny;						// Number of grid points in y-dimension
	public double width;				// Width, in SI
	public double ds;					// Spatial discretization
	public double dt;					// Timestep
	public double depth;				// Extent of circuit in z-dimension, only used to get reasonable values for current probe
	
	public int absorber_width;
	public double absorbing_coeff;
	public double boundary_stretch_factor;
	public double Hz_dissipation;
	public double dt_maximum;
	public int parity = -1;			// Negative sign resulting from flipped y-axis in graphics coordinate system

	public double time;
	public double AC_phase;
	public long stepnumber;
	public long frame;
	public int lastsimspeed = 0;
	public int iteration_multiplier = 0;
	
	public void resetTime() {
		time = 0;
		AC_phase = 0;
		stepnumber = 0;
		frame = 0;
		numerical_overflow = false;
		sign_violation_timer = 0;
	}

	public double error_detection_threshold = 1e-3;
	public int sign_violation_timer;
	public boolean numerical_overflow;
	public String CFL_text = "";
	
	public boolean advsettings_tweaked = false;
	
	public double AC_freq;
	public double AC_amplitude;
	public boolean AC_source_exists;

	
	/* Physical constants */

	public double eps0;			// Vaccuum permittivity
	public double mu0;			// Vaccum permeability
	public double c;			// Speed of light
	public double k;			// Boltzmann constant
	public double e_charge;		// Elementary charge
	public double e_mass;		// Electron mass
	public double eVtoJ = 1.6e-19;		// eV to J conversion factor
	public double q_n;			// Charge of single electron
	public double q_p;			// Charge of single hole

	
	/* Material properties */

	public double T;					// System temperature
	public double beta;					// Inverse temperature

	public double gv_semi;				// Semiconductor effective valence band DOS
	public double gc_semi;				// Semiconductor effective conduction DOS
	public double Eg_semi;				// Semiconductor band gap
	public double chi_semi;				// Semiconductor electron affinity

	public double mu_electron_semi;			// Electron mobility
	public double mu_hole_semi;				// Hole mobility, slightly lower than electron
	public double D_electron_semi;			// Diffusion constant, determined by Einstein relation
	public double D_hole_semi;				// Hole diffusion constant
	public double v_sat_n_semi;				// Velocity at which carrier velocity saturates
	public double v_sat_p_semi;				// Velocity at which carrier velocity saturates
	
	public double d_crit_n;				// Dopant density at which mobility drops by 1/2
	public double d_crit_p;
	public double eps_r_semi;			// Semiconductor dielectric constant

	public double k_rad_semi;			// Radiative recombination rate constant
	public double k_aug_n_semi;			// Auger recombination rate
	public double k_aug_p_semi;
	public double k_SRH_n_semi;			// Shockley-Read-Hall recombination rate, inverse electron lifetime
	public double k_SRH_p_semi;			// Inverse hole lifetime

	public double g_metal;				// Metal density of states
	public double g_metal_high;
	public double g_metal_low;
	public double W_metal_default;		// Metal work function, same as semiconductor
	public double W_metal_high;
	public double W_metal_low;
	public double Eg_metal;			// Same as semiconductor
	public double k_rad_metal;			// Radiative recombination rate constant

	public double g_currentsource;			// Current source carrier concentration
	public double currentsource_mobility;	// Current source mobility
	public double currentsource_sigma;	// Current source mobility
	
	public double n_default_doping_concentration;
	public double p_default_doping_concentration;
	public double n_light_doping_concentration;
	public double p_light_doping_concentration;
	public double n_heavy_doping_concentration;
	public double p_heavy_doping_concentration;

	public double dielectric_eps_r;
	public double ferromagnet_mu_r;
	public double staticcharge_density;
	public double max_EMF;
	public double max_current;
	public double default_AC_freq;

	public double default_flashlight_strength;
	
	public int junction_size;			// Free energy smoothing distance
	public int dopant_smoothing_distance;
	
	public double global_voltage_offset;

	public HashMap<MaterialType, String> default_names;
	public HashMap<MaterialType, String> modified_names;
	
	public boolean lock_resolution = false;

	public void setDefaultParameters() {
		if (!lock_resolution) {
			default_resolution = 256;
			default_width = 2.56e-5;
		}
		depth = 1e-3;

		eps0 = 8.85e-12;
		mu0 = 1.257e-6;
		k = 1.381e-23;
		e_charge = 1.6e-19;
		e_mass = 9e-31;
		T = 297.61;
		
		
		mu_electron_semi = 0.1400*1500;
		mu_hole_semi = 0.5*mu_electron_semi;
		
		v_sat_n_semi = 1.05e8;
		v_sat_p_semi = 0.525e8;

		gc_semi = 2.9366e25;
		gv_semi = 2.9366e25;
		chi_semi = 4.14*eVtoJ;
		Eg_semi = 1.12*eVtoJ;

		k_rad_semi = 1e-9;
		k_aug_n_semi = 0;
		k_aug_p_semi = 0;
		k_SRH_n_semi = 0;
		k_SRH_p_semi = 0;
		

		g_metal = 1.4683e30;
		g_metal_high = 2.5*g_metal;
		g_metal_low = 0.25*g_metal;
		W_metal_default = 4.7*eVtoJ;
		W_metal_high = W_metal_default + 0.3*eVtoJ;
		W_metal_low = W_metal_default - 0.3*eVtoJ;
		Eg_metal = 1.12*eVtoJ;
		k_rad_metal = 10*k_rad_semi;

		
		n_default_doping_concentration = 5e19;
		p_default_doping_concentration = 5e19;
		n_light_doping_concentration = 1e19;
		p_light_doping_concentration = 1e19;
		n_heavy_doping_concentration = 2.5e20;
		p_heavy_doping_concentration = 2.5e20;

		dielectric_eps_r = 25.0;
		ferromagnet_mu_r = 250.0;
		staticcharge_density = 10.0;
		currentsource_mobility = 0.0002;
		max_EMF = 5e5;
		max_current = 5e7;
		default_AC_freq = 1e13;
		
		default_flashlight_strength = 1e31;

		junction_size = 1;
		dopant_smoothing_distance = 0;
		
		d_crit_n = 1e100;
		d_crit_p = 1e100;
		eps_r_semi = 1;
		
		calculateDependentConstants();
	}
	
	public void calculateDependentConstants() {
		c = 1/Math.sqrt(eps0*mu0);
		q_n = -e_charge;
		q_p = e_charge;
		
		beta = 1/(k*T);
		D_electron_semi = mu_electron_semi/(beta*e_charge);
		D_hole_semi = mu_hole_semi/(beta*e_charge);

		g_currentsource = g_metal;
		
		global_voltage_offset = (chi_semi+0.5*Eg_semi)/eVtoJ;

		Material cur_source = new Material();
		initializeMaterial(cur_source, MaterialType.CURRENT);
		currentsource_sigma = cur_source.calcConductivity(e_charge, k*T);
		
		for (MaterialType type : MaterialType.values()) {
			String name = modified_names.get(type);
			if (name != null) type.name = name;
		}
	}
	
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

	public double[][] Jx_n;		// Electron current
	public double[][] Jy_n;
	public double[][] Jx_p;		// Hole current
	public double[][] Jy_p;
	public double[][] Jx_abs;		// Absorber current
	public double[][] Jy_abs;
	public double[][] Jx_free;		// Total free current
	public double[][] Jy_free;

	public double[][] Ey_xgrid;			// y component of E field on x grid
	public double[][] Ex_ygrid;			// x component of E field on y grid

	/* Miscellaneous fields used for display */

	public double[][] Dx;				// Electric displacement field
	public double[][] Dy;
	public double[][] Bz;				// B field
	public double[][] phi;				// Electric scalar potential

	public double[][] Sx;				// Electric displacement field
	public double[][] Sy;

	public double[][] G;				// Generation rate
	public double[][] R;				// Recombination rate
	public double[][] mu_n;				// Total chemical potential of electrons with electrostatic contribution, equal to quasi-fermi level
	public double[][] mu_p;				// Total chemical potential of holes, equal to negative of quasi-fermi level
	public double[][] grad_x_n;			// Gradient
	public double[][] grad_y_n;
	public double[][] grad_x_p;
	public double[][] grad_y_p;
	public double[][] V_avg;				// Average voltage, as measured by voltmeter

	public double[][] heat;
	public double[][] entropy;

	public double[][] vel_x_n;
	public double[][] vel_y_n;
	public double[][] vel_x_p;
	public double[][] vel_y_p;
	
	public double[][] drift_x_n;
	public double[][] drift_y_n;
	public double[][] drift_x_p;
	public double[][] drift_y_p;
	
	public double[][] diff_x_n;
	public double[][] diff_y_n;
	public double[][] diff_x_p;
	public double[][] diff_y_p;

	public double[][] sqrt_D_eff_n;
	public double[][] sqrt_D_eff_p;

	public boolean updateMiscFields = false;
	public boolean store_diff_drift = false;
	public boolean need_diff_drift = false;
	public double[][] debug;

	
	/* Material parameters */

	public Material[][] materials;

	public double[][] mu0_n;		// Chemical potential of a single electron
	public double[][] mu0_p;		// Chemical potential of a single hole
	public double[][] E0_n;			// Energy of a single electron
	public double[][] E0_p;			// Energy of a single hole
	public double[][] n_i;			// Square root of charge carrier equilibrium constant
	public double[][] k_rad;			// Radiative recombination rate constant
	public double[][] k_SRH_n;			// Shockley-Read-Hall recombination rate
	public double[][] k_SRH_p;
	public double[][] k_aug_n;			// Auger recombination rate
	public double[][] k_aug_p;
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

	public double[][] D_n;					// Electron diffusion constant
	public double[][] D_p;					// Hole diffusion constant
	public double[][] v_sat_n;				// Saturation velocity
	public double[][] v_sat_p;				// Saturation velocity
	public int[][] conducting;		// Does the material have partially filled bands
	public int[][] conducting_x;
	public int[][] conducting_y;

	public int[][] semiconducting;

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

	public int[][] abs_depth;
	public int[][] distance;
	public boolean[][] visited;
	public boolean[][] needs_smoothing;
	public double[][] smooth_arr;
	
	
	/* Probes */
	public List<Probe> probes = new CopyOnWriteArrayList<>();
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
		default_names = new HashMap<MaterialType, String>();
		for (MaterialType type : MaterialType.values()) {
			default_names.put(type, type.name);
		}
		
		controls = new Controls(this);
		prefs = new Preferences(this);
		renderer = new Renderer(this);
		savemanager = new SaveManager(this);
		canvas = renderer.new RenderCanvas(this);
		opts = new MainWindow(this);
		canvas.setFocusable(true);
		adv_opts = new AdvancedOptions(this);
		materialmanager = new MaterialManager(this);
		materialviewer = new MaterialViewer(this);
		datafile = new File(datafilename);

		bandplot = new BandPlot(); plots.add(bandplot);
		scalarplot = new ScalarPlot(); plots.add(scalarplot);
		carrierplot = new CarrierPlot(); plots.add(carrierplot);
		xyplot = new XYPlot(); plots.add(xyplot);
		plots.add(new ProbePlot("Voltage probe plot", "Voltage", "V", 1e-3, Quantity.ELECTRIC_POTENTIAL, (p) -> (p instanceof VoltageProbe && !(p instanceof Ground)), 0));
		plots.add(new ProbePlot("Current probe plot", "Current", "I", 1e-3, Quantity.ELECTRIC_CURRENT, (p) -> p instanceof CurrentProbe, 200));
		plots.add(new ProbePlot("Charge probe plot", "Charge", "Q", 1e-15, Quantity.CHARGE, (p) -> p instanceof ChargeProbe, 400));
		plots.add(new ProbePlot("Flux probe plot", "Magnetic flux", "\u03A6", 1e-15, Quantity.MAGNETIC_FLUX, (p) -> p instanceof FluxProbe, 600));
		
		for (Plot p : plots)
			p.initialize();

		SemiSim.detect64Bit();
		
		reset(true, Preset.DEFAULT);

		opts.initialize();
		adv_opts.initialize();
		prefs.initialize();
		materialmanager.initialize();
		materialviewer.initialize();

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
    					reset(false, null);
    				});
    				controls.clear = false;
    			}

    			if (controls.reset) {
    				SwingUtilities.invokeLater(() -> {
        				int result = JOptionPane.showConfirmDialog(opts, "Do you wish to reset the entire simulation?", "Message", JOptionPane.YES_NO_OPTION);
        				if (result == JOptionPane.OK_OPTION)
        				{
        					reset(true, null);
            				opts.setDefaults(this);
            				description = "Description of simulation";
            				SaveManager.currentfile = null;
        				}
    				});
    				controls.reset = false;
    			}

    			if (controls.updateimagesize) {

    				SwingUtilities.invokeLater(() -> {
    					rwLock.writeLock().lock();
    					try {
    						renderer.setCanvasSize();
    						//opts.pack();
    					}
    					finally {
    						rwLock.writeLock().unlock();
    					}
    				});
    				controls.updateimagesize = false;
    			}

    			lastsimspeed = opts.gui_simspeed.getValue();
    			dt = dt_maximum*(lastsimspeed/20.0);

    			if (controls.undo) {
    				SwingUtilities.invokeLater(() -> {
    					controls.undoredo.undo(this);
    				});
    				controls.undo = false;
    			}

    			if (controls.redo) {
    				SwingUtilities.invokeLater(() -> {
    					controls.undoredo.redo(this);
    				});
    				controls.redo = false;
    			}
    			
    			controls.handleMouseInput();

    			if (!opts.gui_paused.isSelected() || controls.advanceframe) {
    				for (int i = 0; i < iteration_multiplier ; i++) {
    					store_diff_drift = (i == iteration_multiplier-1 && need_diff_drift);
    					start_barrier.await();
    					stop_barrier.await();
    				}
    				calcMiscFields(false);
    				frame++;
    			} else if (updateMiscFields) {
    				calcMiscFields(false);
    			}
    			performProbeMeasurements();

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
    			if (!opts.gui_paused.isSelected())
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
	
	public void calculateMaxTimestep() {
		GeneralMaterialType[] types = materialmanager.makelist();
		
		double D_electron_max = 0;
		double D_hole_max = 0;
		double sigma_epsr_max = 0;
		
		String D_electron_name = "";
		String D_hole_name = "";
		String sigma_epsr_name = "";
		
		for (GeneralMaterialType type : types) {
			Material mat = new Material();
			initializeMaterial(mat, type);
			
			if (mat.D_n > D_electron_max) {
				D_electron_max = mat.D_n;
				D_electron_name = mat.getDisplayName();
			}
			
			if (mat.D_p > D_hole_max) {
				D_hole_max = mat.D_p;
				D_hole_name = mat.getDisplayName();
			}
			

			double rho_n = mat.calcEquilibriumElectronCharge(e_charge, k*T);
			double rho_p = mat.calcEquilibriumHoleCharge(e_charge, k*T);
			double sigma_epsr = e_charge*(-mat.D_n*beta*rho_n + mat.D_p*beta*rho_p)/mat.eps_r;
			if (sigma_epsr > sigma_epsr_max) {
				sigma_epsr_max = sigma_epsr;
				sigma_epsr_name = mat.getDisplayName();
			}
		}

		ds = width/default_resolution;
		dt_maximum = 0.9*Math.min(Math.min(ds/(Math.sqrt(2)*c), 4*eps0/sigma_epsr_max), Math.min(ds*ds/(4*D_electron_semi), ds*ds/(4*D_hole_semi)));
		
		Hz_dissipation = 0.01*ds*ds/dt_maximum;
		absorber_width = (int)Math.ceil(0.045*nx);
		absorbing_coeff = 50*c/(ds*nx);
		boundary_stretch_factor = 10;
		renderer.tau = 5000*dt_maximum;
		renderer.tau_events = renderer.tau*0.2;
		renderer.cc_default_dot_density = 1/(25*e_charge*n_default_doping_concentration*ds*ds);
		
		CFL_text = "";
		CFL_text += "Wave equation stability ratio = " + units.toString((Math.sqrt(2)*c)/(ds/dt_maximum), Quantity.DIMENSIONLESS) + "\n";
		CFL_text += "Electron diffusion stability ratio (" + D_electron_name + ") = " + units.toString(dt_maximum/(ds*ds/(4*D_electron_max)), Quantity.DIMENSIONLESS) + "\n";
		CFL_text += "Hole diffusion stability ratio (" + D_hole_name + ") = " + units.toString(dt_maximum/(ds*ds/(4*D_hole_max)), Quantity.DIMENSIONLESS) + "\n";
		CFL_text += "Conduction stability ratio (" + sigma_epsr_name + ") = " + units.toString(dt_maximum*sigma_epsr_max/(4*eps0), Quantity.DIMENSIONLESS) + "\n";
		System.out.print(CFL_text);
	}

	public boolean reset(boolean resetall, Preset preset) {
		
		if (resetall || preset != null) {
			modified_names = default_names;
			if (preset == null)
				Preset.DEFAULT.applyPreset(this);
			else
				preset.applyPreset(this);
			
			materialmanager.resetMaterialList();
		}

		boolean size_changed = (default_resolution != this.resolution);
		this.resolution = default_resolution;
		
		// Simulation and graphics threads should not be active when variables are initialized
		rwLock.writeLock().lock();
		try {
			if (size_changed) {
				System.out.println("Initializing simulation grid with resolution " + resolution);
				if(resolution < 4) {
					resolution = 4;
				}

				log2_resolution = (int) Math.round(Math.log(resolution)/Math.log(2));
				if(1 << log2_resolution != resolution) {
					resolution = 1 << log2_resolution;
				}

				nx = resolution;
				ny = resolution;
			}
			
			this.width = default_width;
			
			calculateMaxTimestep();

			if (size_changed) {
				Ex = new double[nx][ny];		Ey = new double[nx][ny];
				Bz = new double[nx][ny];
				Hz_laplacian = new double[nx][ny];
				rho_abs = new double[nx][ny];
				rho_n = new double[nx][ny];		rho_p = new double[nx][ny];
				rho_back = new double[nx][ny];
				rho_free = new double[nx][ny];

				Jx_abs = new double[nx][ny];	Jy_abs = new double[nx][ny];
				Jx_n = new double[nx][ny];		Jy_n = new double[nx][ny];
				Jx_p = new double[nx][ny];		Jy_p = new double[nx][ny];
				Jx_free = new double[nx][ny];	Jy_free = new double[nx][ny];

				Ey_xgrid = new double[nx][ny];	Ex_ygrid = new double[nx][ny];

				materials = new Material[nx][ny];
				n_i = new double[nx][ny];
				mu0_n = new double[nx][ny];		mu0_p = new double[nx][ny];
				mu_n = new double[nx][ny];		mu_p = new double[nx][ny];
				E0_n = new double[nx][ny];		E0_p = new double[nx][ny];
				k_rad = new double[nx][ny];
				k_SRH_n = new double[nx][ny];	k_SRH_p = new double[nx][ny];
				k_aug_n = new double[nx][ny];	k_aug_p = new double[nx][ny];
				L = new double[nx][ny];
				cmfx_n = new double[nx][ny];	cmfy_n = new double[nx][ny];
				cmfx_p = new double[nx][ny];	cmfy_p = new double[nx][ny];
				D_n = new double[nx][ny];		D_p = new double[nx][ny];
				v_sat_n = new double[nx][ny];	v_sat_p = new double[nx][ny];
				conducting = new int[nx][ny];
				semiconducting = new int[nx][ny];
				conducting_x = new int[nx][ny];	conducting_y = new int[nx][ny];
				ac_x = new int[nx][ny];			ac_y = new int[nx][ny];
				absorptivity = new double[nx][ny];
				absorptivity_x = new double[nx][ny]; absorptivity_y = new double[nx][ny];
				emfx = new double[nx][ny];		emfy = new double[nx][ny];
				epsx = new double[nx][ny];		epsy = new double[nx][ny];
				mu_z = new double[nx][ny];

				vel_x_n = new double[nx][ny];	vel_y_n = new double[nx][ny];
				vel_x_p = new double[nx][ny];	vel_y_p = new double[nx][ny];
				
				drift_x_n = new double[nx][ny];	drift_y_n = new double[nx][ny];
				drift_x_p = new double[nx][ny];	drift_y_p = new double[nx][ny];
				
				diff_x_n = new double[nx][ny];	diff_y_n = new double[nx][ny];
				diff_x_p = new double[nx][ny];	diff_y_p = new double[nx][ny];

				sqrt_D_eff_n = new double[nx][ny];		sqrt_D_eff_p = new double[nx][ny];
				
				heat = new double[nx][ny];
				entropy = new double[nx][ny];

				Dx = new double[nx][ny];		Dy = new double[nx][ny];
				Hz = new double[nx][ny];
				phi = new double[nx][ny];

				Sx = new double[nx][ny];		Sy = new double[nx][ny];

				G = new double[nx][ny];
				R = new double[nx][ny];
				mu_n = new double[nx][ny];			mu_p = new double[nx][ny];
				grad_x_n = new double[nx][ny];	grad_y_n = new double[nx][ny];
				grad_x_p = new double[nx][ny];	grad_y_p = new double[nx][ny];
				V_avg = new double[nx][ny];
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
				needs_smoothing = new boolean[nx][ny];
				abs_depth = new int[nx][ny];
				smooth_arr = new double[nx][ny];

				opts.setVisible(true);

				controls.setResolution(resolution);
				renderer.setResolution(resolution);
			}

			resetTime();
			
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{

					if (resetall || size_changed) {

						if (materials[i][j] == null)
							materials[i][j] = new Material();
						else
							materials[i][j].initialize();

						mu0_n[i][j] = 0;		mu0_p[i][j] = 0;
						mu_n[i][j] = 0;		mu_p[i][j] = 0;
						E0_n[i][j] = 0;		E0_p[i][j] = 0;
						n_i[i][j] = 0;
						k_rad[i][j] = 0;
						k_SRH_n[i][j] = 0;	k_SRH_p[i][j] = 0;
						k_aug_n[i][j] = 0;	k_aug_p[i][j] = 0;
						L[i][j] = 0;

						cmfx_n[i][j] = 0;	cmfy_n[i][j] = 0;
						cmfx_p[i][j] = 0;	cmfy_p[i][j] = 0;

						emfx[i][j] = 0.0;	emfy[i][j] = 0.0;

						epsx[i][j] = eps0;	epsy[i][j] = eps0;
						mu_z[i][j] = mu0;

						D_n[i][j] = 0.0;	D_p[i][j] = 0.0;
						v_sat_n[i][j] = 0.0;	v_sat_p[i][j] = 0.0;

						conducting[i][j] = 0;
						conducting_x[i][j] = 0;	conducting_y[i][j] = 0;
						
						semiconducting[i][j] = 0;

						ac_x[i][j] = 0; ac_y[i][j] = 0;

						absorptivity[i][j] = 0;
						absorptivity_x[i][j] = 0; absorptivity_y[i][j] = 0;

						controls.selected[i][j] = false;
						controls.selected_EMF[i][j] = false;
					}

					Ex [i][j] = 0.0;		Ey [i][j] = 0.0;
					Bz [i][j] = 0.0;
					Hz_laplacian [i][j] = 0.0;

					rho_abs[i][j] = 0.0;
					rho_n[i][j] = 0.0;		rho_p[i][j] = 0.0;
					rho_back[i][j] = 0.0;
					rho_free[i][j] = 0.0;

					Jx_abs[i][j] = 0.0;		Jy_abs[i][j] = 0.0;
					Jx_n[i][j] = 0.0;		Jy_n[i][j] = 0.0;
					Jx_p[i][j] = 0.0;		Jy_p[i][j] = 0.0;
					Jx_free[i][j] = 0.0;	Jy_free[i][j] = 0.0;

					Ey_xgrid[i][j] = 0.0;	Ex_ygrid[i][j] = 0.0;

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
					needs_smoothing[i][j] = false;

					drift_x_n[i][j] = 0;	drift_y_n[i][j] = 0;
					drift_x_p[i][j] = 0;	drift_y_p[i][j] = 0;
					
					diff_x_n[i][j] = 0;	diff_y_n[i][j] = 0;
					diff_x_p[i][j] = 0;	diff_y_p[i][j] = 0;
					
					heat[i][j] = 0.0;
					entropy[i][j] = 0.0;
					
					Dx[i][j] = 0.0; Dy[i][j] = 0.0;
					Hz[i][j] = 0.0;
					phi[i][j] = 0.0;

					Sx[i][j] = 0.0; Sy[i][j] = 0.0;

					G[i][j] = 0.0;
					R[i][j] = 0.0;
					mu_n[i][j] = 0.0;
					mu_p[i][j] = 0.0;
					grad_x_n[i][j] = 0.0;	grad_y_n[i][j] = 0.0;
					grad_x_p[i][j] = 0.0;	grad_y_p[i][j] = 0.0;
					grad_x_n[i][j] = 0.0;	grad_y_n[i][j] = 0.0;
					grad_x_p[i][j] = 0.0;	grad_y_p[i][j] = 0.0;
					V_avg[i][j] = 0.0;
					
					abs_depth[i][j] = 0;
					
					debug[i][j] = 0.0;
				}
			}

			controls.selection.clear();
			controls.clipboard.clear();
			
			if (resetall || size_changed) {
				controls.EMF_selected = false;
				controls.changesmade = false;
				opts.gui_bc.setSelectedItem(BoundaryCondition.DISSIPATIVE);
				controls.prev_boundary = BoundaryCondition.DISSIPATIVE;
				
				probes.clear();
				probe_index = 0;

				for (Plot p: plots) {
					p.frame.setVisible(false);
				}

				controls.undoredo.resetUndoHistory(this);
				controls.resetZoom();
			}
			
			for (Probe p : probes) {
				p.reset();
				p.data.resetData();
			}
			
			renderer.resetChargeDots();

			initializeAllMaterials();
			updateAllMaterials(true);
			controls.undoredo.captureState(this);
			multigridSolve(true, false);
		} finally {
			rwLock.writeLock().unlock();
		}
		
		return true;
	}

	public void constructBoundary() {
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
							initializeMaterial(materials[i][j], MaterialType.ABSORBER);
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

		double coeff = Math.log(boundary_stretch_factor);
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
			System.out.println("Simulation thread " + n_thread + ": " + i_min + " < i <= " + i_max + " initialized.");
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
					
					for (int i = 0; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-1; j++)
							{
								Ey_xgrid[i][j] = 0.25*(Ey[i][j] + Ey[i+1][j] + Ey[i][j-1] + Ey[i+1][j-1]);
							}
						}
					}
					

					for (int i = 1; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 0; j < ny-1; j++)
							{
								Ex_ygrid[i][j] = 0.25*(Ex[i][j] + Ex[i][j+1] + Ex[i-1][j] + Ex[i-1][j+1]);
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

								double sigma_n = 0;
								double sigma_p = 0;

								if (conducting_x[i][j] == 1) {

									double d_n = 0.5*(D_n[i+1][j] + D_n[i][j]);
									double d_p = 0.5*(D_p[i+1][j] + D_p[i][j]);
									
									double rho_n_avg = 0.5*(rho_n[i+1][j] + rho_n[i][j]);
									double rho_p_avg = 0.5*(rho_p[i+1][j] + rho_p[i][j]);
									
									double Esat_n = 0.5*(v_sat_n[i+1][j] + v_sat_n[i][j])/(d_n*e_charge*beta);
									double Esat_p = 0.5*(v_sat_p[i+1][j] + v_sat_p[i][j])/(d_p*e_charge*beta);

									double mr_n = 1/Math.sqrt((ex_prev*ex_prev + Ey_xgrid[i][j]*Ey_xgrid[i][j])/(Esat_n*Esat_n) + 1);
									double mr_p = 1/Math.sqrt((ex_prev*ex_prev + Ey_xgrid[i][j]*Ey_xgrid[i][j])/(Esat_p*Esat_p) + 1);

									double emf_phase = ac_x[i][j]*AC_amplitude+(1-ac_x[i][j]);

									double w_n = (emf_phase*emfx[i][j]*q_n + cmfx_n[i][j] + ex_prev*q_n)*beta*ds/2;
									double w_p = (emf_phase*emfx[i][j]*q_p + cmfx_p[i][j] + ex_prev*q_p)*beta*ds/2;
									
									// whether to use taylor series approximation around x=0
									boolean approx_n = Math.abs(w_n) < 0.2;
									boolean approx_p = Math.abs(w_p) < 0.2;

									double exp_2wn = approx_n? 1 : FastExp.exp(2*w_n);
									double exp_2wp = approx_p? 1 : FastExp.exp(2*w_p);

									sigma_n = -d_n*rho_n_avg*e_charge*beta*(1 + Ey_xgrid[i][j]*Ey_xgrid[i][j]/(Esat_n*Esat_n))*(mr_n*mr_n*mr_n);
									sigma_p = d_p*rho_p_avg*e_charge*beta*(1 + Ey_xgrid[i][j]*Ey_xgrid[i][j]/(Esat_p*Esat_p))*(mr_p*mr_p*mr_p);

									Jx_n[i][j] = mr_n*d_n/ds*(-(rho_n[i+1][j] - rho_n[i][j])*Utils.xtanhxm1(w_n, exp_2wn, approx_n)
									+ 2*rho_n_avg*w_n) - sigma_n*ex_prev;
									Jx_p[i][j] = mr_p*d_p/ds*(-(rho_p[i+1][j] - rho_p[i][j])*Utils.xtanhxm1(w_p, exp_2wp, approx_p)
									+ 2*rho_p_avg*w_p) - sigma_p*ex_prev;
									
									if (store_diff_drift) {
										diff_x_n[i][j] = mr_n*d_n/ds*(-(rho_n[i+1][j] - rho_n[i][j])*Utils.xtanhxm1(w_n, exp_2wn, approx_n));
										diff_x_p[i][j] = mr_p*d_p/ds*(-(rho_p[i+1][j] - rho_p[i][j])*Utils.xtanhxm1(w_p, exp_2wp, approx_p));

										drift_x_n[i][j] = 2*rho_n_avg*w_n*mr_n*d_n/ds;
										drift_x_p[i][j] = 2*rho_p_avg*w_p*mr_p*d_p/ds;

										vel_x_n[i][j] = 2*w_n*mr_n*d_n/ds;
										vel_x_p[i][j] = 2*w_p*mr_p*d_p/ds;
									}
								} else {
									Jx_n[i][j] = 0;
									Jx_p[i][j] = 0;

									if (store_diff_drift) {
										diff_x_n[i][j] = 0;
										diff_x_p[i][j] = 0;

										drift_x_n[i][j] = 0;
										drift_x_p[i][j] = 0;

										vel_x_n[i][j] = 0;
										vel_x_p[i][j] = 0;
									}
								}

								Jx_abs[i][j] = 0;
								
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
								
								double sigma_n = 0;
								double sigma_p = 0;

								if (conducting_y[i][j] == 1) {

									double d_n = 0.5*(D_n[i][j+1] + D_n[i][j]);
									double d_p = 0.5*(D_p[i][j+1] + D_p[i][j]);
									
									double rho_n_avg = 0.5*(rho_n[i][j+1] + rho_n[i][j]);
									double rho_p_avg = 0.5*(rho_p[i][j+1] + rho_p[i][j]);
									
									double Esat_n = 0.5*(v_sat_n[i][j+1] + v_sat_n[i][j])/(d_n*e_charge*beta);
									double Esat_p = 0.5*(v_sat_p[i][j+1] + v_sat_p[i][j])/(d_p*e_charge*beta);

									double mr_n = 1/Math.sqrt((ey_prev*ey_prev + Ex_ygrid[i][j]*Ex_ygrid[i][j])/(Esat_n*Esat_n)+1);
									double mr_p = 1/Math.sqrt((ey_prev*ey_prev + Ex_ygrid[i][j]*Ex_ygrid[i][j])/(Esat_p*Esat_p)+1);

									double emf_phase = ac_y[i][j]*AC_amplitude+(1-ac_y[i][j]);

									double w_n = (emf_phase*emfy[i][j]*q_n + cmfy_n[i][j] + ey_prev*q_n)*beta*ds/2;
									double w_p = (emf_phase*emfy[i][j]*q_p + cmfy_p[i][j] + ey_prev*q_p)*beta*ds/2;

									boolean approx_n = Math.abs(w_n) < 0.2;
									boolean approx_p = Math.abs(w_p) < 0.2;

									double exp_2wn = approx_n? 1 : FastExp.exp(2*w_n);
									double exp_2wp = approx_p? 1 : FastExp.exp(2*w_p);

									sigma_n = -d_n*rho_n_avg*e_charge*beta*(1 + Ex_ygrid[i][j]*Ex_ygrid[i][j]/(Esat_n*Esat_n))*(mr_n*mr_n*mr_n);
									sigma_p = d_p*rho_p_avg*e_charge*beta*(1 + Ex_ygrid[i][j]*Ex_ygrid[i][j]/(Esat_p*Esat_p))*(mr_p*mr_p*mr_p);

									Jy_n[i][j] = mr_n*d_n/ds*(-(rho_n[i][j+1] - rho_n[i][j])*Utils.xtanhxm1(w_n, exp_2wn, approx_n)
										+ 2*rho_n_avg*w_n) - sigma_n*ey_prev;
									Jy_p[i][j] = mr_p*d_p/ds*(-(rho_p[i][j+1] - rho_p[i][j])*Utils.xtanhxm1(w_p, exp_2wp, approx_p)
										+ 2*rho_p_avg*w_p) - sigma_p*ey_prev;
									
									if (store_diff_drift) {
										diff_y_n[i][j] = mr_n*d_n/ds*(-(rho_n[i][j+1] - rho_n[i][j])*Utils.xtanhxm1(w_n, exp_2wn, approx_n));
										diff_y_p[i][j] = mr_p*d_p/ds*(-(rho_p[i][j+1] - rho_p[i][j])*Utils.xtanhxm1(w_p, exp_2wp, approx_p));
										
										drift_y_n[i][j] = 2*rho_n_avg*w_n*mr_n*d_n/ds;
										drift_y_p[i][j] = 2*rho_p_avg*w_p*mr_p*d_p/ds;

										vel_y_n[i][j] = 2*w_n*mr_n*d_n/ds;
										vel_y_p[i][j] = 2*w_p*mr_p*d_p/ds;
									}

								} else {
									Jy_n[i][j] = 0;
									Jy_p[i][j] = 0;
									
									if (store_diff_drift) {
										diff_y_n[i][j] = 0;
										diff_y_p[i][j] = 0;
										
										drift_y_n[i][j] = 0;
										drift_y_p[i][j] = 0;
										
										vel_y_n[i][j] = 0;
										vel_y_p[i][j] = 0;
									}
								}

								Jy_abs[i][j] = 0;

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
								double recomb_rate = 0;
								if (conducting[i][j] == 1) {
									double n = rho_n[i][j]/q_n;
									double p = rho_p[i][j]/q_p;
									double ni = n_i[i][j];
									double rate_const = (k_aug_n[i][j]*n+k_aug_p[i][j]*p)
										+ (k_SRH_n[i][j]*k_SRH_p[i][j])/(k_SRH_n[i][j]*(n+ni) + k_SRH_p[i][j]*(p+ni) + Double.MIN_VALUE) // Avoid divide by zero
										+ k_rad[i][j];
									recomb_rate = rate_const*(n*p - ni*ni) - L[i][j];
								}

								rho_n[i][j] = rho_n[i][j] - (Jx_n[i][j]-Jx_n[i-1][j] + Jy_n[i][j]-Jy_n[i][j-1])*dt/ds - dt*q_n*recomb_rate;
								rho_p[i][j] = rho_p[i][j] - (Jx_p[i][j]-Jx_p[i-1][j] + Jy_p[i][j]-Jy_p[i][j-1])*dt/ds - dt*q_p*recomb_rate;
								rho_abs[i][j] = rho_abs[i][j] - (Jx_abs[i][j]-Jx_abs[i-1][j] + Jy_abs[i][j]-Jy_abs[i][j-1])*dt/ds;

								rho_free[i][j] = rho_abs[i][j]+rho_n[i][j]+rho_p[i][j]+rho_back[i][j];
							}
						}
					}

					mid_barrier.await();
					
					if (n_thread == 0) {
						t6.stop();
						
						/* Enforce Gauss law constraint */
						if (stepnumber%500 == 0) {
							multigridSolve(true, false);
						}

						time += dt;
						AC_phase += 2*Math.PI*AC_freq*dt;
						AC_amplitude = Math.cos(AC_phase);

						controls.advanceframe = false;
					}

					stop_barrier.await();
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void calcMiscFields(boolean updatePhi) {
		ScalarView view_scalar = controls.scalarview.getOption();
		VectorView view_vector = controls.vectorview.getOption();

		if (updatePhi)
			multigridSolve(false, true);

		t8.start();

		//sign_violation = 0;
		double rho_n_max_tmp = 0;
		double rho_p_max_tmp = 0;
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				// Prevent errors when taking log of 0
				if (conducting[i][j] == 0) {
					mu_n[i][j] = Double.NaN;
					mu_p[i][j] = Double.NaN;
					V_avg[i][j] = Double.NaN;
				} else {
					mu_n[i][j] = mu0_n[i][j] + k*T*Math.log(rho_n[i][j]/q_n) + phi[i][j]*q_n;
					mu_p[i][j] = mu0_p[i][j] + k*T*Math.log(rho_p[i][j]/q_p) + phi[i][j]*q_p;
					// Think about why we should take average


					//double probe_n = -e_charge*ni_metal;
					//double probe_p = e_charge*ni_metal;
					//double sigma_n = mu_electron*logmean(-rho_n[i][j], -probe_n);
					//double sigma_p = mu_hole*logmean(rho_p[i][j], probe_p);
					double sigma_n = D_n[i][j]*beta*e_charge*(-rho_n[i][j]);
					double sigma_p = D_p[i][j]*beta*e_charge*rho_p[i][j];

					V_avg[i][j] = (sigma_n*mu_n[i][j]/q_n + sigma_p*mu_p[i][j]/q_p)
							/(sigma_n + sigma_p) - global_voltage_offset;
				}

				if (conducting[i][j] == 1) {
					double n = rho_n[i][j]/q_n;
					double p = rho_p[i][j]/q_p;
					double ni = n_i[i][j];
					double rate_const = (k_aug_n[i][j]*n+k_aug_p[i][j]*p)
						+ (k_SRH_n[i][j]*k_SRH_p[i][j])/(k_SRH_n[i][j]*(n+ni) + k_SRH_p[i][j]*(p+ni) + Double.MIN_VALUE) // Avoid divide by zero
						+ k_rad[i][j];
					
					G[i][j] = rate_const*ni*ni + L[i][j];
					R[i][j] = rate_const*n*p;
				} else {
					G[i][j] = 0;
					R[i][j] = 0;
				}

				if (rho_n[i][j] > error_detection_threshold || rho_p[i][j] < -error_detection_threshold) {
					sign_violation_timer = 20;
					debug[i][j] = 1;
				} else {
					debug[i][j] = 0;
				}

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
				Sx[i][j] = -Ex[i][j]*0.5*(Hz[i][j] + Hz[i][j-1]);
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				Jy_free[i][j] = Jy_n[i][j] + Jy_p[i][j];
				Dy[i][j] = epsy[i][j]*Ey[i][j];
				Sy[i][j] = Ey[i][j]*0.5*(Hz[i][j] + Hz[i-1][j]);
			}
		}

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				Bz[i][j] = Hz[i][j]*mu_z[i][j];
				
				if (!Double.isFinite(Hz[i][j]))
					numerical_overflow = true;
			}
		}

		if (view_scalar == ScalarView.HEAT) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					grad_x_n[i][j] = conducting_x[i][j] == 1? -conducting_x[i][j]*(E0_n[i+1][j] - E0_n[i][j])/ds : 0;
					grad_x_p[i][j] = conducting_x[i][j] == 1? -conducting_x[i][j]*(E0_p[i+1][j] - E0_p[i][j])/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					grad_y_n[i][j] = conducting_x[i][j] == 1? -conducting_y[i][j]*(E0_n[i][j+1] - E0_n[i][j])/ds : 0;
					grad_y_p[i][j] = conducting_x[i][j] == 1? -conducting_y[i][j]*(E0_p[i][j+1] - E0_p[i][j])/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					double n_contrib = -(G[i][j]-R[i][j])*E0_n[i][j];
					double jn_contrib = 0.5*(grad_x_n[i-1][j]*Jx_n[i-1][j]+grad_x_n[i][j]*Jx_n[i][j]
							+grad_y_n[i][j-1]*Jy_n[i][j-1]+grad_y_n[i][j]*Jy_n[i][j])/q_n;

					double p_contrib = -(G[i][j]-R[i][j])*E0_p[i][j];
					double jp_contrib = 0.5*(grad_x_p[i-1][j]*Jx_p[i-1][j]+grad_x_p[i][j]*Jx_p[i][j]
							+grad_y_p[i][j-1]*Jy_p[i][j-1]+grad_y_p[i][j]*Jy_p[i][j])/q_p;

					double ohm_contrib = 0.5*(Ex[i-1][j]*Jx_free[i-1][j]+Ex[i][j]*Jx_free[i][j]+Ey[i][j-1]*Jy_free[i][j-1]+Ey[i][j]*Jy_free[i][j]);

					heat[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
				}
			}
		}

		if (view_scalar == ScalarView.ENTROPY) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					grad_x_n[i][j] = conducting_x[i][j] == 1? -(mu_n[i+1][j] - mu_n[i][j] - phi[i+1][j]*q_n + phi[i][j]*q_n)/ds : 0;
					grad_x_p[i][j] = conducting_x[i][j] == 1? -(mu_p[i+1][j] - mu_p[i][j] - phi[i+1][j]*q_p + phi[i][j]*q_p)/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					grad_y_n[i][j] = conducting_y[i][j] == 1? -(mu_n[i][j+1] - mu_n[i][j] - phi[i][j+1]*q_n + phi[i][j]*q_n)/ds : 0;
					grad_y_p[i][j] = conducting_y[i][j] == 1? -(mu_p[i][j+1] - mu_p[i][j] - phi[i][j+1]*q_p + phi[i][j]*q_p)/ds : 0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					double n_contrib = -(G[i][j]-R[i][j])*mu_n[i][j];
					double jn_contrib = 0.5*(grad_x_n[i-1][j]*Jx_n[i-1][j]+grad_x_n[i][j]*Jx_n[i][j]
							+grad_y_n[i][j-1]*Jy_n[i][j-1]+grad_y_n[i][j]*Jy_n[i][j])/q_n;

					double p_contrib = -(G[i][j]-R[i][j])*mu_p[i][j];
					double jp_contrib = 0.5*(grad_x_p[i-1][j]*Jx_p[i-1][j]+grad_x_p[i][j]*Jx_p[i][j]
							+grad_y_p[i][j-1]*Jy_p[i][j-1]+grad_y_p[i][j]*Jy_p[i][j])/q_p;

					double ohm_contrib = 0.5*(Ex[i-1][j]*Jx_free[i-1][j]+Ex[i][j]*Jx_free[i][j]+Ey[i][j-1]*Jy_free[i][j-1]+Ey[i][j]*Jy_free[i][j]);

					entropy[i][j] = (n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib)/T;
				}
			}
		}
		
		if (view_scalar == ScalarView.ELECTRON_VEL || view_scalar == ScalarView.HOLE_VEL
			|| view_vector == VectorView.ELECTRON_DIFFUSION || view_vector == VectorView.ELECTRON_DRIFT || view_vector == VectorView.ELECTRON_VELOCITY
			|| view_vector == VectorView.HOLE_DIFFUSION || view_vector == VectorView.HOLE_DRIFT || view_vector == VectorView.HOLE_VELOCITY
			|| (opts.gui_carriers.isSelected() && opts.menu_carrier_diffusion.isSelected()))
		{
			need_diff_drift = true;
			
			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					if (conducting[i][j] == 1) {
						double E2 = 0.25*(Ex[i][j]+Ex[i][j+1])*(Ex[i][j]+Ex[i][j+1]) + 0.25*(Ey[i][j]+Ey[i+1][j])*(Ey[i][j]+Ey[i+1][j]);
						double Esat_n = v_sat_n[i][j]/(D_n[i][j]*e_charge*beta);
						double Esat_p = v_sat_p[i][j]/(D_p[i][j]*e_charge*beta);
						sqrt_D_eff_n[i][j] = Math.sqrt(D_n[i][j]/Math.sqrt(E2/(Esat_n*Esat_n) + 1));
						sqrt_D_eff_p[i][j] = Math.sqrt(D_p[i][j]/Math.sqrt(E2/(Esat_p*Esat_p) + 1));
					} else {
						sqrt_D_eff_n[i][j] = 0;
						sqrt_D_eff_p[i][j] = 0;
					}
				}
			}
		} else {
			need_diff_drift = false;
		}
		
		updateMiscFields = false;

		t8.stop();
	}
	
	public void performProbeMeasurements() {
		for (Probe p: probes) {
			p.measure(this, frame%controls.plotinterval == 0);
		}

		if (controls.logdata) {
			logProbeData();
			controls.logdata = false;
		}
	}
	
	public void computeScalarField(double[][] scalarfield, int i1, int j1, ScalarView scalarview) {
		int sx = scalarfield.length;
		int sy = scalarfield[0].length;
		
		switch (scalarview) {
		case NONE:
			break;
		case B_FIELD:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = parity*0.25*(Bz[i+i1][j+j1]+Bz[i+i1-1][j+j1]+Bz[i+i1][j+j1-1]+Bz[i+i1-1][j+j1-1]);
				}
			}
			break;
		case E_FIELD:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = Utils.length(0.5*(Ex[i+i1][j+j1]+Ex[i+i1-1][j+j1]), 0.5*(Ey[i+i1][j+j1]+Ey[i+i1][j+j1-1]));
				}
			}
			break;
		case H_FIELD:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = parity*0.25*(Hz[i+i1][j+j1]+Hz[i+i1-1][j+j1]+Hz[i+i1][j+j1-1]+Hz[i+i1-1][j+j1-1]);
				}
			}
			break;
		case CURRENT:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = Utils.length(0.5*(Jx_free[i+i1][j+j1]+Jx_free[i+i1-1][j+j1]), 0.5*(Jy_free[i+i1][j+j1]+Jy_free[i+i1][j+j1-1]));
				}
			}
			break;
		case POTENTIAL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = phi[i+i1][j+j1];
				}
			}
			break;
		case CHARGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_free[i+i1][j+j1];
				}
			}
			break;
		case BACKGROUND_CHARGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_back[i+i1][j+j1];
				}
			}
			break;
		case ELECTRON_CHARGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_n[i+i1][j+j1];
				}
			}
			break;
		case HOLE_CHARGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_p[i+i1][j+j1];
				}
			}
			break;
		case COMBINED_CHARGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = Double.NaN;
				}
			}
			break;
		case HEAT:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = heat[i+i1][j+j1];
				}
			}
			break;
		case ELECTRON_POTENTIAL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = mu_n[i+i1][j+j1]/eVtoJ;
				}
			}
			break;
		case HOLE_POTENTIAL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = -mu_p[i+i1][j+j1]/eVtoJ;
				}
			}
			break;
		case DEBUG:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = debug[i+i1][j+j1];
				}
			}
			break;
		case GENERATION:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = G[i+i1][j+j1];
				}
			}
			break;
		case RECOMBINATION:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = R[i+i1][j+j1];
				}
			}
			break;
		case AVERAGE_POTENTIAL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = V_avg[i+i1][j+j1]; // -offset
				}
			}
			break;
		case ELECTRON_DENSITY:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_n[i+i1][j+j1]/q_n;
				}
			}
			break;
		case ELECTRON_VEL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = Utils.length(0.5*(vel_x_n[i+i1][j+j1]+vel_x_n[i+i1-1][j+j1]), 0.5*(vel_y_n[i+i1][j+j1]+vel_y_n[i+i1][j+j1-1]));
				}
			}
			break;
		case ELECTRON_VOLTAGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = -mu_n[i+i1][j+j1]/eVtoJ;
				}
			}
			break;
		case ENERGY:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1) {
						double E_energy = 0.25*(Ex[i+i1][j+j1]*Ex[i+i1][j+j1]*epsx[i+i1][j+j1] + Ex[i+i1-1][j+j1]*Ex[i+i1-1][j+j1]*epsx[i+i1-1][j+j1]
								+ Ey[i+i1][j+j1]*Ey[i+i1][j+j1]*epsy[i+i1][j+j1] + Ey[i+i1][j+j1-1]*Ey[i+i1][j+j1-1]*epsy[i+i1][j+j1-1]);
						double B_energy = 0.25*(Hz[i+i1][j+j1]*Hz[i+i1][j+j1]*mu_z[i+i1][j+j1] + Hz[i+i1-1][j+j1]*Hz[i+i1-1][j+j1]*mu_z[i+i1-1][j+j1]
								+ Hz[i+i1][j+j1-1]*Hz[i+i1][j+j1-1]*mu_z[i+i1][j+j1-1] + Hz[i+i1-1][j+j1-1]*Hz[i+i1-1][j+j1-1]*mu_z[i+i1-1][j+j1-1]);
						scalarfield[i][j] = E_energy + B_energy;
					}
				}
			}
			break;
		case ENTROPY:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = entropy[i+i1][j+j1];
				}
			}
			break;
		case HOLE_DENSITY:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = rho_p[i+i1][j+j1]/q_p;
				}
			}
			break;
		case HOLE_VEL:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					if (i + i1 > 0 && j + j1 > 0 && i + i1 < nx-1 && j + j1 < ny-1)
						scalarfield[i][j] = Utils.length(0.5*(vel_x_p[i+i1][j+j1]+vel_x_p[i+i1-1][j+j1]), 0.5*(vel_y_p[i+i1][j+j1]+vel_y_p[i+i1][j+j1-1]));
				}
			}
			break;
		case HOLE_VOLTAGE:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					scalarfield[i][j] = mu_p[i+i1][j+j1]/eVtoJ;
				}
			}
			break;
		case LIGHT:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					double n = rho_n[i+i1][j+j1]/q_n;
					double p = rho_p[i+i1][j+j1]/q_p;
					scalarfield[i][j] = semiconducting[i+i1][j+j1]*k_rad[i+i1][j+j1]*n*p;
					if (!(scalarfield[i][j] > 0))
						scalarfield[i][j] = 0;
				}
			}
			break;
		case RECOMB_AUGER:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					double n = rho_n[i+i1][j+j1]/q_n;
					double p = rho_p[i+i1][j+j1]/q_p;
					scalarfield[i][j] = (k_aug_n[i+i1][j+j1]*n+k_aug_p[i+i1][j+j1]*p)*n*p;
				}
			}
			break;
		case RECOMB_RAD:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					double n = rho_n[i+i1][j+j1]/q_n;
					double p = rho_p[i+i1][j+j1]/q_p;
					scalarfield[i][j] = k_rad[i+i1][j+j1]*n*p;
				}
			}
			break;
		case RECOMB_SRH:
			for (int i = 0; i < sx; i++) {
				for (int j = 0; j < sy; j++) {
					double n = rho_n[i+i1][j+j1]/q_n;
					double p = rho_p[i+i1][j+j1]/q_p;
					double ni = n_i[i+i1][j+j1];
					double rate_const = (k_SRH_n[i+i1][j+j1]*k_SRH_p[i+i1][j+j1])/(k_SRH_n[i+i1][j+j1]*(n+ni) + k_SRH_p[i+i1][j+j1]*(p+ni) + Double.MIN_VALUE);
					scalarfield[i][j] = rate_const*n*p;
				}
			}
			break;
		default:
			break;
		}
	}
	
	public void computeVectorField(double[][][] vectorfield, VectorView vector_view) {
		switch (vector_view) {
		case NONE:
			break;
		case D_FIELD:
			vectorfield[0] = Dx;
			vectorfield[1] = Dy;
			break;
		case E_FIELD:
			vectorfield[0] = Ex;
			vectorfield[1] = Ey;
			break;
		case ELECTRON_CURRENT:
			vectorfield[0] = Jx_n;
			vectorfield[1] = Jy_n;
			break;
		case HOLE_CURRENT:
			vectorfield[0] = Jx_p;
			vectorfield[1] = Jy_p;
			break;
		case TOTAL_CURRENT:
			vectorfield[0] = Jx_free;
			vectorfield[1] = Jy_free;
			break;
		case EMF:
			vectorfield[0] = emfx;
			vectorfield[1] = emfy;
			break;
		case ELECTRON_DIFFUSION:
			vectorfield[0] = diff_x_n;
			vectorfield[1] = diff_y_n;
			break;
		case ELECTRON_DRIFT:
			vectorfield[0] = drift_x_n;
			vectorfield[1] = drift_y_n;
			break;
		case HOLE_DIFFUSION:
			vectorfield[0] = diff_x_p;
			vectorfield[1] = diff_y_p;
			break;
		case HOLE_DRIFT:
			vectorfield[0] = drift_x_p;
			vectorfield[1] = drift_y_p;
			break;
		case POYNTING:
			vectorfield[0] = Sx;
			vectorfield[1] = Sy;
			break;
		case ELECTRON_VELOCITY:
			vectorfield[0] = vel_x_n;
			vectorfield[1] = vel_y_n;
			break;
		case HOLE_VELOCITY:
			vectorfield[0] = vel_x_p;
			vectorfield[1] = vel_y_p;
			break;
		}
	}

	public void computeChemicalForces()
	{

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 1) {
					E0_n[i][j] = materials[i][j].Ec;
					E0_p[i][j] = -materials[i][j].Ev;
					mu0_n[i][j] = materials[i][j].Ec - k*T*Math.log(materials[i][j].gc);
					mu0_p[i][j] = -materials[i][j].Ev - k*T*Math.log(materials[i][j].gv);
				} else {
					E0_n[i][j] = 0;
					E0_p[i][j] = 0;
					mu0_n[i][j] = 0;
					mu0_p[i][j] = 0;
				}
			}
		}

		// Smoothing out free energy in space makes simulation more stable (No longer required in v2.0)
		if (junction_size > 0) {
			markJunctions(junction_size);
			smoothArray(mu0_n, junction_size, false);
			smoothArray(mu0_p, junction_size, false);
			smoothArray(E0_n, junction_size, false);
			smoothArray(E0_p, junction_size, false);
			smoothArray(k_rad, junction_size, true);
			smoothArray(k_SRH_n, junction_size, true);
			smoothArray(k_SRH_p, junction_size, true);
			smoothArray(k_aug_n, junction_size, true);
			smoothArray(k_aug_p, junction_size, true);
		}
		
		if (dopant_smoothing_distance > 0) {
			markJunctions(dopant_smoothing_distance);
			smoothArray(rho_back, dopant_smoothing_distance, false);
		}

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 1) {
					n_i[i][j] = Math.exp(-0.5*beta*(mu0_n[i][j]+mu0_p[i][j]));
				} else {
					n_i[i][j] = 0;
				}
			}
		}

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				cmfx_n[i][j] = -conducting_x[i][j]*(mu0_n[i+1][j] - mu0_n[i][j])/ds;
				cmfx_p[i][j] = -conducting_x[i][j]*(mu0_p[i+1][j] - mu0_p[i][j])/ds;
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				cmfy_n[i][j] = -conducting_y[i][j]*(mu0_n[i][j+1] - mu0_n[i][j])/ds;
				cmfy_p[i][j] = -conducting_y[i][j]*(mu0_p[i][j+1] - mu0_p[i][j])/ds;
			}
		}
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (conducting[i][j] == 0) {
					mu0_n[i][j] = Double.NaN;
					mu0_p[i][j] = Double.NaN;
					E0_n[i][j] = Double.NaN;
					E0_p[i][j] = Double.NaN;
				}
			}
		}
	}
	
	public void smoothArray(double[][] arr, int smoothing_radius, boolean logarithmic) {

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				smooth_arr[i][j] = logarithmic? Math.log(arr[i][j]) : arr[i][j];
			}
		}

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (needs_smoothing[i][j]) {
					double sum = 0;
					double neighbors = 0;

					markNeighborhood(i, j, smoothing_radius);

					for (int di = -smoothing_radius; di <= smoothing_radius; di++) {
						for (int dj = -smoothing_radius; dj <= smoothing_radius; dj++) {
							if (i+di >= 0 && j+dj >= 0 && i+di < nx && j+dj < ny && conducting[i+di][j+dj] == 1 && visited[i+di][j+dj]) {
								sum += smooth_arr[i+di][j+dj];
								neighbors += 1;
							}
							distance[i+di][j+dj] = Integer.MAX_VALUE;
							visited[i+di][j+dj] = false;
						}
					}

					if (logarithmic) {
						arr[i][j] = (sum == Double.NaN)? 0 : Math.exp(sum/neighbors);
					} else {
						arr[i][j] = sum/neighbors;
					}
				}
			}
		}
	}

	public void updateAllMaterials(boolean updateRho) {
		constructBoundary();
		
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
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
				conducting[i][j] = materials[i][j].conducting;
				semiconducting[i][j] = materials[i][j].semiconducting;
				k_rad[i][j] = materials[i][j].k_rad;
				k_SRH_n[i][j] = materials[i][j].k_SRH_n;
				k_SRH_p[i][j] = materials[i][j].k_SRH_p;
				k_aug_n[i][j] = materials[i][j].k_aug_n;
				k_aug_p[i][j] = materials[i][j].k_aug_p;
				D_n[i][j] = materials[i][j].D_n * ((materials[i][j].activated==1)? 1:0.00001);
				D_p[i][j] = materials[i][j].D_p * ((materials[i][j].activated==1)? 1:0.00001);
				v_sat_n[i][j] = materials[i][j].v_sat_n;
				v_sat_p[i][j] = materials[i][j].v_sat_p;
			}
		}
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				conducting_x[i][j] = Math.min(materials[i+1][j].conducting, materials[i][j].conducting);
				emfx[i][j] = 0.5*(materials[i+1][j].emf*Math.cos(materials[i+1][j].emf_direction) + materials[i][j].emf*Math.cos(materials[i][j].emf_direction));
				epsx[i][j] = eps0*0.5*(materials[i+1][j].eps_r + materials[i][j].eps_r);
				absorptivity_x[i][j] = Math.min(materials[i+1][j].absorptivity, materials[i][j].absorptivity);
				ac_x[i][j] = (materials[i][j].type == MaterialType.AC_EMF || materials[i+1][j].type == MaterialType.AC_EMF) ? 1:0;
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				conducting_y[i][j] = Math.min(materials[i][j+1].conducting, materials[i][j].conducting);
				emfy[i][j] = 0.5*(materials[i][j+1].emf*Math.sin(materials[i][j+1].emf_direction) + materials[i][j].emf*Math.sin(materials[i][j].emf_direction));
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
				if (n_i[i][j] > 0 || materials[i][j].type == MaterialType.SWITCH) {
					if (updateRho || (rho_n[i][j] == 0 && rho_p[i][j] == 0)) {
						rho_n[i][j] = calcEquilibriumElectronCharge(rho_back[i][j], n_i[i][j]*n_i[i][j]);
						rho_p[i][j] = calcEquilibriumHoleCharge(rho_back[i][j], n_i[i][j]*n_i[i][j]);
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
				emfx[i][j] = 0.5*(materials[i+1][j].emf*Math.cos(materials[i+1][j].emf_direction)
						+ materials[i][j].emf*Math.cos(materials[i][j].emf_direction));
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				emfy[i][j] = 0.5*(materials[i][j+1].emf*Math.sin(materials[i][j+1].emf_direction)
						+ materials[i][j].emf*Math.sin(materials[i][j].emf_direction));
			}
		}
	}

	public double calcEquilibriumElectronCharge(double rho_back, double K) {
		double B = rho_back/e_charge;
		return -e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}

	public double calcEquilibriumHoleCharge(double rho_back, double K) {
		double B = -rho_back/e_charge;
		return e_charge*0.5*(B+Math.sqrt(B*B+4*K));
	}

	public void markJunctions(int smoothing_radius) {

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				needs_smoothing[i][j] = false;
				Material mat = materials[i][j];
				if (conducting[i][j] == 1) {
					
					markNeighborhood(i, j, smoothing_radius);

					for (int di = -smoothing_radius; di <= smoothing_radius; di++) {
						for (int dj = -smoothing_radius; dj <= smoothing_radius; dj++) {
							if (i+di >= 0 && j+dj >= 0 && i+di < nx && j+dj < ny) {
								if (conducting[i+di][j+dj] == 1) {
									if (!visited[i+di][j+dj]) {
										needs_smoothing[i][j] = true;
									}

									if (mat.type != materials[i+di][j+dj].type || mat.cust_id != materials[i+di][j+dj].cust_id) {
										needs_smoothing[i][j] = true;
									}
								}

								distance[i+di][j+dj] = Integer.MAX_VALUE;
								visited[i+di][j+dj] = false;
							}
						}
					}
				}
			}
		}
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
		materials[i][j].initialize();
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
		if (i < 0 || j < 0 || i >= nx || j >= ny)
			return;
		
		if (materials[i][j].cust_id == -1) {
			initializeMaterial(materials[i][j], materials[i][j].type);
		} else {
			if (materialmanager.mat_map.containsKey(materials[i][j].cust_id))
				materials[i][j].copyFrom(materialmanager.mat_map.get(materials[i][j].cust_id));
			else
				materials[i][j].initialize();
		}
	}

	public void initializeMaterial(int i, int j, GeneralMaterialType material) {
		if (i < 0 || j < 0 || i >= nx || j >= ny)
			return;

		if (material.cust_id == -1) {
			initializeMaterial(materials[i][j], material.type);
		} else {
			if (materialmanager.mat_map.containsKey(material.cust_id))
				materials[i][j].copyFrom(materialmanager.mat_map.get(material.cust_id));
			else
				materials[i][j].initialize();
		}
	}


	public void initializeMaterial(Material mat, GeneralMaterialType material) {
		if (material.cust_id == -1) {
			initializeMaterial(mat, material.type);
		} else {
			if (materialmanager.mat_map.containsKey(material.cust_id))
				mat.copyFrom(materialmanager.mat_map.get(material.cust_id));
			else
				mat.initialize();
		}
	}

	public void initializeMaterial(Material mat, MaterialType material) {
		if (mat.modified)
			return;

		mat.type = material;

		if (material == MaterialType.DIELECTRIC) mat.eps_r = dielectric_eps_r;
		else if (material == MaterialType.FERROMAGNET) mat.mu_r = ferromagnet_mu_r;
		else if (material == MaterialType.POS_CHARGE) mat.rho_back = staticcharge_density;
		else if (material == MaterialType.NEG_CHARGE) mat.rho_back = -staticcharge_density;

		if (material.isConducting())
		{
			mat.conducting = 1;
			double Phi = W_metal_default;
			double g = g_metal;
			mat.k_rad = k_rad_metal;
			mat.k_SRH_n = k_SRH_n_semi;
			mat.k_SRH_p = k_SRH_p_semi;
			mat.k_aug_n = k_aug_n_semi;
			mat.k_aug_p = k_aug_p_semi;
			mat.D_n = D_electron_semi;
			mat.D_p = D_hole_semi;
			mat.v_sat_n = v_sat_n_semi;
			mat.v_sat_p = v_sat_p_semi;
			mat.eps_r = eps_r_semi;
			//mat.mu_r = 1/mat.eps_r;

			if (material == MaterialType.METAL_HIGH_W) Phi = W_metal_high;
			else if (material == MaterialType.METAL_LOW_W) Phi = W_metal_low;
			else if (material == MaterialType.METAL_HIGH_C) g = g_metal_high;
			else if (material == MaterialType.METAL_LOW_C) g = g_metal_low;
			else if (material == MaterialType.CURRENT) {
				g = g_currentsource;
				mat.D_n *= currentsource_mobility;
				mat.D_p *= currentsource_mobility;
			}
			
			mat.computeBandstructurePhiG(Phi, Eg_metal, g, g, k*T);
		}

		if (material.isSemiconducting())
		{
			mat.conducting = 1;
			mat.semiconducting = 1;
			mat.k_rad = k_rad_semi;
			mat.k_SRH_n = k_SRH_n_semi;
			mat.k_SRH_p = k_SRH_p_semi;
			mat.k_aug_n = k_aug_n_semi;
			mat.k_aug_p = k_aug_p_semi;
			mat.D_n = D_electron_semi;
			mat.D_p = D_hole_semi;
			mat.v_sat_n = v_sat_n_semi;
			mat.v_sat_p = v_sat_p_semi;
			mat.eps_r = eps_r_semi;
			//mat.mu_r = 1/mat.eps_r;

			if (material == MaterialType.SEMI_P_TYPE) mat.rho_back = -p_default_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_N_TYPE) mat.rho_back = n_default_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_HEAVY_P_TYPE) mat.rho_back = -p_heavy_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_HEAVY_N_TYPE) mat.rho_back = n_heavy_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_LIGHT_P_TYPE) mat.rho_back = -p_light_doping_concentration*e_charge;
			else if (material == MaterialType.SEMI_LIGHT_N_TYPE) mat.rho_back = n_light_doping_concentration*e_charge;

			mat.computeBandstructureChiG(chi_semi, Eg_semi, gc_semi, gv_semi, k*T);
			
			mat.D_n *= calculateMobilityFactor(d_crit_n, mat.rho_back);
			mat.D_p *= calculateMobilityFactor(d_crit_p, mat.rho_back);
		}

		mat.auto_placed = false;
	}
	
	public double calculateMobilityFactor(double d_crit, double donorDensity) {
		double density = Math.abs(donorDensity/e_charge);
		return 1/(Math.sqrt(density/d_crit)+1);
		
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
				
				if (controls.debugging)
					System.out.println("Starting poisson residual: " + Math.sqrt(num/((denom == 0)? 1 : denom)));
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

				if (controls.debugging)
					System.out.println("Poisson residual: " + Math.sqrt(num/((denom == 0)? 1 : denom)));
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
		data += ("t = " + units.toString(time, Quantity.TIME));

		for (Probe p: probes) {
			if (!(p instanceof Ground || p instanceof Ruler)) {
				data += ", ";
				data += p.getText(units);
			}
		}

		data += "\n";

		datastream.print(data);
		datastream.flush();

		renderer.probetexttimer = 30;
	}
	
	public void setDebugInfo() {
		if (!controls.debugging) {
			return;
		}
		
		int mx = controls.mx;
		int my = controls.my;
		Material mat = materials[mx][my];
		
		String str = "";
		str += "Mouse\n";
		str += ("x\t"  						+	units.toString(mx*ds, Quantity.LENGTH) + "\n");
		str += ("y\t"  						+	units.toString(ds*ny-(my+1)*ds, Quantity.LENGTH) + "\n");
		str += "\nFields\n";
		str += ("E\t"  						+	units.toString(Utils.bilinearinterp_length(Ex, Ey, mx, my, nx, ny), Quantity.ELECTRIC_FIELD) + "\n");
		str += ("D\t"  						+	units.toString(Utils.bilinearinterp_length(Dx, Dy, mx, my, nx, ny), Quantity.ELECTRIC_FLUX_DENSITY) + "\n");
		str += ("B\t"  						+	units.toString(parity*Utils.bilinearinterp(Bz, mx-0.5, my-0.5, nx, ny), Quantity.MAGNETIC_FLUX_DENSITY) + "\n");
		str += ("H\t"  						+	units.toString(parity*Utils.bilinearinterp(Hz, mx-0.5, my-0.5, nx, ny), Quantity.MAGNETIC_FIELD_STRENGTH) + "\n");
		str += "\nCharge/current\n";
		str += ("\u03c1\u2099\t"  			+	units.toString(rho_n[mx][my], Quantity.CHARGE_DENSITY) + "\n");
		str += ("\u03c1\u209A\t"  			+	units.toString(rho_p[mx][my], Quantity.CHARGE_DENSITY) + "\n");
		str += ("\u03c1\u2080\t"			+	units.toString(rho_back[mx][my], Quantity.CHARGE_DENSITY) + "\n");
		str += ("\u03c1 abs\t"				+	units.toString(rho_abs[mx][my], Quantity.CHARGE_DENSITY) + "\n");
		str += ("\u03c1\t"  				+	units.toString(rho_free[mx][my], Quantity.CHARGE_DENSITY) + "\n");
		str += ("J\u2099\t" 				+	units.toString(Utils.bilinearinterp_length(Jx_n, Jy_n, mx, my, nx, ny), Quantity.CURRENT_DENSITY) + "\n");
		str += ("J\u209A\t"  				+	units.toString(Utils.bilinearinterp_length(Jx_p, Jy_p, mx, my, nx, ny), Quantity.CURRENT_DENSITY) + "\n");
		str += ("J abs\t"  					+	units.toString(Utils.bilinearinterp_length(Jx_abs, Jy_abs, mx, my, nx, ny), Quantity.CURRENT_DENSITY) + "\n");
		str += ("J\t"  						+	units.toString(Utils.bilinearinterp_length(Jx_free, Jy_free, mx, my, nx, ny), Quantity.CURRENT_DENSITY) + "\n");
		str += "\nThermodynamic quantities\n";
		str += ("\u03d5\t"  				+	units.toString(phi[mx][my], Quantity.ELECTRIC_POTENTIAL) + "\n");
		str += ("F\u2099\t"  				+	units.toString(mu_n[mx][my]/eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "\n");
		str += ("F\u209a\t"  				+	units.toString(-mu_p[mx][my]/eVtoJ, Quantity.ELECTRIC_POTENTIAL) + "\n");
		str += ("V\t"  						+	units.toString(V_avg[mx][my], Quantity.ELECTRIC_POTENTIAL) + "\n");
		str += ("G\t"  						+	units.toString(G[mx][my], Quantity.RATE_DENSITY) + "\n");
		str += ("R\t"  						+	units.toString(R[mx][my], Quantity.RATE_DENSITY) + "\n");
		str += ("Electron CMF\t"  			+	units.toString(Utils.bilinearinterp_length(cmfy_n, cmfy_n, mx, my, nx, ny), Quantity.FORCE) + "\n");
		str += ("Hole CMF\t"  				+	units.toString(Utils.bilinearinterp_length(cmfy_p, cmfy_p, mx, my, nx, ny), Quantity.FORCE) + "\n");
		str += ("EMF\t"  					+	units.toString(Utils.bilinearinterp_length(emfx, emfy, mx, my, nx, ny), Quantity.ELECTRIC_FIELD) + "\n");
		str += ("Electron vel.\t"  			+	units.toString(Utils.bilinearinterp_length(Jx_n, Jy_n, mx, my, nx, ny)/rho_n[mx][my], Quantity.VELOCITY) + "\n");
		str += ("Hole vel.\t"  				+	units.toString(Utils.bilinearinterp_length(Jx_p, Jy_p, mx, my, nx, ny)/rho_p[mx][my], Quantity.VELOCITY) + "\n");
		str += "\nMaterial properties\n";
		str += ("Type\t"  					+	mat.type.name + "\n");
		str += ("\u2130\t"  				+	units.toString(mat.emf, Quantity.ELECTRIC_FIELD) + "\n");
		str += ("\u03b5/\u03b5\u2080\t" 	+	units.toString(mat.eps_r, Quantity.DIMENSIONLESS) + "\n");
		str += ("\u03bc/\u03bc\u2080\t"  	+	units.toString(mat.mu_r, Quantity.DIMENSIONLESS) + "\n");
		str += ("Absorptivity\t"  			+	units.toString(mat.absorptivity, Quantity.DIMENSIONLESS) + "\n");
		str += ("Conducting\t"  			+	mat.conducting + "\n");
		str += ("Semicond.\t"  				+	mat.semiconducting + "\n");
		str += "\nNumerical stability ratios\n";
		str += CFL_text;
		
		opts.textPane.setText(str);
	}

	public boolean hasGround() {
		for (Probe p : probes) {
			if (p instanceof Ground)
				return true;
		}
		
		return false;
	}
	
	public Ground getGround() {
		for (Probe p : probes) {
			if (p instanceof Ground)
				return (Ground)p;
		}
		
		return null;
	}
	
	public boolean hasRuler() {
		for (Probe p : probes) {
			if (p instanceof Ruler)
				return true;
		}
		
		return false;
	}
	
	public Ruler getRuler() {
		for (Probe p : probes) {
			if (p instanceof Ruler)
				return (Ruler)p;
		}
		
		return null;
	}
	
	public void relabelProbes() {
		probe_index = 0;
		for (Probe p : probes) {
			p.name = getProbeName();
			probe_index++;
		}
	}
	
	public void addProbe(Probe p) {
		probes.add(p);
		relabelProbes();
	}
	
	public void removeProbe(Probe p) {
		probes.remove(p);
		relabelProbes();
	}
	
	int probe_index = 0;
	public String getProbeName() {
		if (probe_index < 26) return String.valueOf((char)('a'+probe_index));
		return String.valueOf(probe_index - 26);
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
}