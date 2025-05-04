package electrodynamics;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.TimerTask;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.CyclicBarrier;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.stream.JsonReader;

public class Electrodynamics extends TimerTask implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener, ActionListener {
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
	
	//TODO:
	// Option: change vf density and toggle vf/lines
	// Change size & sim parameters
	// Interactions with light
	// Band structure diagrams
	// Undo/redo
	
	/* Dynamical simulation variables */
	
	double[][] Ex;			// x component of E field
	double[][] Ey;			// y component of E field
	double[][] Hz;			// z component of B field
	double[][] Hz_laplacian;
	
	double[][] rho_n;		// Electrons
	double[][] rho_p;		// Holes
	double[][] rho_back;	// Background
	double[][] rho_abs;		// Absorber charge
	double[][] rho_free;	// Total free charges
	double[][] mobility_factor;
	
	double[][] Jx_n;		// Electron current
	double[][] Jy_n;
	double[][] Jx_p;		// Hole current
	double[][] Jy_p;
	double[][] Jx_abs;		// Absorber current
	double[][] Jy_abs;
	double[][] Jx_free;		// Total free current
	double[][] Jy_free;

	/* Miscellaneous fields used for display */
	
	double[][] Dx;				// Electric displacement field
	double[][] Dy;
	double[][] Bz;				// B field
	double[][] Sx;				// Poynting vector
	double[][] Sy;
	double[][] u;				// EM energy density
	double[][] phi;				// Electric scalar potential

	double[][] G;				// Generation rate
	double[][] F_n;				// Total chemical potential of electrons (quasi-Fermi level minus the electrostatic potential)
	double[][] F_p;				// Total chemical potential of holes
	double[][] grad_E0x_n;		// Gradient of chemical energy
	double[][] grad_E0y_n;
	double[][] grad_E0x_p;
	double[][] grad_E0y_p;
	double[][] grad_Fx_n;		// Gradient of total chemical potential
	double[][] grad_Fy_n;
	double[][] grad_Fx_p;
	double[][] grad_Fy_p;
	double[][] Q;				// Heat dissipation rate
	double[][] S;				// Free energy dissipation rate
	double[][] F;				// Average total electrochemical potential

	boolean updateMiscFields = false;
	double[][] debug;
	
	/* Material parameters */
	
	Material[][] materials;
	Material[][] selection;
	Material[][] clipboard;
	
	double[][] F0_n;		// Standard chemical potential of electrons
	double[][] F0_p;		// Standard chemical potential of holes
	double[][] E0_n;		// Standard chemical energy of electrons
	double[][] E0_p;		// Standard chemical energy of holes
	double[][] K;			// Charge carrier equilibrium constant
	double[][] E_a;			// Activation energy of recombination
	double[][] R;			// Recombination rate constant
	
	double[][] cmfx_n;		// Chemical-motive force for electrons
	double[][] cmfy_n;

	double[][] cmfx_p;		// Chemical-motive force for holes
	double[][] cmfy_p;
	
	double[][] emfx;		// External electromotive force
	double[][] emfy;
	double[][] epsx;		// Dielectric constant
	double[][] epsy;
	double[][] mu_z;		// Relative permeability
	
	int[][] conducting;		// Does the material have partially filled bands
	int[][] conducting_x;
	int[][] conducting_y;
	
	double[][] absorptivity;
	double[][] absorptivity_x;
	double[][] absorptivity_y;

	/* Multigrid Poisson eq solver */
	
	int log2_resolution;
	double[][] MG_rho0;
	double[][][] MG_rho;
	double[][][] MG_epsx;
	double[][][] MG_epsy;
	double[][] MG_eps_avg;
	double[][] MG_phi1;
	double[][] MG_phi2;

	/* Pathfinding */
	
	int[][] distance;
	boolean[][] visited;
	
	
	/* Physical constants */
	
	double eps0 = 8.85e-12;
	double mu0 = 1.257e-6;
	double c = 1/Math.sqrt(eps0*mu0);
	double c_squared = 1/(eps0*mu0);
	double kT = 4.11e-21;
	double beta = 1/kT;
	double e_charge = 1.6e-19;
	double m_electron = 9e-31;

	
	/* Material properties */
	
	double mu_electron = 0.1400*1500;					// Electron mobility
	double D_electron = mu_electron/(beta*e_charge);	// Diffusion constant, determined by Einstein relation
	
	double mu_hole = 0.5*mu_electron;					// Hole mobility, slightly lower than electron
	double D_hole = mu_hole/(beta*e_charge);

	double eVtoJ = e_charge;				// eV to J conversion factor

	double ni_semi = 1e16;					// Semiconductor equilibrium concentration
	double W_semi = 4.7*eVtoJ;				// Semiconductor work function
	double E_b_semi = 1.12*eVtoJ;			// Semiconductor band gap

	double ni_metal = 5e20;				// Metal charge carrier concentration
	double W_metal_default = W_semi;		// Metal work function
	double W_metal_high = W_semi + 0.3*eVtoJ;
	double W_metal_low = W_semi - 0.3*eVtoJ;
	double E_b_metal = E_b_semi;

	double recombination_rate_semi = 1e7/ni_semi;
	double recombination_rate_metal = 1e8/ni_semi;
	double K_semi = ni_semi*ni_semi;
	
	double recombination_cross_section = 1e-10;
	double arrhenius_prefactor = recombination_cross_section*Math.sqrt(8*kT/(Math.PI*m_electron/2));
	double E_a_semi = -Math.log(recombination_rate_semi/arrhenius_prefactor);
	double E_a_metal = -Math.log(recombination_rate_metal/arrhenius_prefactor);
	
	double q_n = -e_charge;				// Charge of single electron
	double q_p = e_charge;				// Charge of single hole
	
	double n_default_doping_concentration = 5e19;
	double p_default_doping_concentration = 5e19;
	double n_light_doping_concentration = 1e19;
	double p_light_doping_concentration = 1e19;
	double n_heavy_doping_concentration = 2.5e20;
	double p_heavy_doping_concentration = 2.5e20;
	
	double dielectric_eps_r = 5.0;
	double ferromagnet_mu_r = 5.0;
	double staticcharge_density = 10.0;
	
	double E_sat = 5e5;					// Maximum electric field before velocity saturates
	
	int junction_size = 3;
	
	
	/* Probes */
	
	List<VoltageProbe> voltageprobes;
	List<CurrentProbe> currentprobes;
	VoltageProbe ground = null;
	
	
	/* Domain parameters */
	
	double default_width = 2.56e-5;
	double width;				// Width, in SI
	double ds;					// Spatial discretization
	double dt;					// Timestep
	double depth = 1e-3;		// Extent of circuit in z-dimension, only used to get reasonable values for current probe
	int default_resolution = 256;
	int nx;						// Number of grid points in x-dimension
	int ny;						// Number of grid points in y-dimension
	int absorber_width;
	double absorbing_coeff;
	double Hz_dissipation;
	double dt_maximum;
	int parity = -1;			// Negative sign resulting from flipped y-axis in graphics coordinate system
	
	double time = 0.0;
	long stepnumber = 0;
	long frame = 0;
	int lastsimspeed = 0;
	
	double error_detection_threshold = 1e-5;
	boolean sign_violation = false;

	
	/* Multithreading */
	
	int n_threads = Runtime.getRuntime().availableProcessors();
	//int n_threads = 5;
	CyclicBarrier start_barrier = new CyclicBarrier(n_threads + 1);
	CyclicBarrier stop_barrier = new CyclicBarrier(n_threads + 1);
	CyclicBarrier mid_barrier = new CyclicBarrier(n_threads);
	ArrayList<SimulationThread> sim_threads = new ArrayList<SimulationThread>();
	
	
	/* Graphics */
	
	RenderCanvas r;
	MainWindow opts;
	HelpDialog help;
	BufferedImage screen;
	float[][] image_r;
	float[][] image_g;
	float[][] image_b;
	float col_r = 0;
	float col_g = 0;
	float col_b = 0;
	double alphaBG = 0;
	double alphaFG = 0;
	Random rand = new Random();

	ArrayList<Text> texts = new ArrayList<Text>();
	Font bigfont = new Font(Font.SANS_SERIF, Font.PLAIN, 15);
	Font regularfont = new Font(Font.SANS_SERIF, Font.PLAIN, 12);
	int scalefactor;
	int imgwidth = 0;
	int imgheight = 0;
	int targetframerate = 60;
	int frameduration = 1000/targetframerate;
	
	
	/* Performance profiling */
	
	Timer t4 = new Timer("Poisson constraint solver", true);
	Timer t7 = new Timer("Poisson potential solver", true);
	Timer t5 = new Timer("Graphics", true);
	Timer t6 = new Timer("Iterate simulation", true);
	Timer t8 = new Timer("Calc misc fields", true);
	Timer t9 = new Timer("Debug", false);

	
	/* Keyboard controls */
	
	boolean advanceframe = false;
	boolean clear = false;
	boolean reset = false;
	boolean save = false;
	boolean load = false;
	boolean debugging = false;
	boolean cut = false;
	boolean copy = false;
	boolean paste = false;
	boolean delete = false;
    boolean shift_down = false;
    boolean ctrl_down = false;
    boolean alt_down = false;
    
	
	/* Mouse controls */

	PointerInfo pointerinfo = MouseInfo.getPointerInfo();
	boolean mouse_pressed = false;
	boolean mouse_pressed_prev = false;
	boolean modifier_pressed = false;
	boolean moving_selection = false;
	boolean dragging_selection = false;
	boolean brush_changed = false;

	int mousebutton = 0;
	int mx = 0;
	int my = 0;
	int mx_start = 0;
	int my_start = 0;
	
	int mx_index = 0;
	int my_index = 0;
	int mx_start_index = 0;
	int my_start_index = 0;

	double mx_realspace = 0;
	double my_realspace = 0;
	double mxp_realspace = 0;
	double myp_realspace = 0;
	double mx_start_realspace = 0;
	double my_start_realspace = 0;

	int delta_mx_index = 0;
	int delta_my_index = 0;
	
	boolean EMF_selected = false;
	double max_EMF = 5e5;
	
	Brush prev_brush;
	double brushsize = 0;
	int prev_EMF_setting = 0;
	BoundaryCondition prev_boundary = null;
	
	boolean[][] under_brush;
	boolean[][] selected;
	boolean[][] selected_EMF;
	
	int text_x = 0;
	int text_y = 0;
	boolean texting = false;

	Cursor HAND_CURSOR = new Cursor(Cursor.HAND_CURSOR);
	Cursor DEFAULT_CURSOR = new Cursor(Cursor.DEFAULT_CURSOR);
	
	
	/* Saving and loading */
	
	public static File infile;
	public static File outfile;
	int saveversion = 1;
	String fileextension = ".semisim";
	String startingpath = ".";
	

	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(
					UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}

		Electrodynamics w = new Electrodynamics();
		java.util.Timer t = new java.util.Timer();
		for (int i = 0; i < w.n_threads; i++) {
			w.sim_threads.add(w.new SimulationThread(i, w.n_threads, w.nx));
		}

		for (int i = 0; i < w.n_threads; i++) {
			w.sim_threads.get(i).start();
		}
		
		t.schedule(w, 0, w.frameduration);
	}

	public Electrodynamics() {
		opts = new MainWindow();
		r = new RenderCanvas(this);
		r.setFocusable(true);
		
		detect64Bit();

		initializeGrid(default_width, default_resolution);
		opts.add(r, BorderLayout.CENTER);
		r.addMouseListener(this);
		r.addMouseMotionListener(this);
		r.addMouseWheelListener(this);
		r.addKeyListener(this);
		opts.gui_reset.addActionListener(this);
		opts.gui_resetall.addActionListener(this);
		opts.gui_save.addActionListener(this);
		opts.gui_open.addActionListener(this);
		opts.gui_help.addActionListener(this);
		opts.gui_editdesc.addActionListener(this);
		opts.gui_view.addActionListener(this);
		opts.gui_view_vec.addActionListener(this);
		opts.gui_brush.addActionListener(this);
		opts.gui_material.addActionListener(this);
		
		opts.gui_material.removeItem(MaterialType.ABSORBER);
		//opts.gui_material.removeItem(MaterialType.SWITCH);
		
		opts.gui_view.removeItem(ScalarView.DEBUG);
		
		opts.pack();
		
		addKeyBinds(r);
		
		InputMap im = (InputMap)UIManager.get("Button.focusInputMap");
		im.put(KeyStroke.getKeyStroke("pressed SPACE"), "none");
		im.put(KeyStroke.getKeyStroke("released SPACE"), "none");

		
		help = new HelpDialog();
		help.setVisible(false);
	}
	
	public void detect64Bit() {
		if (!System.getProperty("sun.arch.data.model").equals("64"))
		{
			int result = JOptionPane.showConfirmDialog(opts, "Running this application on a 32-bit platform may cause some issues. Do you still wish to proceed?", "Warning", JOptionPane.YES_NO_OPTION);
			if (result != JOptionPane.OK_OPTION)
			{
				System.exit(0);
			}
		}
	}
	
	public void run() {
		try {
			int iterationmultiplier = opts.gui_simspeed_2.getValue();

			if (clear) {
				resetFields(false);
				multigridSolve(true, false);
				time = 0.0;
				clear = false;
			}

			if (reset) {
				int result = JOptionPane.showConfirmDialog(opts, "Do you wish to reset the entire simulation?", "Reset", JOptionPane.YES_NO_OPTION);
				if (result == JOptionPane.OK_OPTION)
				{
					resetFields(true);
					time = 0.0;
				}
				reset = false;
			}

			if (opts.gui_simspeed.getValue() != lastsimspeed) {
				lastsimspeed = opts.gui_simspeed.getValue();
				dt = dt_maximum*(lastsimspeed/20.0);
			}

			handleMouseInput();

			if (!opts.gui_paused.isSelected() || advanceframe) {
				for (int i = 0; i < iterationmultiplier ; i++) {
					start_barrier.await();
					stop_barrier.await();
				}

				if (frame % 2 == 0)
					calcMiscFields(true);
				else
					calcMiscFields(false);

				frame++;
			} else if (updateMiscFields) {
				calcMiscFields(true);
			}


			if (save) {
				writeFile();
				save = false;
			}

			if (load) {
				readFile();
				load = false;
			}

			r.repaint();
		} catch (Exception e) {
			JOptionPane.showConfirmDialog(opts, e.getMessage(), "Error", JOptionPane.OK_OPTION);
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public void initializeGrid(double width, int resolution) {
		ds = width/resolution;
		nx = resolution;
		ny = resolution;

		this.width = nx*ds;
		
		log2_resolution = (int) Math.round(Math.log(resolution)/Math.log(2));
		assert(1 << log2_resolution == resolution);

		dt_maximum = 0.9*ds/(Math.sqrt(2)*c);

		Hz_dissipation = 0.01*ds*ds/dt_maximum;

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
		selection = new Material[nx][ny];
		clipboard = new Material[nx][ny];
		K = new double[nx][ny];
		F0_n = new double[nx][ny];
		F0_p = new double[nx][ny];
		F_n = new double[nx][ny];
		F_p = new double[nx][ny];
		E0_n = new double[nx][ny];
		E0_p = new double[nx][ny];
		E_a = new double[nx][ny];
		R = new double[nx][ny];
		cmfx_n = new double[nx][ny];
		cmfy_n = new double[nx][ny];
		cmfx_p = new double[nx][ny];
		cmfy_p = new double[nx][ny];
		conducting = new int[nx][ny];
		conducting_x = new int[nx][ny];
		conducting_y = new int[nx][ny];
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
		
		under_brush = new boolean[nx][ny];
		selected = new boolean[nx][ny];
		selected_EMF = new boolean[nx][ny];

		voltageprobes = new CopyOnWriteArrayList<VoltageProbe>();		
		currentprobes = new CopyOnWriteArrayList<CurrentProbe>();

		resetFields(true);
		multigridSolve(true, false);
		calcMiscFields(true);

		text_x = 0;
		text_y = 0;
		texting = false;
		
		image_r = new float[nx][ny];
		image_g = new float[nx][ny];
		image_b = new float[nx][ny];
		scalefactor = (int)(768.0/ny);
		
		imgwidth = (int)Math.ceil(scalefactor*nx);
		imgheight = (int)Math.ceil(scalefactor*ny);
		r.setPreferredSize(new Dimension(imgwidth, imgheight));

		opts.setVisible(true);
		
		screen = (BufferedImage) opts.createImage(imgwidth, imgheight);
	}
	
	public void resetFields(boolean resetall) {
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{

				if (resetall) {
					if (materials[i][j] == null)
						materials[i][j] = new Material();
					else
						materials[i][j].erase();

					if (selection[i][j] == null)
						selection[i][j] = new Material();
					else
						selection[i][j].erase();

					if (clipboard[i][j] == null)
						clipboard[i][j] = new Material();
					
					F0_n[i][j] = 0;
					F0_p[i][j] = 0;
					F_n[i][j] = 0;
					F_p[i][j] = 0;
					E0_n[i][j] = 0;
					E0_p[i][j] = 0;
					K[i][j] = 0;
					E_a[i][j] = 0;
					R[i][j] = 0;
					
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

				MG_phi1 [i][j] = 0.0;
				MG_phi2 [i][j] = 0.0;
				
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

				selected[i][j] = false;
				selected_EMF[i][j] = false;
			}
		}

		if (resetall) {
			voltageprobes.clear();
			currentprobes.clear();
			ground = null;

			opts.setTitle("Brandon's semiconductor simulator");
		}
		

		constructBoundary();

		initializeAllMaterials();
		updateAllMaterials();
		checkCFL();
	}
	
	public void constructBoundary() {
		absorber_width = (int)(0.05*nx);
		absorbing_coeff = 50*c/this.width;
		double max_stretch = 10;
		double coeff = Math.log(max_stretch);
		
		if ((BoundaryCondition)opts.gui_bc.getSelectedItem() == BoundaryCondition.DISSIPATIVE) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					double dx = (double)Math.max(0, absorber_width-Math.min(i, nx-1-i))/absorber_width;
					double dy = (double)Math.max(0, absorber_width-Math.min(j, ny-1-j))/absorber_width;
					double depth = Math.sqrt(dx*dx+dy*dy);
					if (depth > 0) {
						double stretchfactor = Math.exp(coeff*depth);
						
						materials[i][j].erase();
						materials[i][j].type = MaterialType.ABSORBER;
						materials[i][j].eps_r = stretchfactor;
						materials[i][j].mu_r = stretchfactor;
						materials[i][j].absorptivity = 1;
					}
				}
			}
		} else {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					if (materials[i][j].type == MaterialType.ABSORBER) {
						materials[i][j].erase();
					}
				}
			}
		}
	}
	
	class SimulationThread extends Thread {
		
		int i_min;
		int i_max;
		int n_thread;
		
		public SimulationThread(int n, int n_threads, int nx) {
			i_min = (n*nx)/n_threads;
			i_max = (n+1)*nx/n_threads-1;
			n_thread = n;
			System.out.println("Thread " + n_thread + ": " + i_min + " < i <= " + i_max);
		}
		
		@Override
		public void run() {
			try {
				while (true) {
					start_barrier.await();
					
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

					for (int i = 1; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-1; j++)
							{	
								double generation_rate = conducting[i][j]*R[i][j]*(K[i][j] - rho_n[i][j]*rho_p[i][j]/(q_n*q_p));

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

					for (int i = 0; i < nx-1; i++)
					{
						if (i >= i_min && i <= i_max) {
							for (int j = 1; j < ny-1; j++)
							{
								double ex_prev = Ex[i][j];

								double mf = Math.min(mobility_factor[i+1][j], mobility_factor[i][j]);
								//double mobility_factor = Math.min(1, E_sat/Math.abs(Ex[i][j]));

								double sigma_n = conducting_x[i][j]*mf*mu_electron*logmean(-rho_n[i+1][j],-rho_n[i][j]);
								double sigma_p = conducting_x[i][j]*mf*mu_hole*logmean(rho_p[i+1][j], rho_p[i][j]);

								Jx_abs[i][j] = 0;

								Jx_n[i][j] = conducting_x[i][j]*(-mf*D_electron*(rho_n[i+1][j] - rho_n[i][j])/ds
										+ sigma_n*(emfx[i][j] + cmfx_n[i][j]/q_n));

								Jx_p[i][j] = conducting_x[i][j]*(-mf*D_hole*(rho_p[i+1][j] - rho_p[i][j])/ds
										+ sigma_p*(emfx[i][j] + cmfx_p[i][j]/q_p));

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

								double sigma_n = conducting_y[i][j]*mf*mu_electron*logmean(-rho_n[i][j+1],-rho_n[i][j]);
								double sigma_p = conducting_y[i][j]*mf*mu_hole*logmean(rho_p[i][j+1], rho_p[i][j]);

								Jy_abs[i][j] = 0;

								Jy_n[i][j] = conducting_y[i][j]*(-mf*D_electron*(rho_n[i][j+1] - rho_n[i][j])/ds
										+ sigma_n*(emfy[i][j] + cmfy_n[i][j]/q_n));

								Jy_p[i][j] = conducting_y[i][j]*(-mf*D_hole*(rho_p[i][j+1] - rho_p[i][j])/ds
										+ sigma_p*(emfy[i][j] + cmfy_p[i][j]/q_p));

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

						if (stepnumber%500 == 0) {
							multigridSolve(true, false);
						}

						time += dt;

						advanceframe = false;
					}
					
					stop_barrier.await();
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}
	
	public double logmean(double x, double y)
	{
		if (x <= 0 || y <= 0)
			return 0;
		
		//My approximation
		if (Math.abs((x-y)/(x+y)) <  1e-3)
			return (2/3.0)*Math.sqrt(x*y) + (1/6.0)*(x+y);

		return (x-y)/lut.log(x/y);
		//return (x-y)/Math.log(x/y);
	}

	LogLUT lut = new LogLUT();
	
	// Quick and dirty way to calculate log precisely (relative error < 10^-9)
	class LogLUT {
		//int min_exponent = -1022;
		//int orders = 2045;
		int min_exponent = -100;
		int orders = 200;
		//int divisions = 4;
		int length;
		long A = (long)0b1111111111 << 52;
		
		double[] x;
		double[] log_x;
		
		public LogLUT() {
			length = orders*4;
			x = new double[length];
			log_x = new double[length];
			
			for (int i = 0; i < length; i++) {
				x[i] = Math.pow(2, min_exponent+i/4)*(1+(i%4)/4.0);
				log_x[i] = Math.log(x[i]);
			}
		}
		
		public double log(double y) {
			int exp = Math.getExponent(y);
			double mantissa = Double.longBitsToDouble(A|Double.doubleToRawLongBits(y)&(~0x7ff0000000000000l));
			int i = (int)(((exp-min_exponent)<<2) + 4*mantissa - 3.5);
			double w = y/x[i];
			return (w-1)/(0.66666666666666667*Math.sqrt(w) + 0.16666666666666667*(w+1)) + log_x[i];
		}
	}
	
	public void updateAllMaterials() {
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
			}
		}
		
		computeChemicalForces();
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (K[i][j] > 0 || materials[i][j].type == MaterialType.SWITCH) {
					if (rho_n[i][j] == 0 && rho_p[i][j] == 0) {
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
			}
		}
		
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				Jx_n[i][j] *= conducting_x[i][j];
				Jx_p[i][j] *= conducting_x[i][j];
				Jx_abs[i][j] *= (materials[i][j].absorptivity > 0)? 1:0;
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
						if (distance[i][j] < smallest_length && visited[i][j] == false && conducting[i][j] == 1) {
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
					R[i][j] = arrhenius_prefactor*Math.exp(-E_a[i][j]);
				} else {
					E_a[i][j] = 0;
					K[i][j] = 0;
					R[i][j] = 0;
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
	
	public double[][] copyArray(double[][] array) {
		double [][] newarray = new double[array.length][];
		for(int i = 0; i < array.length; i++)
			newarray[i] = array[i].clone();
		return newarray;
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
			else if (material == MaterialType.METAL_HIGH_C) materials[i][j].ni = 2.5*ni_metal;
			else if (material == MaterialType.METAL_LOW_C) materials[i][j].ni = 0.25*ni_metal;
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
	}
	
	public void handleMouseInput() {
		Brush brush = (Brush) opts.gui_brush.getSelectedItem();
		BrushShape brushshape = (BrushShape) opts.gui_brush_1.getSelectedItem();
		
		boolean pressing = false;
		boolean releasing = false;
		
		if (mouse_pressed) {
			if (!mouse_pressed_prev) {
				pressing = true;
				r.requestFocus();
				if (Brush.isMaterialModifyingBrush(brush) || brush == Brush.SELECT)
					opts.gui_paused.setSelected(true);
			}
		} else {
			if (mouse_pressed_prev) {
				releasing = true;
			}
		}
		mouse_pressed_prev = mouse_pressed;
		

		if (!mouse_pressed && !releasing) {
			
			if (ctrl_down || shift_down) {
				if (!modifier_pressed && Brush.isMaterialModifyingBrush(brush)) {
					modifier_pressed = true;
					prev_brush = brush;
					
					if (ctrl_down) {
						opts.gui_brush.setSelectedItem(Brush.FILL);
						brush = (Brush) opts.gui_brush.getSelectedItem();
					}
					else if (shift_down) {
						opts.gui_brush.setSelectedItem(Brush.LINE);
						brush = (Brush) opts.gui_brush.getSelectedItem();
					}
					//r.requestFocus();
				}
			} else {
				if (modifier_pressed) {
					modifier_pressed = false;
					opts.gui_brush.setSelectedItem(Brush.DRAW);
					brush = (Brush) opts.gui_brush.getSelectedItem();
					//r.requestFocus();
				}
			}
		}
		
		mx_realspace = Math.round((mx-1)/(double)scalefactor - 0.5)*ds;
		my_realspace = Math.round((my-1)/(double)scalefactor - 0.5)*ds;

		mx_start_realspace =  Math.round((mx_start-1)/(double)scalefactor - 0.5)*ds;
		my_start_realspace = Math.round((my_start-1)/(double)scalefactor - 0.5)*ds;

		mx_index = (int)Math.round((mx-1)/(double)scalefactor - 0.5);
		my_index = (int)Math.round((my-1)/(double)scalefactor - 0.5);

		mx_start_index = (int)Math.round((mx_start-1)/(double)scalefactor - 0.5);
		my_start_index = (int)Math.round((my_start-1)/(double)scalefactor - 0.5);
		
		if (mx_index < 0) mx_index = 0;
		if (my_index < 0) my_index = 0;
		if (mx_index >= nx) mx_index = nx-1;
		if (my_index >= ny) my_index = ny-1;

		if (mx_start_index < 0) mx_start_index = 0;
		if (my_start_index < 0) my_start_index = 0;
		if (mx_start_index >= nx) mx_start_index = nx-1;
		if (my_start_index >= ny) my_start_index = ny-1;

		opts.gui_stepsizelbl.setText("Step size: " + getSI(dt, "s"));
		opts.gui_stepslbl.setText("Steps/frame: " + opts.gui_simspeed_2.getValue());
		
		brushsize = (width/500)*(Math.pow(10.0, opts.gui_brushsize.getValue()/500.0) + opts.gui_brushsize.getValue()/100.0);
		opts.lblBrushSize.setText("Brush size: " + (int)Math.ceil(brushsize/ds));
		
		if (!Brush.isMaterialModifyingBrush(brush))
		{
			opts.gui_parameter2.setVisible(false);
			opts.gui_parameter2_text.setVisible(false);
			opts.gui_parameter2_text.setText("");
		}
		
		
		if (Brush.isMaterialModifyingBrush(brush) && brush != Brush.FILL) {
			opts.gui_brush_1.setVisible(true);
			opts.gui_brush_highlight.setVisible(true);
			opts.gui_brushsize.setVisible(true);
			opts.lblBrushSize.setVisible(true);
		} else {
			opts.gui_brush_1.setVisible(false);
			opts.gui_brush_highlight.setVisible(false);
			opts.gui_brushsize.setVisible(false);
			opts.lblBrushSize.setVisible(false);
		}
		

		if (Brush.isMaterialModifyingBrush(brush) && brush != Brush.ERASE) {
			opts.gui_material.setVisible(true);
		} else {
			opts.gui_material.setVisible(false);
		}
		
		if (brush_changed) {
			if (!(brush == Brush.SELECT || brush == Brush.FLOODSELECT)) {
				for (int i = 0; i < nx; i++)
				{
					for (int j = 0; j < ny; j++)
					{
						selected[i][j] = false;
					}
				}
			}
			
			if (brush != Brush.INTERACT) {
				for (int i = 0; i < nx; i++)
				{
					for (int j = 0; j < ny; j++)
					{
						selected_EMF[i][j] = false;
					}
				}
				EMF_selected = false;
			}
			
			brush_changed = false;
		}
		
		boolean update = false;
		
		if (cut || copy) {
			int i_min = nx-1;
			int j_min = ny-1;
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					clipboard[i][j].erase();
					if (selected[i][j] && materials[i][j].type != MaterialType.VACUUM) {
						if (i < i_min) i_min = i;
						if (j < j_min) j_min = j;
					}
				}
			}
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					if (selected[i][j] && materials[i][j].type != MaterialType.VACUUM) {
						clipboard[i-i_min][j-j_min] = materials[i][j].clone();
						if (cut) {
							materials[i][j].erase();
						}
					}
					if (cut) {
						selected[i][j] = false;
					}
				}
			}
			
			if (cut) update = true;
			
			cut = false;
			copy = false;
		}
		
		if (paste) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					selected[i][j] = false;
				}
			}
			
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					selection[i][j] = clipboard[i][j].clone();
				}
			}
			moving_selection = true;
			dragging_selection = false;
			paste = false;
			update = true;
		}
		
		if (delete) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					if (selected[i][j] && materials[i][j].type != MaterialType.VACUUM) {
						materials[i][j].erase();
					}
				}
			}
			delete = false;
			update = true;
		}
		
		switch(brush) {
		case DRAW:
		case LINE:
		case REPLACE:
		case ERASE:
		case FILL:

			if ((mousebutton == MouseEvent.BUTTON2 || alt_down) && pressing) {
				opts.gui_material.setSelectedItem(materials[mx_index][my_index].type);
			}

			double angle = 0;
			MaterialType mat = (MaterialType) opts.gui_material.getSelectedItem();
			
			if (mat == MaterialType.EMF) {
				opts.gui_parameter2.setVisible(true);
				opts.gui_parameter2_text.setVisible(true);

				int directionval = (int)(opts.gui_parameter2.getValue()/6);
				if (directionval == 0) {
					opts.gui_parameter2_text.setText("EMF direction: Up");
					angle = -Math.PI/2;
				}
				if (directionval == 1) {
					opts.gui_parameter2_text.setText("EMF direction: Right");
					angle = 0;
				}
				if (directionval == 2) {
					opts.gui_parameter2_text.setText("EMF direction: Down");
					angle = Math.PI/2;
				}
				if (directionval == 3) {
					opts.gui_parameter2_text.setText("EMF direction: Left");
					angle = Math.PI;
				}
				//opts.gui_parameter2_text.setText("Brush orientation: " + directionval*(360/24) + " deg");
				//angle = Math.PI * directionval/12.0;
			} else {
				opts.gui_parameter2.setVisible(false);
				opts.gui_parameter2_text.setVisible(false);
				opts.gui_parameter2_text.setText("");
			}
			

			if (mousebutton == MouseEvent.BUTTON3 || brush == Brush.ERASE)
				mat = MaterialType.VACUUM;

			if (!(mousebutton == MouseEvent.BUTTON2 || alt_down)) {
				if (brush == Brush.LINE) {
					if (releasing) {
						drawMaterialLine(mx_start_realspace, my_start_realspace, mx_realspace, my_realspace, brush, brushshape, mat, brushsize, angle);
					}
				} else if (brush == Brush.FILL) {
					if (pressing) {
						floodFillSet(mx_index, my_index, materials[mx_index][my_index].type, mat, angle);
					}
				} else if (mouse_pressed) {
					drawMaterialLine(mxp_realspace, myp_realspace, mx_realspace, my_realspace, brush, brushshape, mat, brushsize, angle);
				}
			}

			if (Brush.isBrushShapeImportant(brush)) {
				for (int i = 0; i < nx; i++)
				{
					for (int j = 0; j < ny; j++)
					{
						if (opts.gui_brush_highlight.isSelected()) {
							double cx = 0;
							double cy = 0;
							cx = i*ds;
							cy = j*ds;

							double px = (cx-mx_realspace);
							double py = (cy-my_realspace);
							double r = 0;

							if (brushshape == BrushShape.CIRCLE)
								r = Math.sqrt(px*px+py*py);
							else if (brushshape == BrushShape.SQUARE)
								r = Math.max(Math.abs(px), Math.abs(py));
							under_brush[i][j] = (r <= brushsize);
						}
						else {
							under_brush[i][j] = false;
						}
					}
				}
			}

			break;

		case INTERACT:
			if (materials[mx_index][my_index].type == MaterialType.EMF || materials[mx_index][my_index].type == MaterialType.SWITCH)
				r.setCursor(HAND_CURSOR);
			else
				r.setCursor(DEFAULT_CURSOR);
			
			if (pressing) {
				boolean turn_on_EMF = !selected_EMF[mx_index][my_index];
				
				for (int i = 0; i < nx; i++)
				{
					for (int j = 0; j < ny; j++)
					{
						selected_EMF[i][j] = false;
					}
				}
				EMF_selected = false;

				if (materials[mx_index][my_index].type == MaterialType.EMF && turn_on_EMF) {
					floodFillSelectEMF(mx_index, my_index, true);
					EMF_selected = true;
					int setting = (int)(Math.round(50*materials[mx_index][my_index].emf/max_EMF));
					opts.gui_parameter3.setValue(setting);
					prev_EMF_setting = setting;
				}
				
				if (materials[mx_index][my_index].type == MaterialType.SWITCH) {
					this.floodFillToggleSwitch(mx_index, my_index, 1-materials[mx_index][my_index].activated);
					update = true;
				}
			}
			break;
		case FLOODSELECT:
		case SELECT:
			if (pressing) {
				if (brush == Brush.FLOODSELECT && !moving_selection) {
					floodFillSelect(mx_index, my_index, materials[mx_index][my_index].type, !selected[mx_index][my_index]);
				}
				else if (moving_selection && !dragging_selection) {
					for (int i = 0; i < nx; i++)
					{
						for (int j = 0; j < ny; j++)
						{
							int si = i-delta_mx_index;
							int sj = j-delta_my_index;
							if (si >= 0 && sj >= 0 && si < nx && sj < ny && selection[si][sj].type != MaterialType.VACUUM) {
								materials[i][j].erase();
								materials[i][j] = selection[si][sj].clone();
								selected[i][j] = true;
							}
						}
					}
					moving_selection = false;
					dragging_selection = false;
				} else if (!moving_selection && selected[mx_index][my_index]) {
					for (int i = 0; i < nx; i++)
					{
						for (int j = 0; j < ny; j++)
						{
							selection[i][j].erase();
							if (selected[i][j]) {
								selection[i][j] = materials[i][j].clone();
								selected[i][j] = false;
								materials[i][j].erase();
							}
						}
					}
					moving_selection = true;
					dragging_selection = true;
					delta_mx_index = 0;
					delta_my_index = 0;
				}
			} else if (mouse_pressed) {
				if (brush != Brush.FLOODSELECT) {
					if (dragging_selection) {
						delta_mx_index = mx_index - mx_start_index;
						delta_my_index = my_index - my_start_index;
					} else {
						int mx0 = Math.min(mx_start_index, mx_index);
						int my0 = Math.min(my_start_index, my_index);
						int mx1 = Math.max(mx_start_index, mx_index);
						int my1 = Math.max(my_start_index, my_index);

						for (int i = 0; i < nx; i++)
						{
							for (int j = 0; j < ny; j++)
							{
								if (i >= mx0 && i <= mx1 && j >= my0 && j <= my1)
									selected[i][j] = true;
								else
									selected[i][j] = false;
							}
						}
					}
				}
			} else if (releasing) {
				if (brush != Brush.FLOODSELECT) {
					delta_mx_index = mx_index - mx_start_index;
					delta_my_index = my_index - my_start_index;
					if (dragging_selection) {
						for (int i = 0; i < nx; i++)
						{
							for (int j = 0; j < ny; j++)
							{
								int si = i-delta_mx_index;
								int sj = j-delta_my_index;
								if (si >= 0 && sj >= 0 && si < nx && sj < ny && selection[si][sj].type != MaterialType.VACUUM) {
									materials[i][j].erase();
									materials[i][j] = selection[si][sj].clone();
									selected[i][j] = true;
								}
							}
						}
						moving_selection = false;
						dragging_selection = false;
					} else {
						if (delta_mx_index == 0 && delta_my_index == 0) {
							for (int i = 0; i < nx; i++)
							{
								for (int j = 0; j < ny; j++)
								{
									selected[i][j] = false;
								}
							}
						}
					}
				}
			} else {
				delta_mx_index = mx_index;
				delta_my_index = my_index;
			}
			break;
		case CURRENT:
			if (pressing) {
				CurrentProbe p = new CurrentProbe();
				p.x1 = mx_start_index;
				p.y1 = my_start_index;
				p.x2 = mx_index;
				p.y2 = my_index;
				currentprobes.add(p);
			} else if (mouse_pressed) {
				currentprobes.get(currentprobes.size()-1).x2 = mx_index;
				currentprobes.get(currentprobes.size()-1).y2 = my_index;
			}
			break;
		case VOLTAGE:
			if (pressing) {
				VoltageProbe p = new VoltageProbe();
				p.x = mx_start_index;
				p.y = my_start_index;
				voltageprobes.add(p);
			} else if (mouse_pressed) {
				voltageprobes.get(voltageprobes.size()-1).x = mx_index;
				voltageprobes.get(voltageprobes.size()-1).y = my_index;
			}
			break;
		case DELETEPROBE:
			r.setCursor(DEFAULT_CURSOR);
			int i = 0;
			while(i < voltageprobes.size()) {
				VoltageProbe p = voltageprobes.get(i);
				if (length(p.x-mx_index, p.y-my_index) < 3) {
					r.setCursor(HAND_CURSOR);
					if (pressing) {
						voltageprobes.remove(i);
						i--;
					}
				}
				i++;
			}
			
			i = 0;
			while(i < currentprobes.size()) {
				CurrentProbe p = currentprobes.get(i);
				if (length(p.x1-mx_index, p.y1-my_index) < 3 || length(p.x2-mx_index, p.y2-my_index) < 3) {
					r.setCursor(HAND_CURSOR);
					if (pressing) {
						currentprobes.remove(i);
						i--;
					}
				}
				i++;
			}
			
			if (ground != null && length(ground.x-mx_index, ground.y-my_index) < 3) {
				r.setCursor(HAND_CURSOR);
				if (pressing) {
					ground = null;
				}
			}
			
			break;
		case TEXT:
			r.setCursor(HAND_CURSOR);
			if (mouse_pressed) {
				texting = true;
				text_x = mx_index;
				text_y = my_index;
			}
			break;
		case GROUND:
			if (pressing) {
				if (ground == null)
					ground = new VoltageProbe();
				ground.x = mx_start_index;
				ground.y = my_start_index;
			} else if (mouse_pressed) {
				ground.x = mx_index;
				ground.y = my_index;
			}
			break;
		}
		
		if (brush != Brush.TEXT)
		{
			texting = false;
		}

		setEMFs();
		
		if (releasing || (BoundaryCondition)opts.gui_bc.getSelectedItem() != prev_boundary || update) {

			constructBoundary();
			updateAllMaterials();
			multigridSolve(true, false);
		}
		
		prev_boundary = (BoundaryCondition)opts.gui_bc.getSelectedItem();

		mxp_realspace = mx_realspace;
		myp_realspace = my_realspace;
	}
	
	public void drawMaterialLine(double x1, double y1, double x2, double y2, Brush brush, BrushShape brushshape, MaterialType mat, double brushsize, double EMF_angle) {
		Vector a = new Vector(0, 0);
		Vector b = new Vector(0, 0);
		Vector p = new Vector(0, 0);
		Vector ab = new Vector(0, 0);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				double cx = 0;
				double cy = 0;
				cx = i*ds;
				cy = j*ds;

				a.initialize(x1, y1);
				b.initialize(x2, y2);
				p.initialize(cx, cy);
				p.addmult(a, -1);
				ab.copy(b);
				ab.addmult(a, -1);

				double l2 = ab.dot(ab);
				if (l2 == 0)
					l2 = 1;
				double t = clamp(p.dot(ab)/l2, 0, 1);
				ab.scalarmult(t);
				p.addmult(ab, -1);
				double r = 0;
				if (brushshape == BrushShape.CIRCLE)
					r = Math.sqrt(p.dot(p));
				else if (brushshape == BrushShape.SQUARE)
					r = Math.max(Math.abs(p.x), Math.abs(p.y));
				if (r <= brushsize) {
					if (mat == MaterialType.VACUUM) {
						materials[i][j].erase();
					} else if (materials[i][j].type == MaterialType.VACUUM || brush == Brush.REPLACE) {
						materials[i][j].erase();
						initializeMaterial(i, j, mat);
						if (mat == MaterialType.EMF) materials[i][j].emf_direction = EMF_angle;
					}
				}
			}
		}
	}
	
	public void setEMFs() {
		int EMF_setting = opts.gui_parameter3.getValue();
		double new_EMF = max_EMF*EMF_setting/50.0;
		
		if (opts.gui_brush.getSelectedItem() == Brush.INTERACT && EMF_selected) {
			opts.gui_parameter3.setVisible(true);
			opts.gui_parameter3_text.setVisible(true);
			opts.gui_parameter3_text.setText("EMF: " + getSI(new_EMF, "V/m"));
		} else {
			opts.gui_parameter3.setVisible(false);
			opts.gui_parameter3_text.setVisible(false);
			opts.gui_parameter3_text.setText("");
		}
		
		if (EMF_setting != prev_EMF_setting && EMF_selected) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					if (materials[i][j].type == MaterialType.EMF && selected_EMF[i][j]) {
						materials[i][j].emf = new_EMF;
					}
				}
			}

			updateJustEMFs();
		}
		
		prev_EMF_setting = EMF_setting;
	}


	public void checkCFL() {
		System.out.println("Wave equation CFL ratio = " + (Math.sqrt(2)*c)/(ds/dt_maximum));
		System.out.println("Electron diffusion CFL ratio = " + dt_maximum/(ds*ds/(4*D_electron)));
		System.out.println("Hole diffusion CFL ratio = " + dt_maximum/(ds*ds/(4*D_hole)));
	}

	public void calcMiscFields(boolean updatePhi) {
		ScalarView view_scalar = (ScalarView) opts.gui_view.getSelectedItem();
		
		//Always find potentials
		if (updatePhi)
			multigridSolve(false, true);
		
		t8.start();
		
		
		
		sign_violation = false;
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
					F_n[i][j] = F0_n[i][j] + kT*Math.log(rho_n[i][j]/q_n);
					F_p[i][j] = F0_p[i][j] + kT*Math.log(rho_p[i][j]/q_p);
					// Think about why we should take average
					

					double probe_n = -e_charge*ni_metal;
					double probe_p = e_charge*ni_metal;
					//double sigma_n = mu_electron*logmean(-rho_n[i][j], -probe_n);
					//double sigma_p = mu_hole*logmean(rho_p[i][j], probe_p);
					double sigma_n = mu_electron*(-rho_n[i][j]);
					double sigma_p = mu_hole*rho_p[i][j];
					
					F[i][j] = (sigma_n*(F_n[i][j]/q_n+phi[i][j]) + sigma_p*(F_p[i][j]/q_p+phi[i][j]))
							/(sigma_n + sigma_p) - W_semi/eVtoJ;
				}
				
				G[i][j] = conducting[i][j]*R[i][j]*(K[i][j] - rho_n[i][j]*rho_p[i][j]/(q_n*q_p));
				
				if (rho_n[i][j] > error_detection_threshold || rho_p[i][j] < -error_detection_threshold)
					sign_violation = true;
			}
		}

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				Jx_free[i][j] = Jx_abs[i][j] + Jx_n[i][j] + Jx_p[i][j];
				Dx[i][j] = epsx[i][j]*Ex[i][j];
				Sy[i][j] = -Ex[i][j]*0.5*(Hz[i][j] + Hz[i][j-1]);
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				Jy_free[i][j] = Jy_abs[i][j] + Jy_n[i][j] + Jy_p[i][j];
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
					double n_contrib = -G[i][j]*E0_n[i][j];
					double jn_contrib = 0.5*(grad_E0x_n[i-1][j]*Jx_n[i-1][j]+grad_E0x_n[i][j]*Jx_n[i][j]
							+grad_E0y_n[i][j-1]*Jy_n[i][j-1]+grad_E0y_n[i][j]*Jy_n[i][j])/q_n;
					
					double p_contrib = -G[i][j]*E0_p[i][j];
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
					double n_contrib = -G[i][j]*F_n[i][j];
					double jn_contrib = 0.5*(grad_Fx_n[i-1][j]*Jx_n[i-1][j]+grad_Fx_n[i][j]*Jx_n[i][j]
							+grad_Fy_n[i][j-1]*Jy_n[i][j-1]+grad_Fy_n[i][j]*Jy_n[i][j])/q_n;
					
					double p_contrib = -G[i][j]*F_p[i][j];
					double jp_contrib = 0.5*(grad_Fx_p[i-1][j]*Jx_p[i-1][j]+grad_Fx_p[i][j]*Jx_p[i][j]
							+grad_Fy_p[i][j-1]*Jy_p[i][j-1]+grad_Fy_p[i][j]*Jy_p[i][j])/q_p;

					double ohm_contrib = 0.5*(Ex[i-1][j]*Jx_free[i-1][j]+Ex[i][j]*Jx_free[i][j]+Ey[i][j-1]*Jy_free[i][j-1]+Ey[i][j]*Jy_free[i][j]);

					S[i][j] = n_contrib + jn_contrib + p_contrib + jp_contrib + ohm_contrib;
				}
			}
		}
		
		for (VoltageProbe p: voltageprobes) {
			p.potential = F[p.x][p.y];
		}
		
		if (ground != null)
			ground.potential = F[ground.x][ground.y];
		
		for (CurrentProbe p: currentprobes) {
			p.current = calcCurrent(p.x1, p.y1, p.x2, p.y2);
		}

		updateMiscFields = false;
		
		t8.stop();
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
			return Math.abs(Jy[x0+1][y0])*ds;
		} else if (dx == -1 && dy == 0) {
			return Math.abs(Jy[x0][y0])*ds;
		} else if (dx == 0 && dy == 1) {
			return Math.abs(Jx[x0][y0+1])*ds;
		} else if (dx == 0 && dy == -1) {
			return Math.abs(Jx[x0][y0])*ds;
		} else if (dx == 1 && dy == 1) {
			return (Math.abs(Jx[x0][y0+1]) + Math.abs(Jy[x0+1][y0+1]))*ds;
		} else if (dx == 1 && dy == -1) {
			return (Math.abs(Jx[x0][y0])  + Math.abs(Jy[x0+1][y0-1]))*ds;
		} else if (dx == -1 && dy == 1) {
			return (Math.abs(Jx[x0][y0+1]) + Math.abs(Jy[x0][y0+1]))*ds;
		} else if (dx == -1 && dy == -1) {
			return (Math.abs(Jx[x0][y0]) + Math.abs(Jy[x0][y0-1]))*ds;
		}
		return 0;
	}
	
	public void prescaleDielectric() {
		downscale_x_vector(epsx, MG_epsx, log2_resolution);
		downscale_y_vector(epsy, MG_epsy, log2_resolution);
	}
	
	public void multigridSolve(boolean correctEfield, boolean computePhi) {
		assert(!(correctEfield && computePhi));
		
		boolean debug = false;
		
		if (correctEfield)
			t4.start();
		if (computePhi)
			t7.start();
			
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				MG_phi1[i][j] = 0;
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
		
		downscale(MG_rho0, MG_rho, log2_resolution);
		
		for (int fineness = 2; fineness <= log2_resolution; fineness++) {
			int nx_tmp = (1 << fineness);
			int ny_tmp = (1 << fineness);
			double gridsize = this.default_width/(1 << fineness);
			
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
			
			/*for (int i = 1; i < 2*nx_tmp-1; i++)
			{
				for (int j = 1; j < 2*nx_tmp-1; j++) {
					MG_phi1[i][j] = this.bilinearinterp(MG_phi2, (i-0.5)/2.0, (j-0.5)/2.0);
				}
			}*/
			
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
	
	public double bilinearinterp(double[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-2 || Math.abs(y-Math.round(y)) < 1e-2) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= nx) i = nx - 1;
			if (j >= ny) j = ny - 1;
			return array[i][j];
		}
		
		if (xfloor < 0) {
			xfloor = 0;
			fx = 0.0;
		} else if (xfloor >= nx - 1) {
			xfloor = nx - 2;
			fx = 1.0;
		}
		if (yfloor < 0) {
			yfloor = 0;
			fy = 0.0;
		} else if (yfloor >= ny - 1) {
			yfloor = ny - 2;
			fy = 1.0;
		}
		double va = array[xfloor][yfloor]*(1.0-fx) + array[xfloor+1][yfloor]*fx;
		double vb = array[xfloor][yfloor+1]*(1.0-fx) + array[xfloor+1][yfloor+1]*fx;
		
		return va*(1.0-fy) + vb*fy;
	}
	
	class FloodFillCoordinate {
		int i;
		int j;
		
		public FloodFillCoordinate(int i, int j) {
			this.i = i;
			this.j = j;
		}
	}

	public void floodFillSet(int i, int j, MaterialType old_mat, MaterialType new_mat, double EMF_angle) {
		if (old_mat == new_mat)
			return;
		
		Queue<FloodFillCoordinate> queue = new LinkedList<FloodFillCoordinate>();
		queue.add(new FloodFillCoordinate(i, j));
		
		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < nx && coord.j >= 0 && coord.j < ny && materials[coord.i][coord.j].type == old_mat && materials[coord.i][coord.j].type != new_mat) {
				materials[coord.i][coord.j].erase();
				initializeMaterial(coord.i, coord.j, new_mat);
				if (new_mat == MaterialType.EMF) materials[coord.i][coord.j].emf_direction = EMF_angle;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}
	
	public void floodFillSelect(int i, int j, MaterialType mat, boolean select) {
		Queue<FloodFillCoordinate> queue = new LinkedList<FloodFillCoordinate>();
		queue.add(new FloodFillCoordinate(i, j));
		
		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < nx && coord.j >= 0 && coord.j < ny && materials[coord.i][coord.j].type == mat && selected[coord.i][coord.j] != select) {
				selected[coord.i][coord.j] = select;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}
	
	public void floodFillSelectEMF(int i, int j, boolean select) {
		Queue<FloodFillCoordinate> queue = new LinkedList<FloodFillCoordinate>();
		queue.add(new FloodFillCoordinate(i, j));
		
		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < nx && coord.j >= 0 && coord.j < ny && materials[coord.i][coord.j].type == MaterialType.EMF && selected_EMF[coord.i][coord.j] != select) {
				selected_EMF[coord.i][coord.j] = select;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}
	
	public void floodFillToggleSwitch(int i, int j, int active) {
		Queue<FloodFillCoordinate> queue = new LinkedList<FloodFillCoordinate>();
		queue.add(new FloodFillCoordinate(i, j));
		
		while (queue.size() > 0) {
			FloodFillCoordinate coord = queue.remove();
			if (coord.i >= 0 && coord.i < nx && coord.j >= 0 && coord.j < ny && materials[coord.i][coord.j].type == MaterialType.SWITCH && materials[coord.i][coord.j].activated != active) {
				materials[coord.i][coord.j].activated = active;
				queue.add(new FloodFillCoordinate(coord.i-1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i+1, coord.j));
				queue.add(new FloodFillCoordinate(coord.i, coord.j-1));
				queue.add(new FloodFillCoordinate(coord.i, coord.j+1));
			}
		}
	}
	
	public void drawPixelRectangle(int x, int y, int w, int h) {
		for (int i = x; i < x+w; i++) {
			for (int j = y; j < y+h; j++) {
				setPixel(i, j);
			}
		}
	}
	
	public void setPixel(int i, int j) {
		if (i < 0 || j < 0 || i >= nx || j >= ny)
			return;
		
		image_r[i][j] = (float)(image_r[i][j]*alphaBG + col_r*alphaFG);
		image_g[i][j] = (float)(image_g[i][j]*alphaBG + col_g*alphaFG);
		image_b[i][j] = (float)(image_b[i][j]*alphaBG + col_b*alphaFG);
	}
	
	public void stampPixelData() {
		t9.start();
		int scansize = nx*scalefactor;
		int[] imgData = ((DataBufferInt)screen.getRaster().getDataBuffer()).getData();
		for (int x = 0; x < nx*scalefactor; x++) {
			for (int y = 0; y < ny*scalefactor; y++) {
				int i = x/scalefactor;
				int j = y/scalefactor;
				double scale = 1f/max(image_r[i][j], image_g[i][j], image_b[i][j], 1f);
				int rgb = clamp((int)(256*image_r[i][j]*scale), 0, 255) << 16
						| clamp((int)(256*image_g[i][j]*scale), 0, 255) << 8
						| clamp((int)(256*image_b[i][j]*scale), 0, 255);
				imgData[x + y*scansize] = rgb;
			}
		}
		t9.stop();
	}
	
	public void drawLine(int x0, int y0, int x1, int y1, boolean draw_starting_point) {
		try {
			int[] imgData = ((DataBufferInt)screen.getRaster().getDataBuffer()).getData();
			int scansize = scalefactor*nx;
			
			int dy = y1 - y0;
			int dx = x1 - x0;
			float t = (float) 0.5;

			int rgb;
			int r;
			int g;
			int b;
			
			if (draw_starting_point) {
				rgb = imgData[x0+y0*scansize];
				r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
				g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
				b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
				if (r > 255)
					r = 255;
				if (g > 255)
					g = 255;
				if (b > 255)
					b = 255;
				imgData[x0+y0*scansize] =  255<<24 | r << 16 | g << 8 | b;
			}


			if (Math.abs(dx) > Math.abs(dy)) {
				float m = (float) dy / (float) dx;
				t += y0;
				dx = (dx < 0) ? -1 : 1;
				m *= dx;
				while (x0 != x1) {
					x0 += dx;
					t += m;

					rgb =  imgData[x0+((int)t)*scansize];
					r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
					g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
					b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					imgData[x0+((int)t)*scansize] = 255<<24 | r << 16 | g << 8 | b;

				}
			} else {
				float m = (float) dx / (float) dy;
				t += x0;
				dy = (dy < 0) ? -1 : 1;
				m *= dy;
				while (y0 != y1) {
					y0 += dy;
					t += m;

					rgb =  imgData[((int)t)+y0*scansize];
					r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
					g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
					b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					imgData[((int)t)+y0*scansize] = 255<<24 | r << 16 | g << 8 | b;
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			return;
		}
	}

	public void drawPixelLine(int x0, int y0, int x1, int y1) {
		int dy = y1 - y0;
		int dx = x1 - x0;
		float t = (float) 0.5;

		setPixel(x0, y0);

		if (Math.abs(dx) > Math.abs(dy)) {
			float m = (float) dy / (float) dx;
			t += y0;
			dx = (dx < 0) ? -1 : 1;
			m *= dx;
			while (x0 != x1) {
				x0 += dx;
				t += m;

				setPixel(x0, (int)t);
			}
		} else {
			float m = (float) dx / (float) dy;
			t += x0;
			dy = (dy < 0) ? -1 : 1;
			m *= dy;
			while (y0 != y1) {
				y0 += dy;
				t += m;

				setPixel((int)t, y0);
			}
		}
	}
	
	public void setColor(int r, int g, int b) {
		col_r = r/255f;
		col_g = g/255f;
		col_b = b/255f;
	}
	
	public void setColorFloat(float r, float g, float b) {
		if (!(r+b+g < Float.MAX_VALUE)) {
			col_r = 0;
			col_g = 0;
			col_b = 0;
			return;
		}
		float scale = 1f/max(r, g, b, 1f);
		col_r = Math.max(r*scale, 0);
		col_g = Math.max(g*scale, 0);
		col_b = Math.max(b*scale, 0);
	}
	
	public float max(float x, float y, float z) {
		return Math.max(Math.max(x, y), z);
	}
	
	public float max(float x, float y, float z, float w) {
		return Math.max(Math.max(Math.max(x, y), z), w);
	}
	
	public float min(float x, float y, float z) {
		return Math.min(Math.min(x, y), z);
	}
	
	public float min(float x, float y, float z, float w) {
		return Math.min(Math.min(Math.min(x, y), z), w);
	}
	
	public int clamp(int val, int min, int max) {
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}

	public double clamp(double val, double min, double max) {
		if (val != val) return min;
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}

	public void setalphaBG(double alpha) {
		alphaBG = alpha;
	}
	
	public void setalphaFG(double alpha) {
		alphaFG = alpha;
	}
	
	public void render() {
		t5.start();
		Graphics2D g = (Graphics2D) screen.getGraphics();
		double scalingconstant = 10.0*Math.pow(10.0, opts.gui_brightness.getValue()/10.0);


		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				image_r[i][j] = 0;
				image_g[i][j] = 0;
				image_b[i][j] = 0;
			}
		}
		
		clearStrings();
		
		/* Draw pixels */
		
		setalphaBG(0);
		setalphaFG(1);

		Brush brush = (Brush) opts.gui_brush.getSelectedItem();

		if (opts.gui_elem_colors.isSelected()) {
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					setColor(materials[i][j].type.color_r, materials[i][j].type.color_g, materials[i][j].type.color_b);
					
					setPixel(i, j);
				}
			}
		} else {
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					setColor(materials[i][j].type.color_grayscale, materials[i][j].type.color_grayscale, materials[i][j].type.color_grayscale);
					setPixel(i, j);
				}
			}
		}
		
		if (moving_selection) {
			if (opts.gui_elem_colors.isSelected()) {
				for (int i = 0; i < nx; i++) {
					for (int j = 0; j < ny; j++) {
						int si = i-delta_mx_index;
						int sj = j-delta_my_index;
						if (si >= 0 && sj >= 0 && si < nx && sj < ny && selection[si][sj].type != MaterialType.VACUUM) {
							setColor(selection[si][sj].type.color_r, selection[si][sj].type.color_g, selection[si][sj].type.color_b);
							setPixel(i, j);
						}
					}
				}
			}else {
				for (int i = 0; i < nx; i++) {
					for (int j = 0; j < ny; j++) {
						int si = i-delta_mx_index;
						int sj = j-delta_my_index;
						if (si >= 0 && sj >= 0 && si < nx && sj < ny && selection[si][sj].type != MaterialType.VACUUM) {
							setColor(selection[si][sj].type.color_grayscale, selection[si][sj].type.color_grayscale, selection[si][sj].type.color_grayscale);
							setPixel(i, j);
						}
					}
				}
			}
		}

		if ((ScalarView) opts.gui_view.getSelectedItem() != ScalarView.NONE) {
			setalphaBG(1.0);
			setalphaFG(1.0);
			for (int i = 1; i < nx-1; i++) {
				for (int j = 1; j < ny-1; j++) {

					switch ((ScalarView) opts.gui_view.getSelectedItem()) {
					case NONE:
						setColorFloat(0, 0, 0);
						break;
					case B_FIELD:
						float v;
						v = (float)(parity*0.25*(Bz[i][j]+Bz[i-1][j]+Bz[i][j-1]+Bz[i-1][j-1])*1e5*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case E_FIELD:
						double vx = 0.5*(Ex[i][j]+Ex[i][j+1])*scalingconstant/1e5;
						double vy = 0.5*(Ey[i][j]+Ey[i+1][j])*scalingconstant/1e5;
						setColorFloat(0, (float)(length(vx, vy)), 0);
						break;
					case H_FIELD:
						v = (float)(parity*Hz[i][j]*mu0*1e5*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case CURRENT:
						vx = 0.5*(Jx_free[i][j]+Jx_free[i-1][j])*scalingconstant;
						vy = 0.5*(Jy_free[i][j]+Jy_free[i][j-1])*scalingconstant;
						setColorFloat(0, (float)(length(vx, vy)/1e7), 0);
						break;
					case POTENTIAL:
						v = (float)(phi[i][j]*scalingconstant);
						setColorFloat(v, 0, -v);
						break;
					case CHARGE:
						v = (float)(rho_free[i][j]*scalingconstant);
						float rc = Math.min(Math.max(v, 0), 1);
						float bc = Math.min(Math.max(-v, 0), 1);
						float gc = Math.min(rc, bc);

						setalphaBG(1-0.5*gc);
						//setalphaFG(gc);
						setColorFloat(rc, gc, bc);
						break;
					case BACKGROUND_CHARGE:
						v = (float)(rho_back[i][j]*scalingconstant);
						setColorFloat(v, 0, -v);
						break;
					case ELECTRON_CHARGE:
						v = (float)(rho_n[i][j]*scalingconstant);
						setColorFloat(v, 0, -v);
						break;
					case HOLE_CHARGE:
						v = (float)(rho_p[i][j]*scalingconstant);
						setColorFloat(v, 0, -v);
						break;
					case COMBINED_CHARGE:
						rc = (float)(clamp(0.2*Math.log(rho_p[i][j]*scalingconstant), 0, 1));
						bc = (float)(clamp(0.2*Math.log(-rho_n[i][j]*scalingconstant), 0, 1));
						gc = Math.min(rc, bc);
						double it = Math.max(rc, bc);
						setalphaBG(1-0.5*gc);
						setalphaFG(0.6*it);
						setColorFloat(rc, gc, bc);
						break;
					case ENERGY:
						v = (float)(u[i][j]*scalingconstant);
						setColorFloat(0, v, 0);
						break;
					case ENTROPY:
						v = (float)(S[i][j]*scalingconstant/1e12);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case HEAT:
						v = (float)(Q[i][j]*scalingconstant/1e12);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case ELECTRON_POTENTIAL:
						v = (float)(conducting[i][j]*(F_n[i][j]/q_n+phi[i][j]-W_semi/eVtoJ)*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case HOLE_POTENTIAL:
						v = (float)(conducting[i][j]*(F_p[i][j]/q_p+phi[i][j]-W_semi/eVtoJ)*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case DEBUG:
						v = (float)((visited[i][j]? 1:0)*scalingconstant);
						setColorFloat(v, 0, -v);
						break;
					case RECOMBINATION:
						v = (float)(-G[i][j]/1e30*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case AVERAGE_POTENTIAL:
						v = (float)(F[i][j]*scalingconstant);
						setColorFloat(v, Math.abs(v), -v);
						break;
					case LIGHT:
						v = (float)(-materials[i][j].semiconducting*this.G[i][j]/1e30*scalingconstant);
						setColorFloat(v, v, v);
						break;
					}
					setPixel(i, j);
				}
			}
		}
		
		boolean highlight = (opts.gui_brush_highlight.isSelected() && Brush.isBrushShapeImportant(brush));
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (materials[i][j].type == MaterialType.EMF && i > 0 && j > 0 && i < nx-1 && j < ny-1) {
					setalphaBG(0.25);
					setalphaFG(0.75);
					
					int offset = 0;
					if (selected_EMF[i][j])
						offset = 60*(2*((i+j)%2)-1);
					
					if ((materials[i+1][j].type != MaterialType.EMF
							|| materials[i-1][j].type != MaterialType.EMF
							|| materials[i][j+1].type != MaterialType.EMF
							| materials[i][j-1].type != MaterialType.EMF))
					{
						offset = -30;
					}
					
					int delta_r = MaterialType.EMF.color_r+offset;
					int delta_g = MaterialType.EMF.color_g+offset;
					int delta_b = MaterialType.EMF.color_b+offset;
					setColor(delta_r, delta_g, delta_b);
					
					setPixel(i, j);
				} else if (!(materials[i][j].activated == 1)) {
					setalphaBG(0.25);
					setalphaFG(0.75);
					setColor(0, 0, 0);
					setPixel(i, j);
				}
				
				if (selected[i][j] || (highlight && under_brush[i][j]))
				{
					setalphaBG(0.75);
					setalphaFG(0.25);
					
					int s = selected[i][j]? 1:0;
					int h = (highlight && under_brush[i][j])? 1:0;
					int delta_r = 256*s + 256*h;
					int delta_g = 100*s + 256*h;
					int delta_b = 256*s + 256*h;
					
					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				}

			}
		}
		
		
		setalphaBG(1);
		setalphaFG(1);
		setColorFloat(0.7f, 0.7f, 0.7f);
		
		if (brush == Brush.LINE && mouse_pressed) {
			drawPixelLine(mx_start_index, my_start_index, mx_index, my_index);
		}


		setalphaFG(1.0);
		setColorFloat(0.5f, 1.0f, 1.0f);

		for (CurrentProbe p: currentprobes) {
			drawPixelRectangle(p.x1-1, p.y1-1, 3, 3);
			drawPixelRectangle(p.x2-1, p.y2-1, 3, 3);
		}
		
		for (VoltageProbe p: voltageprobes) {
			drawPixelRectangle(p.x-1, p.y-1, 3, 3);
		}
		
		if (ground != null) {
			drawPixelRectangle(ground.x-1, ground.y-1, 3, 3);
		}
		
		setalphaFG(0.3);
		setColorFloat(0.5f, 1.0f, 1.0f);
		
		for (CurrentProbe p: currentprobes) {
			drawPixelLine(p.x1, p.y1, p.x2, p.y2);
		}
		
		setalphaFG(0.8);
		setColorFloat(1.0f, 1.0f, 1.0f);
		
		if (texting) {
			drawPixelLine(text_x, text_y, text_x, text_y+7);
		}
				
		stampPixelData();
		
		
		/* Draw vectors */
		
		if ((VectorView) opts.gui_view_vec.getSelectedItem() != VectorView.NONE) {
			setalphaBG(1.0);
			
			double arrowlength = 10.0/scalefactor;

			double vectorscalingconstant = 0.01*Math.pow(10.0, opts.gui_brightness_vec.getValue()/5.0);
			
			VectorMode vector_display_mode = (VectorMode)opts.gui_view_vec_mode.getSelectedItem();
			
			rand.setSeed(4);
			
			int density = 75;

			double randomness = 0;
			
			if (vector_display_mode == VectorMode.ARROWS) {
				randomness = 0.5;
			} else if (vector_display_mode == VectorMode.LINES) {
				randomness = 0.75;
			}

			double[][] vf_x = null;
			double[][] vf_y = null;
			
			switch ((VectorView) opts.gui_view_vec.getSelectedItem()) {
			case NONE:
				break;
			case D_FIELD:
				vf_x = Dx;
				vf_y = Dy;
				break;
			case E_FIELD:
				vf_x = Ex;
				vf_y = Ey;
				break;
			case ELECTRON_CURRENT:
				vf_x = Jx_n;
				vf_y = Jy_n;
				break;
			case HOLE_CURRENT:
				vf_x = Jx_p;
				vf_y = Jy_p;
				break;
			case TOTAL_CURRENT:
				vf_x = Jx_free;
				vf_y = Jy_free;
				break;
			case POYNTING:
				vf_x = Sx;
				vf_y = Sy;
				break;
			case EMF:
				vf_x = emfx;
				vf_y = emfy;
				break;
			}
			
			
			Vector ctr = new Vector(0,0);
			Vector arrow = new Vector(0,0);
			Vector tip1 = new Vector(0,0);
			Vector tip2 = new Vector(0,0);
			Vector body1 = new Vector(0,0);
			Vector body2 = new Vector(0,0);
			
			for (int i = 0; i < density; i++) {
				for (int j = 0; j < density; j++) {

					//double x = (nx-1)*(i+0.5)/50;
					//double y = (ny-1)*(j+0.5)/50;
					double x = nx*(i+randomness*(rand.nextFloat()-0.5))/density;
					double y = ny*(j+randomness*(rand.nextFloat()-0.5))/density;
					
					
					if (vector_display_mode == VectorMode.LINES && vf_x != null) {
						for (int sign = -1; sign <= 1; sign += 2) {
							
							double prevx = x;
							double prevy = y;
							double dx = 0;
							double dy = 0;
							
							for (int k = 0; k < 10; k++) {
									dx = bilinearinterp(vf_x, prevx-0.5, prevy);
									dy = bilinearinterp(vf_y, prevx, prevy-0.5);

								double fieldmagnitude = Math.sqrt(dx*dx+dy*dy);
								setalphaFG(0.1*Math.sqrt(1/(0.1*k*k+1.0)*vectorscalingconstant*fieldmagnitude));
								
								if (fieldmagnitude != 0) {
									dx /= fieldmagnitude;
									dy /= fieldmagnitude;
								}
								
								//setcol(fieldmagnitude, fieldmagnitude, fieldmagnitude, 30);
								setColorFloat(1.0f, 1.0f, 1.0f);

								double nextx = prevx + dx*arrowlength*0.25*sign;
								double nexty = prevy + dy*arrowlength*0.25*sign;

								drawLine((int)((prevx+0.5)*scalefactor), (int)((prevy+0.5)*scalefactor), (int)((nextx+0.5)*scalefactor), (int)((nexty+0.5)*scalefactor), k == 0 && sign == 1);

								prevx = nextx;
								prevy = nexty;
							}
						}
					} else if (vector_display_mode == VectorMode.ARROWS && vf_x != null) {
						ctr.x = x+0.5;
						ctr.y = y+0.5;

						arrow.x = bilinearinterp(vf_x,x-0.5, y);
						arrow.y = bilinearinterp(vf_y,x, y-0.5);

						double fieldmagnitude = Math.max(0.1, vectorscalingconstant*Math.sqrt(arrow.dot(arrow)));
						arrow.normalize();
						tip1.copy(arrow);
						tip2.copy(arrow);
						tip1.rotate(Math.PI*5.0/6.0);
						tip2.rotate(Math.PI*7.0/6.0);
						
						body1.copy(ctr);
						body1.addmult(arrow, -0.5*arrowlength);
						body2.copy(ctr);
						body2.addmult(arrow, 0.5*arrowlength);
						tip1.scalarmult(0.35*arrowlength);
						tip1.add(body2);
						tip2.scalarmult(0.35*arrowlength);
						tip2.add(body2);
						setColorFloat((float)fieldmagnitude, (float)fieldmagnitude, (float)fieldmagnitude);
						setalphaFG(0.1*Math.sqrt(fieldmagnitude));
						drawLine((int)(body1.x*scalefactor), (int)(body1.y*scalefactor), (int)(body2.x*scalefactor), (int)(body2.y*scalefactor), true);
						drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip1.x*scalefactor), (int)(tip1.y*scalefactor), false);
						drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip2.x*scalefactor), (int)(tip2.y*scalefactor), false);
					}
				}
			}
		}

		/* Draw text */

		((Graphics2D)g).setRenderingHint(
				RenderingHints.KEY_TEXT_ANTIALIASING,
				RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		
		for (VoltageProbe p: voltageprobes) {
			if (ground != null)
				drawStringWithBackgroundAndBorder("V = " + getSI(p.potential - ground.potential, "V"), p.x*scalefactor-5, p.y*scalefactor - 12, g);
			else
				drawStringWithBackgroundAndBorder("V = " + getSI(p.potential, "V"), p.x*scalefactor-5, p.y*scalefactor - 12, g);
		}

		if (ground != null)
			drawStringWithBackgroundAndBorder("Ground = " + getSI(ground.potential - ground.potential, "V"), ground.x*scalefactor-5, ground.y*scalefactor - 12, g);
		
		for (CurrentProbe p: currentprobes) {
			double xa = 0.5*(p.x1+p.x2)*scalefactor;
			double ya = 0.5*(p.y1+p.y2)*scalefactor;
			
			double dx = p.x2 - p.x1;
			double dy = p.y2 - p.y1;
			double len = length(dx, dy);
			dx = dx/len;
			dy = dy/len;
			if (Math.abs(dx) > Math.abs(dy))
			{
				dx = -Math.abs(dx);
			} else {
				dy = -2*Math.abs(dy);
			}
			
			drawStringWithBackgroundAndBorder("I = " + getSI(p.current*depth, "A"), (int)(xa-8*dy)-5, (int)(ya+12*dx)+5, g);
		}
		
		drawStringBackgrounds(g);
		drawStrings(g);
		startNewStringLayer();
		
		{
			double mx_t = mx_index;
			double my_t = my_index;
			//double mx_t = (mouseX/(double)scalefactor);
			//double my_t = (mouseY/(double)scalefactor);
			int mi = mx_index;
			int mj = my_index;
			
			if (mi < 0)
				mi = 0;
			if (mi >= nx-1)
				mi = nx-2;
			if (mj < 0)
				mj = 0;
			if (mj >= ny-1)
				mj = ny-2;
			
			Material mat = materials[mi][mj];

			int vspacing = 12;
			int voffset = 1 + my;
			int hoffset = 5 + mx+15;
			
			if (opts.gui_tooltip.isSelected()) {
				if (voffset + 22*vspacing > ny*scalefactor) {
					voffset = voffset - ((voffset + 22*vspacing) - ny*scalefactor);
				}
				if (hoffset + 120 > ny*scalefactor) {
					hoffset = hoffset - ((hoffset + 120) - ny*scalefactor);
				}
			}
			
			String name = "Material: " + mat.type.name + (mat.modified? " (Modified)" : "");

			this.drawBigStringWithBackground(name, hoffset, voffset + 1*vspacing, g);
			if (opts.gui_tooltip.isSelected()) {
				voffset = voffset+3;
				int line = 2;
				drawTwoColumnString("E" , 							getSI(length(bilinearinterp(Ex, mx_t-0.5,my_t), bilinearinterp(Ey, mx_t,my_t-0.5)), "V/m"), hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("B" , 							getSI(parity*bilinearinterp(Bz, mx_t-0.5, my_t-0.5), "T"),	hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03d5" , 						getSI(bilinearinterp(phi,mx_t, my_t), "V"),					hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u2130" , 						getSI(mat.emf, "V/m"),										hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03b5/\u03b5\u2080" , 		getSI(mat.eps_r, ""),										hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03bc/\u03bc\u2080" , 		getSI(mat.mu_r, ""),											hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u2099" , 				getSI(bilinearinterp(rho_n,mx_t, my_t), "C/m^3"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u209A" , 				getSI(bilinearinterp(rho_p,mx_t, my_t), "C/m^3"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u2080",					getSI(bilinearinterp(rho_back,mx_t, my_t), "C/m^3"),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1" , 						getSI(bilinearinterp(rho_free,mx_t, my_t), "C/m^3"),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J\u2099" ,						getSI(length(bilinearinterp(Jx_n,mx_t-0.5,my_t), bilinearinterp(Jy_n,mx_t,my_t-0.5)), "A/m^2"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J\u209A" , 					getSI(length(bilinearinterp(Jx_p,mx_t-0.5,my_t), bilinearinterp(Jy_p,mx_t,my_t-0.5)), "A/m^2"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J" , 							getSI(length(bilinearinterp(Jx_free,mx_t-0.5,my_t), bilinearinterp(Jy_free,mx_t,my_t-0.5)), "A/m^2"),	hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("F\u2099" , 					getSI(bilinearinterp(F_n,mx_t, my_t)/q_n+bilinearinterp(phi,mx_t, my_t)-W_semi/eVtoJ, "V"),	hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("F\u209a" , 					getSI(bilinearinterp(F_p,mx_t, my_t)/q_p+bilinearinterp(phi,mx_t, my_t)-W_semi/eVtoJ, "V"),	hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("F" , 							getSI(bilinearinterp(F,mx_t, my_t), "V"),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("CMF\u2099" ,					getSI(length(bilinearinterp(cmfx_n,mx_t-0.5,my_t), bilinearinterp(cmfy_n,mx_t,my_t-0.5))/q_n, "V/m"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("CMF\u209A" , 					getSI(length(bilinearinterp(cmfx_p,mx_t-0.5,my_t), bilinearinterp(cmfy_n,mx_t,my_t-0.5))/q_p, "V/m"),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("x" , 							getSI(mx_t*ds, "m"),								hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("y" , 							getSI(width-(my_t+1)*ds, "m"),								hoffset, voffset + line*vspacing, g); line++;
			}
		}

		drawStringBackgrounds(g);
		drawStrings(g);
		startNewStringLayer();
		
		int vspacing = 13;
		int voffset = 3;
		int hoffset = 5;
		int line = 1;
		drawStringWithBackground("Time: " + getSI(time, "s"), hoffset, voffset + line*vspacing, g); line++;
		if (opts.gui_paused.isSelected())
		{
			drawStringWithBackground("Paused", hoffset, voffset + line*vspacing, g); line++;
		}
		if (sign_violation) {
			drawStringWithBackground("Warning: Numerical instability detected. Please decrease timestep.", hoffset, voffset + line*vspacing, g); line++;
		}
		if (debugging) {
			long total = Runtime.getRuntime().totalMemory();
			long used  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			drawStringWithBackground("Used memory " + getSI(used, "B"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground("Total memory " + getSI(total, "B"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground(t4.name + " " + getSI(t4.time, "s"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground(t5.name + " " + getSI(t5.avgtime, "s"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground(t6.name + " " + getSI(t6.avgtime*opts.gui_simspeed_2.getValue(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground(t7.name + " " + getSI(t7.avgtime, "s"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground(t8.name + " " + getSI(t8.avgtime, "s"), hoffset, voffset + line*vspacing, g); line++;
		}

		drawStringBackgrounds(g);
		drawStrings(g);

		if (!opts.gui_brush_highlight.isSelected()) {
			int r = (int)(scalefactor*brushsize/ds);
			int brushshape = opts.gui_brush_1.getSelectedIndex();
			if (Brush.isBrushShapeImportant(brush))
				if (brushshape == 0) {
					g.setColor(new Color(50, 50, 50));
					g.drawOval(mx - r, my - r, 2*r, 2*r);
					g.setColor(new Color(200, 200, 200));
					g.drawOval(mx - r-1, my - r-1, 2*r, 2*r);
				} else {
					g.setColor(new Color(50, 50, 50));
					g.drawRect(mx - r, my - r, 2*r-2, 2*r-2);
					g.setColor(new Color(200, 200, 200));
					g.drawRect(mx - r-1, my - r-1, 2*r, 2*r);
				}
		}
		
		t5.stop();
	}

	class Text {
		String text;
		int x;
		int y;
		int minwidth;
		boolean big;
		boolean hasBackground;
		boolean hasBorder;
		
		public Text(String text, int x, int y, boolean big, boolean hasBackground, boolean hasBorder) {
			this.text = text;
			this.x = x;
			this.y = y;
			this.big = big;
			this.hasBackground = hasBackground;
			this.hasBorder = hasBorder;
			minwidth = 0;
		}
	}
	
	public void clearStrings() {
		texts.clear();
	}
	
	public void drawStringBackgrounds(Graphics g) {
		if (!opts.gui_text_bg.isSelected())
			return;

		for (Text text : texts) {
			if (text.hasBackground) {
				if (text.big)
					g.setFont(bigfont);
				else
					g.setFont(regularfont);

				int width = Math.max(text.minwidth, g.getFontMetrics().stringWidth(text.text)+8);
				int height = g.getFontMetrics().getHeight()+4;
				int x = text.x-3;
				int y = text.y-height+6;

				g.setColor(Color.GRAY);
				g.fillRect(x-2, y-2, width+4, height+4);
			}
		}
		
		for (Text text : texts) {
			if (text.hasBackground) {
				if (text.big)
					g.setFont(bigfont);
				else
					g.setFont(regularfont);

				int width = Math.max(text.minwidth, g.getFontMetrics().stringWidth(text.text)+8);
				int height = g.getFontMetrics().getHeight()+4;
				int x = text.x-3;
				int y = text.y-height+6;

				g.setColor(Color.BLACK);
				g.fillRect(x, y, width, height);
			}
		}
	}

	public void drawStrings(Graphics g) {
		for (Text text : texts) {
			if (text.big)
				g.setFont(bigfont);
			else
				g.setFont(regularfont);
			
			g.setColor(Color.DARK_GRAY);
			g.drawString(text.text, text.x+1, text.y+1);
			g.setColor(Color.WHITE);
			g.drawString(text.text, text.x, text.y);
		}
	}
	
	public void startNewStringLayer() {
		texts.clear();
	}

	public void drawString(String str1, int x, int y, Graphics g) {
		texts.add(new Text(str1, x, y, false, false, false));
	}
	
	public void drawStringWithBackground(String str1, int x, int y, Graphics g) {
		texts.add(new Text(str1, x, y, false, true, false));
	}
	
	public void drawStringWithBackgroundAndBorder(String str1, int x, int y, Graphics g) {
		texts.add(new Text(str1, x, y, false, true, true));
	}
	
	public void drawBigString(String str1, int x, int y, Graphics g) {
		texts.add(new Text(str1, x, y, true, false, false));
	}
	
	public void drawBigStringWithBackground(String str1, int x, int y, Graphics g) {
		texts.add(new Text(str1, x, y, true, true, false));
	}
	
	public void drawTwoColumnString(String str1, String str2, int x, int y, Graphics g) {
		texts.add(new Text(String.format("%-10s", str1), x, y, false, true, true));
		texts.add(new Text(str2, x+40, y, false, true, true));
		texts.get(texts.size()-1).minwidth = 80;
	}
	
	public double length(double x, double y) {
		return Math.sqrt(x*x+y*y);
	}

	public static String getSI(double quantity, String unit) {
		if (!Double.isFinite(quantity))
			return Double.toString(quantity) + " " + unit;
		
		double mag = Math.abs(quantity);
		String precision = "%.2f";
		if (mag < 1E-18)
			return "0 " + unit;
		else if (mag < 1E-15)
			return String.format("%.2f", quantity*1e15) + " f" + unit;
		else if (mag < 1E-12)
			return String.format(precision, quantity*1e15) + " f" + unit;
		else if (mag < 1E-9)
			return String.format(precision, quantity*1e12) + " p" + unit;
		else if (mag < 1E-6)
			return String.format(precision, quantity*1e9) + " n" + unit;
		else if (mag < 1E-3)
			return String.format(precision, quantity*1e6) + " \u00b5" + unit;
		else if (mag < 1)
			return String.format(precision, quantity*1e3) + " m" + unit;
		else if (mag < 1E3)
			return String.format(precision, quantity) + " " + unit;
		else if (mag < 1E6)
			return String.format(precision, quantity*1e-3) + " k" + unit;
		else if (mag < 1E9)
			return String.format(precision, quantity*1e-6) + " M" + unit;
		else if (mag < 1E12)
			return String.format(precision, quantity*1e-9) + " G" + unit;
		else
			return String.format(precision, quantity*1e-12) + " T" + unit;
	}
	
	public boolean readFile()
	{
		
		File testfile = new File(startingpath);
		if (!testfile.canRead()) {
			JOptionPane.showMessageDialog(opts,
				    "Error: Java does not have access to this folder. Please see instructions to fix this issue.");
		}
		
		try {
			EventQueue.invokeAndWait(new Runnable() {
			    @Override
			    public void run() {
			    	JFileChooser fd = new JFileChooser(startingpath);
					fd.setFileFilter(new FileFilter(){
						public boolean accept(File f) {
							if (f.isDirectory()) {
								return true;
							}
							if (f.getName().endsWith(fileextension)) return true;
							return false;
						}
						@Override
						public String getDescription() {
							return fileextension;
						}
					});
					fd.setVisible(true);
					int result = fd.showOpenDialog(opts);
					startingpath = fd.getCurrentDirectory().getPath();
					
					if (result == JFileChooser.APPROVE_OPTION)
						infile = fd.getSelectedFile();
					else
						infile = null;
			    }
			});
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		if (infile == null || !infile.exists()) return false;
		try {
			JsonReader fstr = new JsonReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile))));
			Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
			
			fstr.beginObject();
			if (!fstr.nextName().equals("version")) {
				fstr.close();
				throw new IllegalArgumentException();
			}
			
			int version = fstr.nextInt();
			if (version != 1) {
				fstr.close();
				throw new IllegalArgumentException();
			}

			JOptionPane optionPane = new JOptionPane("Loading file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
			JDialog dialog = optionPane.createDialog("Loading");

			dialog.setModal(false);
			dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
			dialog.setVisible(true);


			initializeGrid(width, default_resolution);

			while (fstr.hasNext()) {
                String name = fstr.nextName();
                switch (name){
                	case "time": time = fstr.nextDouble(); break;
                	case "gui_paused": opts.gui_paused.setSelected(fstr.nextBoolean()); break;
                	case "gui_tooltip": opts.gui_tooltip.setSelected(fstr.nextBoolean()); break;
                	case "gui_text_bg": opts.gui_text_bg.setSelected(fstr.nextBoolean()); break;
                	case "gui_view": opts.gui_view.setSelectedIndex(fstr.nextInt()); break;
                	case "gui_view_vec": opts.gui_view_vec.setSelectedIndex(fstr.nextInt()); break;
                	case "gui_view_vec_mode": opts.gui_view_vec_mode.setSelectedIndex(fstr.nextInt()); break;
                	case "gui_simspeed": opts.gui_simspeed.setValue(fstr.nextInt()); break;
                	case "gui_simspeed_2": opts.gui_simspeed_2.setValue(fstr.nextInt()); break;
                	case "gui_brightness": opts.gui_brightness.setValue(fstr.nextInt()); break;
                	case "gui_brightness_vec": opts.gui_brightness_vec.setValue(fstr.nextInt()); break;
                	case "gui_elem_colors": opts.gui_elem_colors.setSelected(fstr.nextBoolean()); break;
                	case "gui_bc": opts.gui_bc.setSelectedIndex(fstr.nextInt()); break;
                	case "description": opts.textPane.setText(fstr.nextString()); break;

                	case "ex": Ex = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "ey": Ey = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "hz": Hz = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

                	case "rho_c": rho_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "rho_n": rho_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "rho_p": rho_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "rho_back": rho_back = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "rho_free": rho_free = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

                	case "jx_c": Jx_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "jy_c": Jy_abs = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "jx_n": Jx_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "jy_n": Jy_n = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "jx_p": Jx_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;
                	case "jy_p": Jy_p = validateArraySize((double[][]) gson.fromJson(fstr, double[][].class)); break;

                	case "materials": materials = validateArraySize((Material[][]) gson.fromJson(fstr, Material[][].class)); break;

                	case "voltageprobes": voltageprobes = new CopyOnWriteArrayList<>(Arrays.asList(
                			(VoltageProbe[]) gson.fromJson(fstr, VoltageProbe[].class))); break;
                	case "currentprobes": currentprobes = new CopyOnWriteArrayList<>(Arrays.asList(
                			(CurrentProbe[]) gson.fromJson(fstr, CurrentProbe[].class))); break;

                	case "ground": ground = (VoltageProbe) gson.fromJson(fstr, VoltageProbe.class); break;

                    default: fstr.skipValue(); break; // skip others
                }
            }
			fstr.endObject();
			fstr.close();
			
			opts.textPane.setEditable(false);
			opts.textPane.setCaretPosition(0);
			constructBoundary();
			initializeAllMaterials();
			updateAllMaterials();
			calcMiscFields(true);

			dialog.dispose();
			opts.setTitle("Brandon's semiconductor simulator - " + infile.getName());
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException | IllegalArgumentException e) {
			JOptionPane.showMessageDialog(opts,
				    "Error: Unable to load file.");
			e.printStackTrace();
			return false;
		}
		return true;
	}

	public Material[][] validateArraySize(Material[][] array) throws RuntimeException {
		Material[][] field = new Material[nx][ny];
		
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				if (i < nx && j < ny) {
					field[i][j] = array[i][j];
				}
			}
		}
		
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (field[i][j] == null || field[i][j].type == MaterialType.ABSORBER)
					field[i][j] = new Material();
			}
		}
		
		return field;
	}
	
	public double[][] validateArraySize(double[][] array) throws RuntimeException {
		double[][] field = new double[nx][ny];
		
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				if (i < nx && j < ny) {
					field[i][j] = array[i][j];
				}
			}
		}
		
		return field;
	}
	
	public boolean writeFile()
	{
		File testfile = new File(startingpath);
		if (!testfile.canWrite()) {
			JOptionPane.showMessageDialog(opts,
				    "Error: Java does not have access to this folder. Please see instructions to fix this issue.");
		}
		
		try {
			EventQueue.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					JFileChooser fd = new JFileChooser(startingpath);
					fd.setFileFilter(new FileFilter(){
						public boolean accept(File f) {
							if (f.isDirectory()) {
								return true;
							}
							if (f.getName().endsWith(fileextension)) return true;
							return false;
						}
						@Override
						public String getDescription() {
							return fileextension;
						}
					});
					int result = fd.showSaveDialog(opts);
					startingpath = fd.getCurrentDirectory().getPath();
					
					if (result == JFileChooser.APPROVE_OPTION)
						outfile = fd.getSelectedFile();
					else
						outfile = null;
				}
			});
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		if (outfile == null) return false;
		if (!outfile.getName().endsWith(fileextension))
			outfile = new File(outfile.getAbsolutePath() + fileextension);
		
		if (outfile.exists()) {
			int result = JOptionPane.showConfirmDialog(opts, "A file with that name already exists. Do you wish to overwrite it?", "Save file", JOptionPane.YES_NO_OPTION);
			if (result != JOptionPane.OK_OPTION)
				return false;
		}
		
		try {
			PrintWriter fstr = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outfile)));
			
			Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
			
	        JsonObject obj = new JsonObject();
	        


			JOptionPane optionPane = new JOptionPane("Saving file, please wait.", JOptionPane.INFORMATION_MESSAGE, JOptionPane.DEFAULT_OPTION, null, new Object[]{}, null);
			JDialog dialog = optionPane.createDialog("Saving");

			dialog.setModal(false);
			dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
			dialog.setVisible(true);
	        
	        // Version should always be first
			obj.addProperty("version", saveversion);
			obj.addProperty("time", time);
			
			obj.addProperty("gui_paused", opts.gui_paused.isSelected());
			obj.addProperty("gui_tooltip", opts.gui_tooltip.isSelected());
			obj.addProperty("gui_text_bg", opts.gui_text_bg.isSelected());
			obj.addProperty("gui_view", opts.gui_view.getSelectedIndex());
			obj.addProperty("gui_view_vec", opts.gui_view_vec.getSelectedIndex());
			obj.addProperty("gui_view_vec_mode", opts.gui_view_vec_mode.getSelectedIndex());
			obj.addProperty("gui_simspeed", opts.gui_simspeed.getValue());
			obj.addProperty("gui_simspeed_2", opts.gui_simspeed_2.getValue());
			obj.addProperty("gui_brightness", opts.gui_brightness.getValue());
			obj.addProperty("gui_brightness_vec", opts.gui_brightness_vec.getValue());
			obj.addProperty("gui_elem_colors", opts.gui_elem_colors.isSelected());
			obj.addProperty("gui_bc", opts.gui_bc.getSelectedIndex());
			obj.addProperty("description", opts.textPane.getText());
			
			obj.add("ex", gson.toJsonTree(Ex));
			obj.add("ey", gson.toJsonTree(Ey));
			obj.add("hz", gson.toJsonTree(Hz));

			obj.add("rho_c", gson.toJsonTree(rho_abs));
			obj.add("rho_n", gson.toJsonTree(rho_n));
			obj.add("rho_p", gson.toJsonTree(rho_p));
			obj.add("rho_back", gson.toJsonTree(rho_back));
			obj.add("rho_free", gson.toJsonTree(rho_free));
			
			obj.add("jx_c", gson.toJsonTree(Jx_abs));
			obj.add("jy_c", gson.toJsonTree(Jy_abs));
			obj.add("jx_n", gson.toJsonTree(Jx_n));
			obj.add("jy_n", gson.toJsonTree(Jy_n));
			obj.add("jx_p", gson.toJsonTree(Jx_p));
			obj.add("jy_p", gson.toJsonTree(Jy_p));
			
			obj.add("materials", gson.toJsonTree(materials));

			obj.add("voltageprobes", gson.toJsonTree((VoltageProbe[]) voltageprobes.toArray(new VoltageProbe[voltageprobes.size()])));
			obj.add("currentprobes", gson.toJsonTree((CurrentProbe[]) currentprobes.toArray(new CurrentProbe[currentprobes.size()])));
			obj.add("ground", gson.toJsonTree(ground));
			
			String json = gson.toJson(obj);
			fstr.print(json);
			fstr.flush();
			fstr.close();
			
			dialog.dispose();
			opts.setTitle("Brandon's semiconductor simulator - " + outfile.getName());
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}

	@Override
	public void mouseClicked(MouseEvent arg0) {}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent e) {
		mouse_pressed = true;
		mousebutton = e.getButton();
		mx = e.getX();
		my = e.getY();
		mx_start = e.getX();
		my_start = e.getY();
	}
	@Override
	public void mouseReleased(MouseEvent e) {
		mouse_pressed = false;
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		mx = e.getX();
		my = e.getY();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		mx = arg0.getX();
		my = arg0.getY();
	}

	@SuppressWarnings("serial")
	private Action key_pause = new AbstractAction(null) {
		@Override
		public void actionPerformed(ActionEvent e) {
			if (texting) return;
			opts.gui_paused.setSelected(!opts.gui_paused.isSelected());
		}
	};

    @SuppressWarnings("serial")
    private Action key_frame = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
			advanceframe = true;
        }
    };

    @SuppressWarnings("serial")
    private Action key_dbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
			debugging = !debugging;
			Timer.allEnabled = debugging;
        }
    };

    @SuppressWarnings("serial")
    private Action key_changebrush = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		opts.gui_brush_1.setSelectedIndex((opts.gui_brush_1.getSelectedIndex()+1)%2);
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_shift = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		shift_down = true;
        }
    };
    @SuppressWarnings("serial")
    private Action key_shift_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		shift_down = false;
        }
    };
    @SuppressWarnings("serial")
    private Action key_ctrl = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		ctrl_down = true;
        }
    };
    @SuppressWarnings("serial")
    private Action key_ctrl_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		ctrl_down = false;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_cut = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		cut = true;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_copy = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		copy = true;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_paste = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		paste = true;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_delete = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		delete = true;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_color = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		opts.gui_elem_colors.setSelected(!opts.gui_elem_colors.isSelected());
        }
    };
    
    ScalarView prev_scalar_view = ScalarView.NONE;
    VectorView prev_vector_view = VectorView.NONE;
    
    @SuppressWarnings("serial")
    private Action key_scalar_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		if (opts.gui_view.getSelectedItem() == ScalarView.NONE)
    			opts.gui_view.setSelectedItem(prev_scalar_view);
    		else
    		{
    			prev_scalar_view = (ScalarView) opts.gui_view.getSelectedItem();
    			opts.gui_view.setSelectedItem(ScalarView.NONE);
    		}
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_vector_view = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		if (opts.gui_view_vec.getSelectedItem() == VectorView.NONE)
    			opts.gui_view_vec.setSelectedItem(prev_vector_view);
    		else
    		{
    			prev_vector_view = (VectorView) opts.gui_view_vec.getSelectedItem();
    			opts.gui_view_vec.setSelectedItem(VectorView.NONE);
    		}
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_tooltip = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		opts.gui_tooltip.setSelected(!opts.gui_tooltip.isSelected());
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_textbg = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			if (texting) return;
    		opts.gui_text_bg.setSelected(!opts.gui_text_bg.isSelected());
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_alt = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
    		alt_down = true;
        }
    };
    
    @SuppressWarnings("serial")
    private Action key_alt_up = new AbstractAction(null) {
		@Override
        public void actionPerformed(ActionEvent e) {
			alt_down = false;
        }
    };

    public void addKeyBinds(JPanel contentPane) {
    	InputMap map = contentPane.getInputMap(JComponent.WHEN_FOCUSED);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_P, 0), key_pause);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0), key_pause);
    	contentPane.getActionMap().put(key_pause, key_pause);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_F, 0), key_frame);
    	contentPane.getActionMap().put(key_frame, key_frame);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, 0), key_dbg);
    	contentPane.getActionMap().put(key_dbg, key_dbg);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_Q, 0), key_changebrush);
    	contentPane.getActionMap().put(key_changebrush, key_changebrush);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, KeyEvent.SHIFT_DOWN_MASK), key_shift);
    	contentPane.getActionMap().put(key_shift, key_shift);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_SHIFT, 0, true), key_shift_up);
    	contentPane.getActionMap().put(key_shift_up, key_shift_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, KeyEvent.CTRL_DOWN_MASK), key_ctrl);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, KeyEvent.META_DOWN_MASK), key_ctrl);
    	contentPane.getActionMap().put(key_ctrl, key_ctrl);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_CONTROL, 0, true), key_ctrl_up);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_META, 0, true), key_ctrl_up);
    	contentPane.getActionMap().put(key_ctrl_up, key_ctrl_up);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_X, KeyEvent.CTRL_DOWN_MASK), key_cut);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_X, KeyEvent.META_DOWN_MASK), key_cut);
    	contentPane.getActionMap().put(key_cut, key_cut);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, KeyEvent.CTRL_DOWN_MASK), key_copy);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, KeyEvent.META_DOWN_MASK), key_copy);
    	contentPane.getActionMap().put(key_copy, key_copy);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, KeyEvent.CTRL_DOWN_MASK), key_paste);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, KeyEvent.META_DOWN_MASK), key_paste);
    	contentPane.getActionMap().put(key_paste, key_paste);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE, 0), key_delete);
    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0), key_delete);
    	contentPane.getActionMap().put(key_delete, key_delete);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_C, 0), key_color);
    	contentPane.getActionMap().put(key_color, key_color);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, 0), key_scalar_view);
    	contentPane.getActionMap().put(key_scalar_view, key_scalar_view);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_V, 0), key_vector_view);
    	contentPane.getActionMap().put(key_vector_view, key_vector_view);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_T, 0), key_tooltip);
    	contentPane.getActionMap().put(key_tooltip, key_tooltip);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_G, 0), key_textbg);
    	contentPane.getActionMap().put(key_textbg, key_textbg);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, KeyEvent.ALT_DOWN_MASK), key_alt);
    	contentPane.getActionMap().put(key_alt, key_alt);

    	map.put(KeyStroke.getKeyStroke(KeyEvent.VK_ALT, 0, true), key_alt_up);
    	contentPane.getActionMap().put(key_alt_up, key_alt_up);
    	
    }
    
	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		opts.gui_brushsize.setValue(opts.gui_brushsize.getValue() - (int)(10*e.getPreciseWheelRotation()));
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == opts.gui_reset)
			clear = true;
		else if (e.getSource() == opts.gui_resetall)
			reset = true;
		else if (e.getSource() == opts.gui_save)
			save = true;
		else if (e.getSource() == opts.gui_open)
			load = true;
		else if (e.getSource() == opts.gui_help)
			help.setVisible(true);
		else if (e.getSource() == opts.gui_editdesc) {
			opts.textPane.setEditable(!opts.textPane.isEditable());
		} else if (e.getSource() == opts.gui_view) {
			updateMiscFields = true;
		} else if (e.getSource() == opts.gui_view_vec) {
			updateMiscFields = true;
		} else if (e.getSource() == opts.gui_brush) {
			brush_changed = true;
		}
	}
	
	
	enum Brush {
		INTERACT("Interact"),
		DRAW("Draw"),
		VOLTAGE("Add voltage probe"),
		CURRENT("Add current probe"),
		GROUND("Add ground"),
		DELETEPROBE("Delete probe"),
		REPLACE("Replace"),
		LINE("Line"),
		FILL("Fill"),
		ERASE("Eraser"),
		SELECT("Select & Move"),
		FLOODSELECT("Select region"),
		TEXT("Add text");
		
		String name;
		Brush(String name)
		{
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
		
		public static boolean isMaterialModifyingBrush(Brush brush) {
			return (brush == Brush.DRAW
					|| brush == Brush.LINE
					|| brush == Brush.REPLACE
					|| brush == Brush.ERASE
					|| brush == Brush.FILL);
		}
		
		public static boolean isBrushShapeImportant(Brush brush) {
			return (brush == Brush.DRAW
					|| brush == Brush.LINE
					|| brush == Brush.REPLACE
					|| brush == Brush.ERASE);
		}
	};

	enum BrushShape {
		CIRCLE("Circle brush"),
		SQUARE("Square brush");
		
		String name;
		BrushShape(String name)
		{
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	};

	enum ScalarView {
		NONE("No scalar overlay"),
		E_FIELD("View E field magnitude"),
		B_FIELD("View B field"),
		CHARGE("View \u03c1: Net charge density"),
		CURRENT("View J: Total current magnitude"),
		H_FIELD("View H field"),
		POTENTIAL("View \u03d5: Electric scalar potential"),
		ENERGY("View u: Electromagnetic energy density"),
		ELECTRON_CHARGE("View \u03c1\u2099: Electron charge density"),
		HOLE_CHARGE("View \u03c1\u209A: Hole charge density"),
		COMBINED_CHARGE("View: Combined electron+hole charge density"),
		BACKGROUND_CHARGE("View \u03c1\u2080: Background charge density"),
		HEAT("View Q: Heat dissipation"),
		ENTROPY("View s: Entropy generation (Free energy dissipation)"),
		ELECTRON_POTENTIAL("View F\u2099: Electron chemical potential"),
		HOLE_POTENTIAL("View F\u209A: Hole chemical potential"),
		AVERAGE_POTENTIAL("View F: Average electrochemical potential"),
		RECOMBINATION("View R: Recombination rate"),
		LIGHT("View: Emitted light"),
		DEBUG("Debug");
		
		String name;
		ScalarView(String name)
		{
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	};

	enum VectorView {
		NONE("No vector overlay"),
		E_FIELD("View E field"),
		D_FIELD("View D field"),
		ELECTRON_CURRENT("View J\u2099: Electron current"),
		HOLE_CURRENT("View J\u209A: Hole current"),
		TOTAL_CURRENT("View J: Total current"),
		EMF("View \u2130: External electromotive force"),
		POYNTING("View S: Poynting vector");
		
		String name;
		VectorView(String name)
		{
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	};
	
	enum VectorMode {
		ARROWS("Show vectors"),
		LINES("Show lines");
		
		String name;
		VectorMode(String name)
		{
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	};
	
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

	@Override
	public void keyTyped(KeyEvent e) {
		if (texting) {
			if (Font7x5.getCharacter(e.getKeyChar()) != null) {
				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 7; j++) {
						if (text_x+i+1 >= 0 && text_x+i+1 < nx && text_y+j >= 0 && text_y+j < ny
								&& Font7x5.getPixel(e.getKeyChar(), 4-i, j) == 1 && materials[text_x+i+1][text_y+j].type == MaterialType.VACUUM) {
							initializeMaterial(text_x+i+1, text_y+j, MaterialType.DECO);
						}
					}
				}
				text_x += 6;
			}
			updateAllMaterials();
		}
	}

	@Override
	public void keyPressed(KeyEvent e) {
		if (texting) {
			if (e.getKeyCode() == KeyEvent.VK_BACK_SPACE || e.getKeyCode() == KeyEvent.VK_DELETE) {
				text_x -= 6;
				for (int i = 0; i < 5; i++) {
					for (int j = 0; j < 7; j++) {
						if (materials[text_x+i+1][text_y+j].type == MaterialType.DECO) {
							materials[text_x+i+1][text_y+j].erase();
						}
					}
				}
			}
			updateAllMaterials();
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {}
}

class CurrentProbe {
	int x1;
	int y1;
	int x2;
	int y2;
	
	double current = 0;
}

class VoltageProbe {
	int x;
	int y;
	
	double potential = 0;
}

class RenderCanvas extends JPanel {
	
	private static final long serialVersionUID = 7369516276529576171L;
	
	Electrodynamics parent;
	@Override
	public void paintComponent(Graphics real) {
		parent.render();
		real.drawImage(parent.screen, 0, 0, parent.opts);
	}
	
	public RenderCanvas(Electrodynamics w) {
		parent = w;
	}
}

class Timer {
	long tstart = 0;
	String name;
	boolean enabled = true;
	boolean outputavg = false;
	double avgtime = 0;
	double time = 0;
	static boolean allEnabled = false;
	
	public Timer(String name, boolean enabled) {
		this.name = name;
		this.enabled = enabled;
	}

	void start() {
		if (enabled) {
			tstart = System.nanoTime();
		}
	}

	void disableOutput() {
		enabled = false;
	}
	void enableOutput() {
		enabled = true;
	}

	void stop() {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.99+time*0.01;
		}
	}

	void stop(String msg) {
		if (allEnabled && enabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			time = diff/1e9;
			avgtime = avgtime*0.95+time*0.05;
		}
	}
}

class Vector {
	double x;
	double y;
	
	public Vector(double x, double y) {
		this.x = x;
		this.y = y;
	}
	
	public Vector copy() {
		return new Vector(x, y);
	}
	
	public void copy(Vector b) {
		this.x = b.x;
		this.y = b.y;
	}
	
	public void initialize(double x, double y) {
		this.x = x;
		this.y = y;
	}
	
	public void add(Vector b) {
		x += b.x;
		y += b.y;
	}
	
	public void scalarmult(double c) {
		x *= c;
		y *= c;
	}

	public void addmult(Vector b, double c) {
		x += b.x * c;
		y += b.y * c;
	}
	
	public void rotate(double theta) {
		double xf = x*Math.cos(theta) + y*Math.sin(theta);
		double yf = -x*Math.sin(theta) + y*Math.cos(theta);
		x = xf;
		y = yf;
	}
	
	public void normalize() {
		double magnitude = Math.sqrt(x*x+y*y);
		if (magnitude != 0) {
			x /= magnitude;
			y /= magnitude;
		}
	}

	public double dot(Vector b) {
		return this.x * b.x + this.y * b.y;
	}
}

enum MaterialType
{
	
	EMF					("Voltage source (Adjustable)",			230, 216, 46, 230),
	SWITCH				("Switch",								194, 194, 194, 120),
	METAL				("Metal",								153, 153, 153, 120),
	METAL_HIGH_C		("Conductive metal",					191, 191, 191, 120),
	METAL_LOW_C			("Resistive metal",						94, 94, 94, 120),
	METAL_HIGH_W		("High workfunction metal",				163, 116, 116, 120),
	METAL_LOW_W			("Low workfunction metal",				116, 121, 163, 120),
	SEMI				("Intrinsic semiconductor",				207,  161, 212, 120),
	SEMI_P_TYPE			("P-type semiconductor",				191,  74,  34, 120),
	SEMI_N_TYPE			("N-type semiconductor",				 84, 123, 191, 120),
	SEMI_HEAVY_P_TYPE	("Heavily doped P-type semiconductor",	204,  41,  41, 120),
	SEMI_HEAVY_N_TYPE	("Heavily doped N-type semiconductor",	 39,  52, 194, 120),
	SEMI_LIGHT_P_TYPE	("Lightly doped P-type semiconductor",	201, 131,  73, 120),
	SEMI_LIGHT_N_TYPE	("Lightly doped N-type semiconductor",	137, 188, 204, 120),
	DIELECTRIC			("Dielectric",							 81, 171,  51, 120),
	FERROMAGNET			("Ferromagnet",							116, 50, 117, 120),
	POS_CHARGE			("Positive static charge",				116, 50, 50, 120),
	NEG_CHARGE			("Negative static charge",				50, 50, 117, 120),
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