package electrodynamics;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.TimerTask;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;

public class Electrodynamics extends TimerTask implements MouseListener, MouseMotionListener, MouseWheelListener, ActionListener {
	
	//different aspect ratios, multithreading, change size
	/*
	 * Demonstrations: Skin effect, LRC circuit, Simple circuit/speed of light demo,
	 * AC circuit (effect of cap/inductance), inductor, capacitor,
	 * mutual inductance/transformer, dielectric,
	 * charge escaping from center of conductor,
	 * optical fiber, antenna communication, faraday cage, ohm's law,
	 * rc circuit, lr circuit, resonant frequency, circuit maze, conductors at equilibrium
	 * Wave diffracting around obstacle
	 */
	
	/* Physical fields */
	double[][] ex; //x component of E field
	double[][] ey; //y component of E field
	double[][] bz; //z component of B field
	double[][] rho; //charge density
	double[][] jx; //x component of current density
	double[][] jy; //y component of current density

	double[][] bzp; //Temporary storage for B field

	/* Used in poisson solver */
	double[][] deltarho;
	double[][] deltarhoMG;
	double[][][] epsxMG;
	double[][][] epsyMG;
	double[][] epsavgMG;
	double[][] poisson1;
	double[][] poisson2;

	/* Material properties */
	Material[][] materials;
	double[][] diffx; //Charge carrier diffusivity
	double[][] diffy;
	double[][] sigmax; //Conductivity
	double[][] sigmay;
	double[][] emfx; //EMF
	double[][] emfy;
	double[][] epsrx; //Dielectric x
	double[][] epsry; //Dielectric y
	double[][] mu_rz; //Relative permeability

	int maxfreq = 21;
	int[][] freqx;
	int[][] freqy;
	double[] amplitudes;
	
	boolean[][] visited; //Used in flood fill

	/* Auxiliary fields */
	double[][] dx; //Displacement x
	double[][] dy; //Displacement y
	double[][] hz; //H field
	double[][] sx; //Poynting x
	double[][] sy; //Poynting y
	double[][] u; //Energy density
	
	/* Domain size parameters */
	double width = 50;
	double height = width;
	double ds = 0.2;
	double dt = 0;
	int nx = (int)(width/ds)+1;
	int ny = (int)(height/ds)+1;
	
	/* Physical and numerical constants */
	//double eps0 = 8.85e-12;
	//double mu0 = 1.257e-6;
	double eps0 = 1/30.0;
	double mu0 = 1/30.0;
	double c = 1/Math.sqrt(eps0*mu0);
	double csq = 1/(eps0*mu0);
	double tauminus = 1-c*dt/ds;
	double tauplus = 1+c*dt/ds;
	double maxdt = 0.95*ds/(Math.sqrt(2)*c);

	double D_max = 2.0; //Diffusivity of charge carriers in material
	double sigma_max = 10.0; //Conductivity of material
	
	double time = 0.0;
	int stepnumber = 0;
	int lastsimspeed = 0;
	double phase = 0;
	
	boolean advanceframe = false;
	boolean clear = false;
	boolean reset = false;
	boolean save = false;
	boolean load = false;

	/* Graphics */
	RenderCanvas r;
	MainWindow opts;
	HelpDialog help;
	BufferedImage screen;
	int scalefactor = 3;
	int imgwidth = 0;
	int imgheight = 0;
	int targetframerate = 60;
	int frameduration = 1000/targetframerate;
	
	/* Timing */
	Timer t1 = new Timer("Poisson 1", false);
	Timer t2 = new Timer("Poisson 2", false);
	Timer t3 = new Timer("Poisson 3 intermediate", false);
	Timer t4 = new Timer("Poisson 3 total", false);
	Timer t5 = new Timer("Graphics", false);
	Timer t6 = new Timer("FD loop", false);
	
	static int BRUSH_INTERACT = 0;
	static int BRUSH_CHARGE = 1;
	static int BRUSH_CONDUCTOR = 2;
	static int BRUSH_EMF = 3;
	static int BRUSH_DIELECTRIC = 4;
	static int BRUSH_PARAMAGNET = 5;
	static int BRUSH_SWITCH = 6;
	
	static int BRUSH_SHAPE_CIRCLE = 0;
	static int BRUSH_SHAPE_SQUARE = 1;
	
	static int VIEWSCALAR_NONE = 0;
	static int VIEWSCALAR_E = 1;
	static int VIEWSCALAR_B = 2;
	static int VIEWSCALAR_RHO = 3;
	static int VIEWSCALAR_J = 4;
	static int VIEWSCALAR_H = 5;
	static int VIEWSCALAR_U = 6;
	
	static int VIEWVEC_NONE = 0;
	static int VIEWVEC_E = 1;
	static int VIEWVEC_J = 2;
	static int VIEWVEC_D = 3;
	static int VIEWVEC_S = 4;

	/* Mouse controls */
	boolean drawingcharge = false;
	boolean linemode = false;
	int lineorientation = 0;
	int mouseX = 0;
	int mouseY = 250;
	int mouseXstart = 0;
	int mouseYstart = 250;
	int mouseButton = 0;
	boolean mouseIsPressed = false;
	double mxp = 0.0;
	double myp = 0.0;

	String fileextension = ".txt";
	Font bigfont = new Font("", Font.PLAIN, 15);
	Font regularfont = new Font("", Font.PLAIN, 12);

	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(
					UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}

		Electrodynamics w = new Electrodynamics();
		java.util.Timer t = new java.util.Timer();
		t.schedule(w, 0, w.frameduration);
	}

	public Electrodynamics() {
		checkCFL();
		
		opts = new MainWindow();

		ex = new double[nx][ny];
		ey = new double[nx][ny];
		bz = new double[nx][ny];
		rho = new double[nx][ny];
		jx = new double[nx][ny];
		jy = new double[nx][ny];

		bzp = new double[nx][ny];

		deltarho = new double[nx][ny];
		deltarhoMG = new double[nx][ny];
		poisson1 = new double[nx][ny];
		poisson2 = new double[nx][ny];
		
		materials = new Material[nx][ny];
		diffx = new double[nx][ny];
		diffy = new double[nx][ny];
		sigmax = new double[nx][ny];
		sigmay = new double[nx][ny];
		emfx = new double[nx][ny];
		emfy = new double[nx][ny];
		epsrx = new double[nx][ny];
		epsry = new double[nx][ny];
		freqx = new int[nx][ny];
		freqy = new int[nx][ny];
		mu_rz = new double[nx][ny];
		amplitudes = new double[maxfreq];
		visited = new boolean[nx][ny];
		epsxMG = new double[11][nx][ny];
		epsyMG = new double[11][nx][ny];
		epsavgMG = new double[nx][ny];

		dx = new double[nx][ny];
		dy = new double[nx][ny];
		hz = new double[nx][ny];
		sx = new double[nx][ny];
		sy = new double[nx][ny];
		u = new double[nx][ny];

		resetfields(true);
		prescaleDielectric();
		correctEfield();
		
		r= new RenderCanvas(this);

		opts.setVisible(true);
		imgwidth = (int)Math.ceil(scalefactor*nx);
		imgheight = (int)Math.ceil(scalefactor*ny);
		screen = (BufferedImage) opts.createImage(imgwidth, imgheight);
		
		opts.add(r, BorderLayout.CENTER);
		r.addMouseListener(this);
		r.addMouseMotionListener(this);
		r.addMouseWheelListener(this);
		r.setPreferredSize(new Dimension(imgwidth-1, imgheight));
		opts.gui_reset.addActionListener(this);
		opts.gui_resetall.addActionListener(this);
		opts.gui_save.addActionListener(this);
		opts.gui_open.addActionListener(this);
		opts.gui_help.addActionListener(this);
		opts.gui_editdesc.addActionListener(this);
		
		opts.pack();
		
		InputMap im = (InputMap)UIManager.get("Button.focusInputMap");
		im.put(KeyStroke.getKeyStroke("pressed SPACE"), "none");
		im.put(KeyStroke.getKeyStroke("released SPACE"), "none");
		
		for (int i = 0; i < keys.length; i++) {
			Action key = keys[i];
			opts.contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
					keystrokes[i], key);
			opts.contentPane.getActionMap().put(key, key);
		}
		
		help = new HelpDialog();
		help.setVisible(false);
	}

	public void run() {
		int iterationmultiplier = opts.gui_simspeed_2.getValue();

		if (clear || reset) {
			resetfields(reset);
			prescaleDielectric();
			correctEfield();
			time = 0.0;
			clear = false;
			reset = false;
		}

		if (opts.gui_simspeed.getValue() != lastsimspeed) {
			lastsimspeed = opts.gui_simspeed.getValue();
			dt = maxdt*(lastsimspeed/20.0);
			tauminus = 1-c*dt/ds;
			tauplus = 1+c*dt/ds;
		}

		handleMouseInput();

		if (!opts.gui_paused.isSelected() || advanceframe) {
			for (int i = 0; i < iterationmultiplier ; i++) {
				this.updateMainFields();
			}
		}

		calcAuxillaryFields();

		if (save) {
			writeFile();
			save = false;
		}

		if (load) {
			readFile();
			load = false;
		}

		r.repaint();
	}
	
	public void checkCFL() {
		if (c > ds/(Math.sqrt(2)*dt)) {
			System.out.println("Error: Wave speed too high. Ratio = " + (c/(ds/(Math.sqrt(2*dt)))));
			System.exit(0);
		}
		if (dt > ds*ds/(4*D_max)) {
			System.out.println("Error: Diffusion coefficient too high. Ratio = " + dt/(ds*ds/(4*D_max)));
			System.exit(0);
		}
	}

	public double sigmoid(double x) {
		return 1.0/(1.0+Math.exp(-4.0*x));
	}
	
	public void resetfields(boolean resetall) {
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				ex [i][j] = 0.0;
				ey [i][j] = 0.0;
				bz [i][j] = 0.0;
				rho [i][j] = 0.0;
				jx [i][j] = 0.0;
				jy [i][j] = 0.0;

				bzp [i][j] = 0.0;

				poisson1 [i][j] = 0.0;
				poisson2 [i][j] = 0.0;
				
				dx[i][j] = 0.0;
				dy[i][j] = 0.0;
				hz[i][j] = 0.0;
				sx[i][j] = 0.0;
				sy[i][j] = 0.0;
				u[i][j] = 0.0;
				visited[i][j] = false;
				
				if (resetall) {
					//double x = i*ds;
					//double y = j*ds;
					//rho[i][j] = (1+Math.signum(200.0 - (x-25)*(x-25) - (y-25)*(y-25)));
					//double dist = Math.sqrt((x-25)*(x-25) + (y-25)*(y-25));
					//rho[i][j] = 100.0*Math.exp(-dist*dist);
					
					if (materials[i][j] == null)
						materials[i][j] = new Material();
					else
						materials[i][j].reset();
					
					sigmax[i][j] = 0.0;
					sigmay[i][j] = 0.0;
					emfx[i][j] = 0.0;
					emfy[i][j] = 0.0;
					epsrx[i][j] = 0.0;
					epsry[i][j] = 0.0;
					diffx[i][j] = 0.0;
					diffy[i][j] = 0.0;
					epsrx[i][j] = 1.0;
					epsry[i][j] = 1.0;
					mu_rz[i][j] = 1.0;

					//double dist2 = Math.sqrt((x-25)*(x-25) + (y-30)*(y-30));
					//epsx[i][j] = 1 + 10.0*sigmoid(-(dist2-2.0));
					//epsy[i][j] = 1 + 10.0*sigmoid(-(dist2-2.0));
				}
			}
		}
	}
	
	public void handleMouseInput() {
		int brush = opts.gui_brush.getSelectedIndex();
		int brushshape = opts.gui_brush_1.getSelectedIndex();
		
		boolean rising = false;
		
		if (mouseIsPressed) {
			if (!drawingcharge) {
				drawingcharge = true;
				rising = true;
				if (brush != BRUSH_INTERACT)
					opts.gui_paused.setSelected(true);
			}
		} else {
			if (drawingcharge) {
				drawingcharge = false;
				prescaleDielectric();
				correctEfield();
			}
		}


		double mx = ((mouseX+1)/(double)scalefactor)*ds;
		double my = ((mouseY+1)/(double)scalefactor)*ds;

		if (brush == BRUSH_INTERACT) {
			int mi = (int)Math.round(((mouseX+1)/(double)scalefactor) - 0.5);
			int mj = (int)Math.round(((mouseY+1)/(double)scalefactor) - 0.5);
			
			if (materials[mi][mj].type == Material.SWITCH || materials[mi][mj].type == Material.EMFSOURCE) {
				r.setCursor(new Cursor(Cursor.HAND_CURSOR));
			} else {
				r.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
			}
		}

		opts.gui_stepsizelbl.setText("Step size: " + getSI(dt, "s"));
		opts.gui_stepslbl.setText("Steps/frame: " + opts.gui_simspeed_2.getValue());
		
		double intensity = 1.0*Math.pow(10.0, opts.gui_brushintensity2.getValue()/10.0);
		double brushsize = 1.5*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
		
		if (brush == BRUSH_CONDUCTOR || brush == BRUSH_SWITCH || brush == BRUSH_INTERACT) {
			opts.gui_brushintensity.setText("Brush Intensity");
			opts.gui_brushintensity.setEnabled(false);
			opts.gui_brushintensity2.setEnabled(false);
		} else if (brush == BRUSH_EMF) {
			opts.gui_brushintensity.setText("EMF Strength: " + getSI(intensity, "V/m"));
			opts.gui_brushintensity.setEnabled(true);
			opts.gui_brushintensity2.setEnabled(true);
		} else if (brush == BRUSH_DIELECTRIC) {
			opts.gui_brushintensity.setText("Dielectric constant: " + getSI(1+intensity, ""));
			opts.gui_brushintensity.setEnabled(true);
			opts.gui_brushintensity2.setEnabled(true);
		} else if (brush == BRUSH_PARAMAGNET) {
			opts.gui_brushintensity.setText("Relative permeability: " + getSI(1+intensity, ""));
			opts.gui_brushintensity.setEnabled(true);
			opts.gui_brushintensity2.setEnabled(true);
		} else if (brush == BRUSH_CHARGE) {
			opts.gui_brushintensity.setText("Charge density: " + getSI(intensity, "C/m^3"));
			opts.gui_brushintensity.setEnabled(true);
			opts.gui_brushintensity2.setEnabled(true);
			
		}
		
		double angle = 0;
		if (brush != BRUSH_EMF) {
			opts.gui_direction.setEnabled(false);
			opts.gui_label_direction.setEnabled(false);
			opts.gui_emf_freq.setEnabled(false);
			opts.gui_acfreqlbl.setEnabled(false);
			opts.gui_label_direction.setText("Brush orientation");
			opts.gui_acfreqlbl.setText("AC Frequency");
		} else {
			opts.gui_direction.setEnabled(true);
			opts.gui_label_direction.setEnabled(true);
			opts.gui_emf_freq.setEnabled(true);
			opts.gui_acfreqlbl.setEnabled(true);

			int directionval = opts.gui_direction.getValue();
			opts.gui_label_direction.setText("Brush orientation: " + directionval*(360/24) + " deg");
			if (opts.gui_emf_freq.getValue() == 0) {
				opts.gui_acfreqlbl.setText("AC Frequency: 0 (DC)");
			} else {
				opts.gui_acfreqlbl.setText("AC Frequency: " + getSI(opts.gui_emf_freq.getValue()/(2*Math.PI), "Hz"));
			}
			angle = Math.PI * directionval/12.0;
		}
		
		if (mouseIsPressed) {

			if (brush == BRUSH_INTERACT) {
				if (rising) {
					int mi = (int)Math.round(((mouseX+1)/(double)scalefactor) - 0.5);
					int mj = (int)Math.round(((mouseY+1)/(double)scalefactor) - 0.5);

					if (materials[mi][mj].type == Material.SWITCH || materials[mi][mj].type == Material.EMFSOURCE) {

						for (int i = 0; i < nx; i++)
						{
							for (int j = 0; j < ny; j++)
							{
								visited[i][j] = false;
							}
						}

						try {
							floodFill(mi, mj, !materials[mi][mj].activated, materials[mi][mj].type);
						}
						catch (StackOverflowError e) {
							JOptionPane.showMessageDialog(opts,
									"Error: Stack overflow in handleMouseInput.");
						}
					}
				}
			} else {
				
				for (int i = 1; i < nx-1; i++)
				{
					for (int j = 1; j < ny-1; j++)
					{
						double cx = 0;
						double cy = 0;
						cx = (i+0.5)*ds;
						cy = (j+0.5)*ds;

						Vector a = new Vector(mxp, myp);
						Vector b = new Vector(mx, my);
						Vector p = new Vector(cx, cy);
						p.addmult(a, -1);
						Vector ab = b.copy();
						ab.addmult(a, -1);

						double l2 = ab.dot(ab);
						if (l2 == 0)
							l2 = 1;
						double t = clamp(p.dot(ab)/l2, 0, 1);
						ab.scalarmult(t);
						p.addmult(ab, -1);
						//double r = Math.sqrt((cx-mx)*(cx-mx) + (cy-my)*(cy-my));
						double r = 0;
						if (brushshape == BRUSH_SHAPE_CIRCLE)
							r = Math.sqrt(p.dot(p));
						else if (brushshape == BRUSH_SHAPE_SQUARE)
							r = Math.max(Math.abs(p.x), Math.abs(p.y));
						if (r <= brushsize) {
							if (brush == BRUSH_CHARGE) {
								double drho = 0.0;
								if (mouseButton == MouseEvent.BUTTON1)
									drho = intensity;
								else if (mouseButton == MouseEvent.BUTTON3)
									drho = -intensity;
								rho[i][j] += drho*Math.exp(-3.0*r/brushsize)/targetframerate;
							} else {
								if (brush == BRUSH_CONDUCTOR) {
									if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.VACUUM) {
										materials[i][j].reset();
										materials[i][j].type = Material.CONDUCTOR;
										materials[i][j].sigma = sigma_max;
										materials[i][j].D = D_max;
									}
								} else if (brush == BRUSH_EMF) {
									if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.VACUUM) {
										materials[i][j].reset();
										materials[i][j].type = Material.EMFSOURCE;
										materials[i][j].emfx = intensity*Math.cos(angle);
										materials[i][j].emfy = intensity*Math.sin(-angle);
										materials[i][j].sigma = sigma_max;
										materials[i][j].D = D_max;
										materials[i][j].frequency = opts.gui_emf_freq.getValue();
									}
								} else if (brush == BRUSH_DIELECTRIC) {
									if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.VACUUM) {
										materials[i][j].reset();
										materials[i][j].type = Material.DIELECTRIC;
										materials[i][j].epsr = 1+intensity;
									}
								} else if (brush == BRUSH_PARAMAGNET) {
									if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.VACUUM) {
										materials[i][j].reset();
										materials[i][j].type = Material.PARAMAGNET;
										materials[i][j].mur = 1+intensity;
									}
								} else if (brush == BRUSH_SWITCH) {
									if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.VACUUM) {
										materials[i][j].reset();
										materials[i][j].type = Material.SWITCH;
										materials[i][j].sigma = sigma_max;
										materials[i][j].D = D_max;
									}
								} 

								if (mouseButton == MouseEvent.BUTTON3) {
									materials[i][j].reset();
								}

								updateMaterial(0, i, j);
								mu_rz[i][j] = materials[i][j].mur;
							}
						}
					}
				}
			}

		}
		mxp = mx;
		myp = my;
	}
	
	public void updateMaterial(int orientation, int i0, int j0) {
		for (int a = 0; a < 4; a++)
		{
			int i = 0;
			int j = 0;
			if (a==0) {
				i = i0;
				j = j0;
				orientation = 0;
			} else if (a==1) {
				i = i0;
				j = j0+1;
				orientation = 0;
			} else if (a==2) {
				i = i0;
				j = j0;
				orientation = 1;
			} else if (a==3) {
				i = i0+1;
				j = j0;
				orientation = 1;
			}
			if (orientation == 0) {
				int i1 = i;
				int j1 = j-1;
				int i2 = i;
				int j2 = j;

				if (j1 < 0)
					j1 = 0;
				if (j2 == ny-1)
					j2 = ny-2;
				
				diffx[i][j] = Math.max(materials[i1][j1].D*materials[i1][j1].diffusion_on, materials[i2][j2].D*materials[i2][j2].diffusion_on);
				emfx[i][j] = Math.max(materials[i1][j1].emfx*materials[i1][j1].emf_on, materials[i2][j2].emfx*materials[i2][j2].emf_on);
				sigmax[i][j] = Math.max(materials[i1][j1].sigma*materials[i1][j1].conductivity_on, materials[i2][j2].sigma*materials[i2][j2].conductivity_on);
				epsrx[i][j] = Math.max(materials[i1][j1].epsr, materials[i2][j2].epsr);
				freqx[i][j] = Math.max(materials[i1][j1].frequency,materials[i2][j2].frequency);
			} else if (orientation == 1) {
				int i1 = i;
				int j1 = j;
				int i2 = i-1;
				int j2 = j;

				if (i1 < 0)
					i1 = 0;
				if (i2 == nx-1)
					i2 = nx-2;

				diffy[i][j] = Math.max(materials[i1][j1].D*materials[i1][j1].diffusion_on, materials[i2][j2].D*materials[i2][j2].diffusion_on);
				emfy[i][j] = Math.max(materials[i1][j1].emfy*materials[i1][j1].emf_on, materials[i2][j2].emfy*materials[i2][j2].emf_on);
				sigmay[i][j] = Math.max(materials[i1][j1].sigma*materials[i1][j1].conductivity_on, materials[i2][j2].sigma*materials[i2][j2].conductivity_on);
				epsry[i][j] = Math.max(materials[i1][j1].epsr, materials[i2][j2].epsr);
				freqy[i][j] = Math.max(materials[i1][j1].frequency,materials[i2][j2].frequency);
			}
		}
	}

	public void updateMainFields() {

		t6.start();
		
		stepnumber++;

		/*Interior B field & charge*/
		
		for (int i = 1; i < nx-2; i++)
		{
			for (int j = 1; j < ny-2; j++)
			{
				bz[i][j] = bz[i][j] + (-(ey[i+1][j] - ey[i][j]) + (ex[i][j+1] - ex[i][j]))*dt/ds;
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				rho[i][j] = rho[i][j] - (jx[i][j]-jx[i-1][j] + jy[i][j]-jy[i][j-1])*dt/ds;
			}
		}

		/*Boundary*/

		for (int j = 1; j < ny-2; j++)
		{
			bz[0][j] = bzp[1][j] + (tauminus/tauplus)*(bzp[0][j]-bz[1][j]) + (dt/(2*tauplus*ds)) * (ex[0][j+1] - ex[0][j] + ex[1][j+1] - ex[1][j]);
			bz[nx-2][j] = bzp[nx-3][j] + (tauminus/tauplus)*(bzp[nx-2][j]-bz[nx-3][j]) + (dt/(2*tauplus*ds)) * (ex[nx-2][j+1] - ex[nx-2][j] + ex[nx-3][j+1] - ex[nx-3][j]);
		}

		for (int i = 1; i < nx-2; i++)
		{
			bz[i][0] = bzp[i][1] + (tauminus/tauplus)*(bzp[i][0]-bz[i][1]) - (dt/(2*tauplus*ds)) * (ey[i+1][0] - ey[i][0] + ey[i+1][1] - ey[i][1]);
			bz[i][ny-2] = bzp[i][ny-3] + (tauminus/tauplus)*(bzp[i][ny-2]-bz[i][ny-3]) - (dt/(2*tauplus*ds)) * (ey[i+1][ny-2] - ey[i][ny-2] + ey[i+1][ny-3] - ey[i][ny-3]);
		}

		/* Corners */
		bz[0][0] = 0.25*(bz[1][0] + bz[0][1]) + 0.5*bzp[0][0];
		bz[nx-2][0] = 0.25*(bz[nx-3][0] + bz[nx-2][1]) + 0.5*bzp[nx-2][0];
		bz[0][ny-2] = 0.25*(bz[0][ny-3] + bz[1][ny-2]) + 0.5*bzp[0][ny-2];
		bz[nx-2][ny-2] = 0.25*(bz[nx-3][ny-2] + bz[nx-2][ny-3]) + 0.5*bzp[nx-2][ny-2];

		for (int i = 0; i <= 1; i++)
			for (int j = 0; j < ny; j++)
				bzp[i][j] = bz[i][j];

		for (int i = nx-3; i <= nx-2; i++)
			for (int j = 0; j < ny; j++)
				bzp[i][j] = bz[i][j];

		for (int i = 0; i < nx; i++)
			for (int j = 0; j <= 1; j++)
				bzp[i][j] = bz[i][j];

		for (int i = 0; i < nx; i++)
			for (int j = ny-3; j <= ny-2; j++)
				bzp[i][j] = bz[i][j];

		
		/*E field and current*/
		//10.0*Math.pow(10.0, opts.gui_emf_freq.getValue()/15.0)
		for (int i = 0; i < maxfreq; i++) {
			amplitudes[i] = Math.cos(i*time);
		}

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				jx[i][j] = -diffx[i][j]*(rho[i+1][j] - rho[i][j])/ds + sigmax[i][j]*(emfx[i][j]*amplitudes[freqx[i][j]]+0.5*ex[i][j]);
				ex[i][j] = (ex[i][j] + (csq*(bz[i][j]/mu_rz[i][j]-bz[i][j-1]/mu_rz[i][j-1])/ds - jx[i][j]/eps0)*(dt/epsrx[i][j]))
						/(1+0.5*dt*sigmax[i][j]/(eps0*epsrx[i][j]));
				jx[i][j] += 0.5*sigmax[i][j]*ex[i][j];
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				jy[i][j] = -diffy[i][j]*(rho[i][j+1] - rho[i][j])/ds + sigmay[i][j]*(emfy[i][j]*amplitudes[freqy[i][j]]+0.5*ey[i][j]);
				ey[i][j] = (ey[i][j] + (-csq*(bz[i][j]/mu_rz[i][j]-bz[i-1][j]/mu_rz[i-1][j])/ds - jy[i][j]/eps0)*(dt/epsry[i][j]))
						/(1+0.5*dt*sigmay[i][j]/(eps0*epsry[i][j]));
				jy[i][j] += 0.5*sigmay[i][j]*ey[i][j];
			}
		}
		
		t6.stop();

		if (stepnumber%40 == 0) {
			correctEfield();
		}

		time += dt;
		advanceframe = false;
	}

	public void calcAuxillaryFields() {
		if (opts.gui_view.getSelectedIndex() == VIEWSCALAR_H) {
			for (int i = 1; i < nx-2; i++)
			{
				for (int j = 1; j < ny-2; j++)
				{
					hz[i][j] = bz[i][j]/(mu_rz[i][j]*mu0);
				}
			}
		}

		else if (opts.gui_view.getSelectedIndex() == VIEWSCALAR_U) {
			for (int i = 1; i < nx-2; i++)
			{
				for (int j = 1; j < ny-2; j++)
				{
					u[i][j] = (ex[i][j]*ex[i][j]*epsrx[i][j]
							+ ex[i][j+1]*ex[i][j+1]*epsrx[i][j+1]
									+ ey[i][j]*ey[i][j]*epsry[i][j]
											+ey[i+1][j]*ey[i+1][j]*epsry[i+1][j])*eps0/4.0
							+ bz[i][j]*bz[i][j]/(2*mu_rz[i][j]*mu0);
				}
			}
		}

		if (opts.gui_view_vec.getSelectedIndex() == VIEWVEC_D) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					dx[i][j] = epsrx[i][j]*ex[i][j]*eps0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					dy[i][j] = epsry[i][j]*ey[i][j]*eps0;
				}
			}
		}

		else if (opts.gui_view_vec.getSelectedIndex() == VIEWVEC_S) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					sy[i][j] = -ex[i][j]*0.5*(bz[i][j]/mu_rz[i][j] + bz[i][j-1]/mu_rz[i][j-1])/mu0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					sx[i][j] = ey[i][j]*0.5*(bz[i][j]/mu_rz[i][j] + bz[i-1][j]/mu_rz[i-1][j])/mu0;
				}
			}
		}

	}
	
	public void prescaleDielectric() {
		int maxfineness = (int)Math.floor(Math.log(nx)/Math.log(2));
		for (int fineness = 2; fineness <= maxfineness; fineness++) {
			int nxcgint = (1 << fineness) - 1;
			int nycgint = (1 << fineness) - 1;
			double gridsize = width/(1 << fineness);
			downscale(epsrx, epsxMG[fineness], nx-2, ny-1, nxcgint, nycgint+1, ds/2.0, 0, ds, gridsize/2.0, 0, gridsize);
			downscale(epsry, epsyMG[fineness], nx-1, ny-2, nxcgint+1, nycgint, 0, ds/2.0, ds, 0, gridsize/2.0, gridsize);
		}
	}
	
	public void correctEfield() {
		boolean debug = false;
		t4.start();
		t3.start();
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				poisson1[i][j] = 0;
				//poisson2[i][j] = 0;
				deltarhoMG[i][j] = 0;
			}
		}

		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				deltarho[i][j] = eps0*((ex[i][j]*epsrx[i][j]-ex[i-1][j]*epsrx[i-1][j] + ey[i][j]*epsry[i][j]-ey[i][j-1]*epsry[i][j-1])/ds) - rho[i][j];
			}
		}
		
		double err = 0;
		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				err += Math.abs(((epsrx[i-1][j]+epsrx[i][j]+epsry[i][j-1]+epsry[i][j])*poisson1[i][j]
						- (poisson1[i-1][j]*epsrx[i-1][j] + poisson1[i+1][j]*epsrx[i][j] + poisson1[i][j-1]*epsry[i][j-1] + poisson1[i][j+1]*epsry[i][j]))
						/(ds*ds) - deltarho[i][j]/eps0);
			}
		}
		System.out.println("Poisson error: " + err);
		
		int nxint = nx-2;
		int nyint = ny-2;
		
		int[] stepsarray = {0, 0, 200, 200, 200, 200, 200, 50, 20, 20, 20};
		
		t3.stop("Init");
		
		int maxfineness = (int)Math.floor(Math.log(nx)/Math.log(2));
		for (int fineness = 2; fineness <= maxfineness; fineness++) {
			int nxcgint = (1 << fineness) - 1;
			int nycgint = (1 << fineness) - 1;
			double gridsize = width/(1 << fineness);
			

			t3.start();
			downscale(deltarho, deltarhoMG, nx-1, ny-1, nxcgint+1, nycgint+1, 0, 0, ds, 0, 0, gridsize);
			t3.stop("Charge computation");
			
			//public void downscale(double[][] source, double[][] dest, int sx, int sy, int dx, int dy, double soffset, double sspacing, double doffset, double dspacing) {
			
			int poissonsteps = stepsarray[fineness];
			double alpha = (gridsize*gridsize)/eps0;
			double beta = 4;

			t3.start();
			JacobiIteration(poissonsteps, nxcgint, nycgint, alpha, beta, gridsize, fineness, debug);
			t3.stop("Iteration");

			t3.start();
			if (fineness < maxfineness) {

				for (int i = 0; i < nxcgint+2; i++) {
					for (int j = 0; j < nycgint+2; j++) {
						poisson2[i][j] = poisson1[i][j];
					}
				}

				nxcgint = (1 << (fineness+1)) - 1;
				nycgint = (1 << (fineness+1)) - 1;
				gridsize = width/(1 << (fineness+1));


				for (int i = 1; i < nxcgint+1; i++)
				{
					for (int j = 1; j < nycgint+1; j++) {
						poisson1[i][j] = bilinearinterp(poisson2, i/2.0, j/2.0);
					}
				}


				for (int i = 0; i < nxcgint+2; i++) {
					poisson2[i][0] = 0;
					poisson2[i][nycgint+1] = 0;
				}
				for (int j = 0; j < nycgint+2; j++) {
					poisson2[0][j] = 0;
					poisson2[nxcgint+1][j] = 0;
				}
			} else {
				for (int i = 0; i < nxcgint+2; i++) {
					for (int j = 0; j < nycgint+2; j++) {
						poisson2[i][j] = poisson1[i][j];
					}
				}

				for (int i = 1; i < nxint+1; i++)
				{
					for (int j = 1; j < nyint+1; j++) {
						poisson1[i][j] = bilinearinterp(poisson2, i*(double)(nxcgint+1.0)/((double)nx-1.0), j*(double)(nycgint+1.0)/((double)ny-1.0));
					}
				}
				
			}
			t3.stop("Upscaling");
		}
		

		int poissonsteps = stepsarray[maxfineness+1];

		double alpha = (ds*ds)/eps0;
		double beta = 4;

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				deltarhoMG[i][j] = deltarho[i][j];
				epsxMG[maxfineness+1][i][j] = epsrx[i][j];
				epsyMG[maxfineness+1][i][j] = epsry[i][j];
			}
		}
		
		t3.start();
		JacobiIteration(poissonsteps, nxint, nyint, alpha, beta, ds, maxfineness+1, debug);
		
		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				ex[i][j] = ex[i][j] + (poisson1[i+1][j]-poisson1[i][j])/ds;
				ey[i][j] = ey[i][j] + (poisson1[i][j+1]-poisson1[i][j])/ds;
			}
		}

		t3.stop("Final iteration");

		err = 0;
		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				err += Math.abs(((epsrx[i-1][j]+epsrx[i][j]+epsry[i][j-1]+epsry[i][j])*poisson1[i][j]
						- (poisson1[i-1][j]*epsrx[i-1][j] + poisson1[i+1][j]*epsrx[i][j] + poisson1[i][j-1]*epsry[i][j-1] + poisson1[i][j+1]*epsry[i][j]))
						/(ds*ds) - deltarho[i][j]/eps0);
			}
		}
		System.out.println("Poisson error: " + err);
		t4.stop();
		//t3.stop();
	}
	

	public void JacobiIteration(int steps, int xbound, int ybound, double alpha, double beta, double gridsize, int fineness, boolean debug) {

		for (int i = 1; i < xbound+1; i++) {
			for (int j = 1; j < ybound+1; j++) {
				epsavgMG[i][j] = (epsxMG[fineness][i-1][j]+epsxMG[fineness][i][j]+epsyMG[fineness][i][j-1]+epsyMG[fineness][i][j]);
			}
		}
		
		for (int poissonit = 0; poissonit < steps; poissonit++) {
			for (int i = 1; i < xbound+1; i++) {
				for (int j = 1; j < ybound+1; j++) {
					poisson2[i][j] = ((poisson1[i-1][j]*epsxMG[fineness][i-1][j]
							+ poisson1[i+1][j]*epsxMG[fineness][i][j]
									+ poisson1[i][j-1]*epsyMG[fineness][i][j-1]
											+ poisson1[i][j+1]*epsyMG[fineness][i][j]) + deltarhoMG[i][j]*alpha)/epsavgMG[i][j];
				}
			}
			for (int i = 1; i < xbound+1; i++) {
				for (int j = 1; j < ybound+1; j++) {
					poisson1[i][j] = ((poisson2[i-1][j]*epsxMG[fineness][i-1][j]
							+ poisson2[i+1][j]*epsxMG[fineness][i][j]
									+ poisson2[i][j-1]*epsyMG[fineness][i][j-1]
											+ poisson2[i][j+1]*epsyMG[fineness][i][j]) + deltarhoMG[i][j]*alpha)/epsavgMG[i][j];
				}
			}

			if (debug) {
				for (int i = 1; i < xbound+1; i++) {
					for (int j = 1; j < ybound+1; j++) {
						bz[i][j] = ((epsxMG[fineness][i-1][j]+epsxMG[fineness][i][j]+epsyMG[fineness][i][j-1]+epsyMG[fineness][i][j])*poisson1[i][j]
								- (poisson1[i-1][j]*epsxMG[fineness][i-1][j] + poisson1[i+1][j]*epsxMG[fineness][i][j] + poisson1[i][j-1]*epsyMG[fineness][i][j-1] + poisson1[i][j+1]*epsyMG[fineness][i][j]))
								/(gridsize*gridsize) - deltarhoMG[i][j]/eps0;
					}
				}

				r.repaint();
				try {
					Thread.sleep(17);
				} catch (InterruptedException e) {}
			}
		}
	}
	
	public void downscale(double[][] source, double[][] dest, int sx, int sy, int dx, int dy, double soffsetx, double soffsety, double sspacing, double doffsetx, double doffsety, double dspacing) {
		for (int i = 0; i <= dx; i++) {
			for (int j = 0; j < dy; j++) {
				dest[i][j] = 0;
			}
		}
		double scalefactor = (sspacing*sspacing/(dspacing*dspacing));
		for (int i = 0; i <= sx; i++) {
			for (int j = 0; j <= sy; j++) {
				double cx = (i*sspacing + soffsetx - doffsetx)/dspacing;
				double cy = (j*sspacing + soffsety - doffsety)/dspacing;
				int xfloor = (int)Math.floor(cx);
				int yfloor = (int)Math.floor(cy);
				double fx = cx - xfloor;
				double fy = cy - yfloor;
				if (xfloor < 0) {
					xfloor = 0;
					fx = 0.0;
				} else if (xfloor >= dx) {
					xfloor = dx - 1;
					fx = 1.0;
				}
				if (yfloor < 0) {
					yfloor = 0;
					fy = 0.0;
				} else if (yfloor >= dy) {
					yfloor = dy - 1;
					fy = 1.0;
				}
				double srcval = source[i][j]*scalefactor;
				dest[xfloor][yfloor] += srcval*(1-fx)*(1-fy);
				dest[xfloor+1][yfloor] += srcval*fx*(1-fy);
				dest[xfloor][yfloor+1] += srcval*(1-fx)*fy;
				dest[xfloor+1][yfloor+1] += srcval*fx*fy;
			}
		}
	}
	
	public double bilinearinterp(double[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		
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

	public void floodFill(int i, int j, boolean active, int type) {
		if (i >= 0 && i < nx && j > 0 && j < ny && !visited[i][j] && materials[i][j].type == type) {
			materials[i][j].activated = active;
			materials[i][j].diffusion_on = ((materials[i][j].type == Material.SWITCH) & !active) ? 0:1;
			materials[i][j].conductivity_on = ((materials[i][j].type == Material.SWITCH) & !active) ? 0:1;
			materials[i][j].emf_on = ((materials[i][j].type == Material.EMFSOURCE) & !active) ? 0:1;
			visited[i][j] = true;
			
			floodFill(i-1, j, active, type);
			floodFill(i+1, j, active, type);
			floodFill(i, j-1, active, type);
			floodFill(i, j+1, active, type);

			updateMaterial(0, i, j);
			updateMaterial(0, i, j+1);
			updateMaterial(1, i, j);
			updateMaterial(1, i+1, j);
		}
	}

	public double getex(double x, double y) {
		return bilinearinterp(ex, x - 0.5, y);
	}
	
	public double getey(double x, double y) {
		return bilinearinterp(ey, x, y - 0.5);
	}
	
	public double getbz(double x, double y) {
		return bilinearinterp(bz, x - 0.5, y - 0.5);
	}
	
	public double getrho(double x, double y) {
		return bilinearinterp(rho, x, y);
	}
	
	public double getjx(double x, double y) {
		return bilinearinterp(jx, x - 0.5, y);
	}
	
	public double getjy(double x, double y) {
		return bilinearinterp(jy, x, y - 0.5);
	}

	public double getdx(double x, double y) {
		return bilinearinterp(dx, x - 0.5, y);
	}
	
	public double getdy(double x, double y) {
		return bilinearinterp(dy, x, y - 0.5);
	}

	public double getsx(double x, double y) {
		return bilinearinterp(sx, x, y - 0.5);
	}
	
	public double getsy(double x, double y) {
		return bilinearinterp(sy, x - 0.5, y);
	}
	
	public void fillRectangle(int x, int y, int w, int h) {
		if (x < 0) {
			x=0;
		} 
		if (x+w > screen.getWidth()) {
			w = screen.getWidth()-x-1;
		}
		if (y < 0) {
			y=0;
		} 
		if (y+h > screen.getHeight()) {
			h = screen.getHeight()-y-1;
		}
		
		for (int i = x; i < x+w; i++) {
			for (int j = y; j < y+h; j++) {
				int rgb = screen.getRGB(i, j);
				int r = ((int)(((rgb>>16)&255)*alphaBG + col_r*alphaFG));
				int g = ((int)(((rgb>>8)&255)*alphaBG + col_g*alphaFG));
				int b = ((int)(((rgb)&255)*alphaBG + col_b*alphaFG));
				if (r > 255)
					r = 255;
				if (g > 255)
					g = 255;
				if (b > 255)
					b = 255;
				screen.setRGB(i, j, 255<<24 | r << 16 | g << 8 | b);
			}
		}
	}
	public void drawLine(int x0, int y0, int x1, int y1) {
		try {
			int dy = y1 - y0;
			int dx = x1 - x0;
			float t = (float) 0.5; // offset for rounding


			int rgb = screen.getRGB(x0, y0);
			int r = ((int)(((rgb>>16)&255)*alphaBG + col_r*alphaFG));
			int g = ((int)(((rgb>>8)&255)*alphaBG + col_g*alphaFG));
			int b = ((int)(((rgb)&255)*alphaBG + col_b*alphaFG));
			if (r > 255)
				r = 255;
			if (g > 255)
				g = 255;
			if (b > 255)
				b = 255;
			screen.setRGB(x0, y0, 255<<24 | r << 16 | g << 8 | b);


			if (Math.abs(dx) > Math.abs(dy)) { // slope < 1
				float m = (float) dy / (float) dx; // compute slope
				t += y0;
				dx = (dx < 0) ? -1 : 1;
				m *= dx;
				while (x0 != x1) {
					x0 += dx; // step to next x value
					t += m; // add slope to y value

					rgb = screen.getRGB(x0, (int)t);
					r = ((int)(((rgb>>16)&255)*alphaBG + col_r*alphaFG));
					g = ((int)(((rgb>>8)&255)*alphaBG + col_g*alphaFG));
					b = ((int)(((rgb)&255)*alphaBG + col_b*alphaFG));
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					screen.setRGB(x0, (int)t, 255<<24 | r << 16 | g << 8 | b);

				}
			} else { // slope >= 1
				float m = (float) dx / (float) dy; // compute slope
				t += x0;
				dy = (dy < 0) ? -1 : 1;
				m *= dy;
				while (y0 != y1) {
					y0 += dy; // step to next y value
					t += m; // add slope to x value

					rgb = screen.getRGB((int)t, y0);
					r = ((int)(((rgb>>16)&255)*alphaBG + col_r*alphaFG));
					g = ((int)(((rgb>>8)&255)*alphaBG + col_g*alphaFG));
					b = ((int)(((rgb)&255)*alphaBG + col_b*alphaFG));
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					screen.setRGB((int)t, y0, 255<<24 | r << 16 | g << 8 | b);
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			return;
		}
	}

	public void setColor(int r, int g, int b) {
		col_r = r;
		col_g = g;
		col_b = b;
	}

	public void setalphaBG(double alpha) {
		alphaBG = alpha;
	}
	public void setalphaFG(double alpha) {
		alphaFG = alpha;
	}
	public int col_r = 0;
	public int col_g = 0;
	public int col_b = 0;
	public double alphaBG = 0;
	public double alphaFG = 0;
	
	public void render() {
		t5.start();
		Graphics2D g = (Graphics2D) screen.getGraphics();
		g.setColor(new Color(50, 50, 50));
		g.fillRect(0, 0, imgwidth, imgheight);
		double scalingconstant = 10.0*Math.pow(10.0, opts.gui_brightness.getValue()/10.0);
		int scalarview = opts.gui_view.getSelectedIndex();
		int vectorview = opts.gui_view_vec.getSelectedIndex();
		setalphaBG(0);
		setalphaFG(1);
		for (int i = 0; i < nx-1; i++) {
			for (int j = 0; j < ny-1; j++) {
				double x = i;
				double y = j;

				if (scalarview == VIEWSCALAR_NONE) {
					setcol(0, 0, 0, 30);
				}
				else if (scalarview == VIEWSCALAR_E) {
					double exg = 0.5*(ex[i][j]+ex[i][j+1])*scalingconstant;
					double eyg = 0.5*(ey[i][j]+ey[i+1][j])*scalingconstant;
					setcol(Math.abs(exg), 0, Math.abs(eyg), 30);
				}
				else if (scalarview == VIEWSCALAR_B) {
					double bzg = this.bz[i][j]*c*scalingconstant;
					setcol(clamp(bzg, 0, 999), 0, clamp(-bzg, 0, 999), 30);
				}
				else if (scalarview == VIEWSCALAR_RHO) {
					//x -= 0.5;
					//y -= 0.5;
					double rhog = (rho[i][j]+rho[i+1][j]+rho[i][j+1]+rho[i+1][j+1])*2*ds*scalingconstant;
					setcol(clamp(rhog, 0, 999), 0, clamp(-rhog, 0, 999), 30);
				}
				else if (scalarview == VIEWSCALAR_J) {
					double jxg = 0.5*(jx[i][j]+jx[i][j+1])*scalingconstant/sigma_max;
					double jyg = 0.5*(jy[i][j]+jy[i+1][j])*scalingconstant/sigma_max;
					setcol(Math.abs(jxg), 0, Math.abs(jyg), 30);
				}
				else if (scalarview == VIEWSCALAR_H) {
					double hzg = this.hz[i][j]*c*mu0*scalingconstant;
					setcol(clamp(hzg, 0, 999), 0, clamp(-hzg, 0, 999), 30);
				}
				else if (scalarview == VIEWSCALAR_U) {
					double ug = this.u[i][j]*scalingconstant*5;
					setcol(0, ug, 0, 30);
				}
				
				fillRectangle((int)(x*scalefactor), (int)(y*scalefactor), scalefactor, scalefactor);
			}
		}

		setalphaBG(1);
		setalphaFG(0.5);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (materials[i][j].type == Material.CONDUCTOR) {
					setColor(180, 180, 180);
				} else if (materials[i][j].type == Material.DIELECTRIC) {
					setColor(96, 163, 77);
				} else if (materials[i][j].type == Material.PARAMAGNET) {
					setColor(74, 167, 147);
				} else if (materials[i][j].type == Material.EMFSOURCE) {
					if (materials[i][j].activated) {
						setColor(255, 118, 140);
					} else {
						setColor(146, 53, 74);
					}
				} else if (materials[i][j].type == Material.SWITCH) {
					if (materials[i][j].activated) {
						setColor(247, 209, 108);
					} else {
						setColor(128, 110, 91);
					}
				}

				if (materials[i][j].type != Material.VACUUM)
					fillRectangle((int)(i*scalefactor), (int)(j*scalefactor), scalefactor, scalefactor);
			}
		}

		setalphaBG(1);
		setalphaFG(0.3);
		
		double arrowlength = 10.0/scalefactor;

		double vectorscalingconstant = 10.0*Math.pow(10.0, opts.gui_brightness_1.getValue()/10.0);
		
		//boolean testFieldLines = false;
		
		for (int i = 0; i < 50; i++) {
			for (int j = 0; j < 50; j++) {
				double x = (nx-1)*(i+0.5)/50;
				double y = (ny-1)*(j+0.5)/50;
				Vector ctr = new Vector(x, y);
				Vector arrow = null;

				if (vectorview == VIEWVEC_NONE)
					arrow = new Vector(0, 0);
				else if (vectorview == VIEWVEC_E)
					arrow = new Vector(getex(x, y), getey(x, y));
				else if (vectorview == VIEWVEC_J)
					arrow = new Vector(getjx(x, y)/sigma_max, getjy(x, y)/sigma_max);
				else if (vectorview == VIEWVEC_D)
					arrow = new Vector(getdx(x, y)/eps0, getdy(x, y)/eps0);
				else if (vectorview == VIEWVEC_S)
					arrow = new Vector(getsx(x, y), getsy(x, y));

				double fieldmagnitude = vectorscalingconstant*Math.sqrt(arrow.dot(arrow));
				arrow.normalize();
				Vector tip1 = arrow.copy();
				Vector tip2 = arrow.copy();
				tip1.rotate(Math.PI*5.0/6.0);
				tip2.rotate(Math.PI*7.0/6.0);
				
				Vector body1 = ctr.copy();
				body1.addmult(arrow, -0.5*arrowlength);
				Vector body2 = ctr.copy();
				body2.addmult(arrow, 0.5*arrowlength);
				tip1.scalarmult(0.2*arrowlength);
				tip1.add(body2);
				tip2.scalarmult(0.2*arrowlength);
				tip2.add(body2);
				setcol(fieldmagnitude, fieldmagnitude, fieldmagnitude, 30);
				drawLine((int)(body1.x*scalefactor), (int)(body1.y*scalefactor), (int)(body2.x*scalefactor), (int)(body2.y*scalefactor));
				drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip1.x*scalefactor), (int)(tip1.y*scalefactor));
				drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip2.x*scalefactor), (int)(tip2.y*scalefactor));
				//tip1.add(ctr);
			}
		}

		((Graphics2D)g).setRenderingHint(
				RenderingHints.KEY_TEXT_ANTIALIASING,
				RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		double brushsize = 1.5*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
		int r = (int)(scalefactor*brushsize/ds);
		int brushshape = opts.gui_brush_1.getSelectedIndex();
		if (opts.gui_brush.getSelectedIndex() != BRUSH_INTERACT)
		if (brushshape == 0) {
			g.setColor(new Color(50, 50, 50));
			g.drawOval(mouseX - r+1, mouseY - r+1, 2*r, 2*r);
			g.setColor(new Color(200, 200, 200));
			g.drawOval(mouseX - r, mouseY - r, 2*r, 2*r);
		} else {
			g.setColor(new Color(50, 50, 50));
			g.drawRect(mouseX - r+1, mouseY - r+1, 2*r-2, 2*r-2);
			g.setColor(new Color(200, 200, 200));
			g.drawRect(mouseX - r, mouseY - r, 2*r, 2*r);
		}
		
		
		g.setColor(Color.WHITE);


		if (opts.gui_tooltip.isSelected()) {
			double mx = (mouseX/(double)scalefactor);
			double my = (mouseY/(double)scalefactor);
			int mi = (int)Math.floor(mx);
			int mj = (int)Math.floor(my);
			
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
			int voffset = 1 + mouseY;
			int hoffset = 5 + mouseX+15;
			double eg = length(getex(mx,my), getey(mx,my));
			double bg = getbz(mx, my);
			double rhog = getrho(mx, my);
			double jg = length(getjx(mx,my), getjy(mx,my));
			
			String name = "Material: " + Material.NAMES[mat.type];
			if (mat.type == Material.EMFSOURCE || mat.type == Material.SWITCH) {
				name = name + " (" + (mat.activated? "on" : "off") + ")";
			}
			
			g.setFont(bigfont);
			drawString(name, hoffset, voffset + 1*vspacing, g);
			voffset = voffset+3;
			g.setFont(regularfont);
			drawString("E: " + getSI(eg, "V/m"), hoffset, voffset + 2*vspacing, g);
			drawString("B: " + getSI(bg, "T"), hoffset, voffset + 3*vspacing, g);
			drawString("\u03c1: " + getSI(rhog, "C/m^3"), hoffset, voffset + 4*vspacing, g);
			drawString("J: " + getSI(jg, "A/m^2"), hoffset, voffset + 5*vspacing, g);
			drawString("\u03c3: " + getSI(mat.sigma, "S"), hoffset, voffset + 6*vspacing, g);
			drawString("emf: " + getSI(length(mat.emfx, mat.emfy), "V/m"), hoffset, voffset + 7*vspacing, g);
			drawString("\u03b5r: " + getSI(mat.epsr, ""), hoffset, voffset + 8*vspacing, g);
			drawString("\u03bcr: " + getSI(mat.mur, ""), hoffset, voffset + 9*vspacing, g);
		}
		if (opts.gui_paused.isSelected())
			drawString("Paused", 5, 15, g);
		t5.stop();
		//g.setColor(.)
	}
	
	public void drawString(String str, int x, int y, Graphics g) {
		g.setColor(Color.DARK_GRAY);
		g.drawString(str, x+1, y+1);
		g.setColor(Color.WHITE);
		g.drawString(str, x, y);
	}

	public void setcol(double r, double g, double b, int base) {
		setColor(clamp((int)(255.0*r), 0, 999) + base, clamp((int)(255.0*g), 0, 999) + base, clamp((int)(255.0*b), 0, 999) + base);
	}
	
	public int clamp(int val, int min, int max) {
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}

	public double clamp(double val, double min, double max) {
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}
	
	public double length(double x, double y) {
		return Math.sqrt(x*x+y*y);
	}

	public static String getSI(double quantity, String unit) {
		double mag = Math.abs(quantity);
		if (mag < 1E-12)
			return String.format("%.1f", quantity*1e15) + " p" + unit;
		else if (mag < 1E-9)
			return String.format("%.1f", quantity*1e12) + " n" + unit;
		else if (mag < 1E-6)
			return String.format("%.1f", quantity*1e9) + " u" + unit;
		else if (mag < 1E-3)
			return String.format("%.1f", quantity*1e6) + " \u00b5" + unit;
		else if (mag < 1)
			return String.format("%.1f", quantity*1e3) + " m" + unit;
		else if (mag < 1E3)
			return String.format("%.1f", quantity) + " " + unit;
		else if (mag < 1E6)
			return String.format("%.1f", quantity*1e-3) + " k" + unit;
		else if (mag < 1E9)
			return String.format("%.1f", quantity*1e-6) + " M" + unit;
		else
			return String.format("%.1f", quantity*1e-9) + " G" + unit;
	}
	
	public boolean readFile()
	{

		JFileChooser fd = new JFileChooser();
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
		fd.showOpenDialog(opts);
		File infile = fd.getSelectedFile();
		if (infile == null) return false;
		try {
			GZIPInputStream fstr = new GZIPInputStream(new FileInputStream(infile));
			ObjectInputStream ostr = new ObjectInputStream(fstr);

			int version = ostr.readInt();

			if (version == 1) {
				time = ostr.readDouble();
				phase = ostr.readDouble();
				opts.gui_paused.setSelected(ostr.readBoolean());
				opts.gui_tooltip.setSelected(ostr.readBoolean());
				opts.gui_view.setSelectedIndex(ostr.readInt());
				opts.gui_view_vec.setSelectedIndex(ostr.readInt());
				opts.gui_simspeed.setValue(ostr.readInt());
				opts.gui_simspeed_2.setValue(ostr.readInt());
				opts.gui_brightness.setValue(ostr.readInt());
				opts.gui_brightness_1.setValue(ostr.readInt());
				opts.gui_emf_freq.setValue(ostr.readInt());
				opts.gui_brush.setSelectedIndex(ostr.readInt());
				opts.gui_brush_1.setSelectedIndex(ostr.readInt());
				opts.gui_brushsize.setValue(ostr.readInt());
				opts.gui_brushintensity2.setValue(ostr.readInt());
				opts.gui_direction.setValue(ostr.readInt());
				opts.textPane.setText(ostr.readUTF());

				(ex) = (double[][]) ostr.readObject();
				(ey) = (double[][]) ostr.readObject();
				(bz) = (double[][]) ostr.readObject();
				(rho) = (double[][]) ostr.readObject();
				(jx) = (double[][]) ostr.readObject();
				(jy) = (double[][]) ostr.readObject();
				(bzp) = (double[][]) ostr.readObject();
				(materials) = (Material[][]) ostr.readObject();
				(diffx) = (double[][]) ostr.readObject();
				(diffy) = (double[][]) ostr.readObject();
				(sigmax) = (double[][]) ostr.readObject();
				(sigmay) = (double[][]) ostr.readObject();
				(emfx) = (double[][]) ostr.readObject();
				(emfy) = (double[][]) ostr.readObject();
				(epsrx) = (double[][]) ostr.readObject();
				(epsry) = (double[][]) ostr.readObject();
				(mu_rz) = (double[][]) ostr.readObject();
				(freqx) = (int[][]) ostr.readObject();
				(freqy) = (int[][]) ostr.readObject();
				(amplitudes) = (double[]) ostr.readObject();
			}

			opts.textPane.setEditable(false);
			opts.textPane.setCaretPosition(0);
			prescaleDielectric();

			ostr.close();
			fstr.close();
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException | ClassNotFoundException | IllegalArgumentException e) {
			JOptionPane.showMessageDialog(opts,
				    "Error: Unable to load file.");
			return false;
		}
		return true;
	}

	public boolean writeFile()
	{
		JFileChooser fd = new JFileChooser();
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
		fd.showSaveDialog(opts);
		File outfile = fd.getSelectedFile();
		if (outfile == null) return false;
		if (!outfile.getName().endsWith(fileextension))
			outfile = new File(outfile.getAbsolutePath() + fileextension);
		try {
			GZIPOutputStream fstr = new GZIPOutputStream(new FileOutputStream(outfile));
			ObjectOutputStream ostr = new ObjectOutputStream(fstr);
			
			ostr.writeInt(1);
			ostr.writeDouble(time);
			ostr.writeDouble(phase);
			
			ostr.writeBoolean(opts.gui_paused.isSelected());
			ostr.writeBoolean(opts.gui_tooltip.isSelected());
			ostr.writeInt(opts.gui_view.getSelectedIndex());
			ostr.writeInt(opts.gui_view_vec.getSelectedIndex());
			ostr.writeInt(opts.gui_simspeed.getValue());
			ostr.writeInt(opts.gui_simspeed_2.getValue());
			ostr.writeInt(opts.gui_brightness.getValue());
			ostr.writeInt(opts.gui_brightness_1.getValue());
			ostr.writeInt(opts.gui_emf_freq.getValue());
			ostr.writeInt(opts.gui_brush.getSelectedIndex());
			ostr.writeInt(opts.gui_brush_1.getSelectedIndex());
			ostr.writeInt(opts.gui_brushsize.getValue());
			ostr.writeInt(opts.gui_brushintensity2.getValue());
			ostr.writeInt(opts.gui_direction.getValue());
			ostr.writeUTF(opts.textPane.getText());

			ostr.writeObject(ex);
			ostr.writeObject(ey);
			ostr.writeObject(bz);
			ostr.writeObject(rho);
			ostr.writeObject(jx);
			ostr.writeObject(jy);
			ostr.writeObject(bzp);
			ostr.writeObject(materials);
			ostr.writeObject(diffx);
			ostr.writeObject(diffy);
			ostr.writeObject(sigmax);
			ostr.writeObject(sigmay);
			ostr.writeObject(emfx);
			ostr.writeObject(emfy);
			ostr.writeObject(epsrx);
			ostr.writeObject(epsry);
			ostr.writeObject(mu_rz);
			ostr.writeObject(freqx);
			ostr.writeObject(freqy);
			ostr.writeObject(amplitudes);
			
			
			ostr.close();
			fstr.close();
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException e) {
			JOptionPane.showMessageDialog(opts,
				    "Error: Unable to save file.");
			return false;
		}
		return true;
	}
	
	
	PointerInfo p = MouseInfo.getPointerInfo();

	@Override
	public void mouseClicked(MouseEvent arg0) {}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent e) {
		mouseIsPressed = true;
		mouseButton = e.getButton();
		mouseX = e.getX();
		mouseY = e.getY();
		if (e.isShiftDown()) {
			linemode = true;
			mxp = ((mouseX+1)/(double)scalefactor)*ds;
			myp = ((mouseY+1)/(double)scalefactor)*ds;
		}else
			linemode = false;
		lineorientation = 0;
		mouseXstart = e.getX();
		mouseYstart = e.getY();
	}
	@Override
	public void mouseReleased(MouseEvent e) {
		mouseIsPressed = false;
		if (e.isShiftDown())
			linemode = true;
		else
			linemode = false;	
		lineorientation = 0;
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		if (!linemode) {
			mouseX = e.getX();
			mouseY = e.getY();
		} else {
			if (lineorientation == 0) {
				int deltax = e.getX() - mouseXstart;
				int deltay = e.getY() - mouseYstart;
				if (Math.abs(deltax) > 10 || Math.abs(deltay) > 10) {
					if (Math.abs(deltax) > Math.abs(deltay))
						lineorientation = 1;
					else
						lineorientation = 2;
				}
			} else {
				if (lineorientation == 1) {
					mouseY = mouseYstart;
					mouseX = e.getX();
				} else if (lineorientation == 2) {
					mouseX = mouseXstart;
					mouseY = e.getY();
				}
			}
		}
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		mouseX = arg0.getX();
		mouseY = arg0.getY();
	}
	
	private Action key_p = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -3637591771664096455L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_paused.setSelected(!opts.gui_paused.isSelected());
        }
    };

	private Action key_space = new AbstractAction(null) {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1490613691285764886L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_paused.setSelected(!opts.gui_paused.isSelected());
        }
    };
    
    private Action key_f = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -4904966391361916427L;

		@Override
        public void actionPerformed(ActionEvent e) {
			advanceframe = true;
        }
    };
    private Action key_1 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 7068876033770630932L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(0);
        }
    };
    private Action key_2 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 3073920795953268113L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(1);
        }
    };
    private Action key_3 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -3058332124536791250L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(2);
        }
    };
    private Action key_4 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -4234202778202512715L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(3);
        }
    };
    private Action key_5 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 3108303844267015392L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(4);
        }
    };
    private Action key_sh1 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 5658413183986668046L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(1);
        }
    };
    private Action key_sh2 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 5406187940315930665L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(2);
        }
    };
    private Action key_sh3 = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 6811823244960746777L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(3);
        }
    };
    private Action key_c = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 1008228098024062489L;

		@Override
        public void actionPerformed(ActionEvent e) {
			clear = true;
        }
    };
    private Action key_r = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 4617767307508428830L;

		@Override
        public void actionPerformed(ActionEvent e) {
			reset = true;
        }
    };
    private Action key_m = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -8185767523763271616L;

		@Override
        public void actionPerformed(ActionEvent e) {
			//measure = !measure;
        }
    };
    private Action key_eq = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -31743209394080263L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brightness.setValue(opts.gui_brightness.getValue()+1);
        }
    };
    private Action key_minus = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = -1057139941631626988L;

		@Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brightness.setValue(opts.gui_brightness.getValue()-1);
        }
    };
    private Action key_tab = new AbstractAction(null) {
        /**
		 * 
		 */
		private static final long serialVersionUID = 3773665211041914687L;

		@Override
        public void actionPerformed(ActionEvent e) {
    		opts.gui_brush_1.setSelectedIndex((opts.gui_brush_1.getSelectedIndex()+1)%2);
        }
    };
    
	public Action[] keys = {key_p, key_space, key_f, key_1, key_2, key_3, key_4, key_5, key_sh1, key_sh2, key_sh3, key_c, key_r, key_m, key_eq, key_minus, key_tab};
	public KeyStroke[] keystrokes = {
			KeyStroke.getKeyStroke(KeyEvent.VK_P, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_F, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_1, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_2, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_3, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_4, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_5, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_1, InputEvent.SHIFT_DOWN_MASK),
			KeyStroke.getKeyStroke(KeyEvent.VK_2, InputEvent.SHIFT_DOWN_MASK),
			KeyStroke.getKeyStroke(KeyEvent.VK_3, InputEvent.SHIFT_DOWN_MASK),
			KeyStroke.getKeyStroke(KeyEvent.VK_C, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_R, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_M, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_EQUALS, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_MINUS, 0),
			KeyStroke.getKeyStroke(KeyEvent.VK_Q, 0)};
    
	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		if (e.getPreciseWheelRotation() < 0)
			opts.gui_brushsize.setValue(opts.gui_brushsize.getValue()+1);
		if (e.getPreciseWheelRotation() > 0)
			opts.gui_brushsize.setValue(opts.gui_brushsize.getValue()-1);
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
			opts.textPane.setEnabled(!opts.textPane.isEnabled());
			opts.textPane.setEditable(opts.textPane.isEnabled());
		}
	}

}

class RenderCanvas extends JPanel {

	/**
	 * 
	 */
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
	boolean outputenabled = true;
	boolean outputavg = false;
	double avgtime = 0;
	
	public Timer(String name, boolean enabled) {
		this.name = name;
		this.outputenabled = enabled;
	}

	void start() {
		if (outputenabled) {
			tstart = System.nanoTime();
		}
	}

	void disableOutput() {
		outputenabled = false;
	}
	void enableOutput() {
		outputenabled = true;
	}

	void stop() {
		if (outputenabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			double mstime = diff/1000000.0;
			avgtime = avgtime*0.98+mstime*0.02;
			if (outputavg)
				System.out.println(name + ", dt: " + String.format("%.2f", avgtime) + " ms");
			else
				System.out.println(name + ", dt: " + String.format("%.2f", mstime) + " ms");
		}
	}

	void stop(String msg) {
		if (outputenabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			double mstime = diff/1000000.0;
			avgtime = avgtime*0.98+mstime*0.02;
			System.out.println(name + ", dt: " + String.format("%.2f", mstime) + " ms");
			System.out.println("Message: " + msg + "\n");
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

class Material implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 4264510250493981234L;
	double D = 0.0;
	double sigma = 0.0;
	double emfx = 0.0;
	double emfy = 0.0;
	int frequency = 0;
	double epsr = 1.0;
	double mur = 1.0;

	int type = 0;
	boolean activated = true;
	int conductivity_on = 1;
	int emf_on = 1;
	int diffusion_on = 1;
	
	static int VACUUM = 0;
	static int CONDUCTOR = 1;
	static int DIELECTRIC = 2;
	static int EMFSOURCE = 3;
	static int PARAMAGNET = 4;
	static int SWITCH = 5;
	
	static String[] NAMES = {"Vacuum", "Conductor", "Dielectric", "EMF Source", "Paramagnet", "Switch"};
	
	public void reset() {
		type = VACUUM;
		D = 0;
		sigma = 0;
		emfx = 0;
		emfy = 0;
		frequency = 0;
		epsr = 1.0;
		mur = 1.0;
		
		activated = true;
		conductivity_on = 1;
		emf_on = 1;
		diffusion_on = 1;
	}
}
