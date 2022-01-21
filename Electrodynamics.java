import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
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
import java.util.ArrayList;
import java.util.Random;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileFilter;

public class Electrodynamics implements Runnable, MouseListener, MouseMotionListener, MouseWheelListener, ActionListener {
	
	//different aspect ratios, multithreading
	//Units: cm, ns
	
	double[][] display;
	
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
	double[][] conductivityx;
	double[][] conductivityy;
	double[][] emfx;
	double[][] emfy;
	double[][] eps_rx;
	double[][] eps_ry;
	double[][] mu_rz;

	/* Auxiliary fields */
	double[][] dx; //displacement
	double[][] dy;
	double[][] hz; //H field
	double[][] sx; //poynting
	double[][] sy;
	double[][] u; //energy
	
	/* Domain size parameters */
	double xmin = 0;
	double xmax = 50;
	double ymin = 0;
	double ymax = 50;
	double width = xmax-xmin;
	double height = ymax-ymin;
	double ds = 0.2;
	double dt = 0;
	int nx = (int)(width/ds)+1;
	int ny = (int)(height/ds)+1;
	//public static int nx = 257;
	//public static int ny = 257;
	
	int interiorpoints = (nx-2)*(ny-2);
	
	/* Physical constants */
	double eps0 = 1/30.0;
	double mu0 = 1/30.0;
	double c = 1/Math.sqrt(eps0*mu0);
	double csq = 1/(eps0*mu0);
	double D = 2.0;
	double Lambda = 5.0;
	double q = 1.0;
	double sigma = c*(dt/ds);
	double sigmasq = sigma*sigma;
	double tauminus = 1-c*dt/ds;
	double tauplus = 1+c*dt/ds;
	double maxdt = 0.95*ds/(Math.sqrt(2)*c);
	
	double time = 0.0;
	int stepnumber = 0;
	int lastsimspeed = 0;
	double phase = 0;
	
	boolean advanceframe = false;
	boolean clear = false;
	boolean reset = false;
	boolean measure = false;

	/* Graphics */
	RenderCanvas r;
	MainWindow opts;
	BufferedImage screen;
	int scalefactor = 3;
	int imgwidth = 0;
	int imgheight = 0;
	
	/* Timing */
	Timer t1 = new Timer("Poisson 1", false);
	Timer t2 = new Timer("Poisson 2", false);
	Timer t3 = new Timer("Poisson 3 intermediate", false);
	Timer t4 = new Timer("Poisson 3 total", true);
	Timer t5 = new Timer("Graphics", false);
	Timer t6 = new Timer("FD loop", false);
	
	/* Mouse controls */
	int brush_addcharge = 0;
	int brush_addinsulator = 1;
	int brush_addemf = 2;
	int brush_adddielectric = 3;
	int brush_addparamagnet = 4;
	boolean drawingcharge = false;
	boolean linemode = false;
	int lineorientation = 0;
	int mouseX = 0;
	int mouseY = 250;
	int mouseXstart = 0;
	int mouseYstart = 250;
	int mouseButton = 0;
	boolean mouseIsPressed = false;

	public static void main(String[] args) {
		Electrodynamics w = new Electrodynamics();
		w.run();
	}
	
	public Electrodynamics() {
		//System.out.println("CFL condition: ds/sqrt(2)dt = " + (ds/(1.414*dt)) + ", c = " + c);
		System.out.println(maxdt);
		System.out.println(ds*ds/(4*D));
		System.out.println(maxdt*Lambda/eps0);
		opts = new MainWindow();
		display = new double[nx][ny];
		

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
		conductivityx = new double[nx][ny];
		conductivityy = new double[nx][ny];
		emfx = new double[nx][ny];
		emfy = new double[nx][ny];
		eps_rx = new double[nx][ny];
		eps_ry = new double[nx][ny];
		mu_rz = new double[nx][ny];
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
		opts.pack();
		
		for (int i = 0; i < keys.length; i++) {
			Action key = keys[i];
			opts.contentPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
					keystrokes[i], key);
			opts.contentPane.getActionMap().put(key, key);
		}

	}

	@Override
	public void run() {
		while (true) {
			try {
				Thread.sleep(16);
			} catch (InterruptedException e) {}
			
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
				sigma = c*(dt/ds);
				sigmasq = sigma*sigma;
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
			
			r.repaint();
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
				display[i][j]=0.0;
				ex [i][j] = 0.0;
				ey [i][j] = 0.0;

				bz [i][j] = 0.0;
				bzp [i][j] = 0.0;


				jx [i][j] = 0.0;
				jy [i][j] = 0.0;


				poisson1 [i][j] = 0.0;
				poisson2 [i][j] = 0.0;
				
				dx[i][j] = 0.0;
				dy[i][j] = 0.0;
				hz[i][j] = 0.0;
				sx[i][j] = 0.0;
				sy[i][j] = 0.0;
				u[i][j] = 0.0;
				
				if (resetall) {
					//double x = i*ds;
					//double y = j*ds;
					//rho[i][j] = (1+Math.signum(200.0 - (x-25)*(x-25) - (y-25)*(y-25)));
					//double dist = Math.sqrt((x-25)*(x-25) + (y-25)*(y-25));
					//rho[i][j] = 100.0*Math.exp(-dist*dist);
					rho [i][j] = 0.0;
					
					if (materials[i][j] == null)
						materials[i][j] = new Material();
					else
						materials[i][j].reset();
					
					conductivityx[i][j] = 0.0;
					conductivityy[i][j] = 0.0;
					emfx[i][j] = 0.0;
					emfy[i][j] = 0.0;
					//double dist2 = Math.sqrt((x-25)*(x-25) + (y-30)*(y-30));
					eps_rx[i][j] = 1.0;
					eps_ry[i][j] = 1.0;
					mu_rz[i][j] = 1.0;
					
					//epsx[i][j] = 1 + 10.0*sigmoid(-(dist2-2.0));
					//epsy[i][j] = 1 + 10.0*sigmoid(-(dist2-2.0));
				}
			}
		}


		if (resetall) {
			/*for (int i = 1; i < nx-2; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					conductivityx[i][j] = 1.0;
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-2; j++)
				{
					conductivityy[i][j] = 1.0;
				}
			}*/
		}
	}
	
	public double mxp = 0.0;
	public double myp = 0.0;
	
	public void handleMouseInput() {
		int brush = opts.gui_brush.getSelectedIndex();
		int brushshape = opts.gui_brush_1.getSelectedIndex();
		if (mouseIsPressed) {
			if (!drawingcharge) {
				drawingcharge = true;
				opts.gui_paused.setSelected(true);
			}
		} else {
			if (drawingcharge) {
				drawingcharge = false;
				prescaleDielectric();
				correctEfield();
			}
		}

		int directionval = opts.gui_direction.getValue();
		opts.gui_label_direction.setText("Brush direction: " + directionval*(360/24) + " deg");
		double angle = Math.PI * directionval/12.0;

		double mx = (mouseX/scalefactor)*ds;
		double my = (mouseY/scalefactor)*ds;
		
		if (mouseIsPressed) {
			
			double brushsize = 1.5*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
			double intensity = 1.0*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
			
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
					if (brushshape == 0)
						r = Math.sqrt(p.dot(p));
					else
						r = Math.max(Math.abs(p.x), Math.abs(p.y));
					if (r <= brushsize) {
						if (brush == brush_addcharge) {
							double drho = 0.0;
							if (mouseButton == MouseEvent.BUTTON1)
								drho = intensity;
							else if (mouseButton == MouseEvent.BUTTON3)
								drho = -intensity;
							rho[i][j] += drho*Math.exp(-3.0*r/brushsize);
						} else {
							if (brush == brush_addinsulator) {
								if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.m_vacuum) {
									materials[i][j].reset();
									materials[i][j].type = Material.m_conductor;
									materials[i][j].conductivity = 1.0;
								}
							} else if (brush == brush_addemf) {
								if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.m_vacuum) {
									materials[i][j].reset();
									materials[i][j].type = Material.m_battery;
									materials[i][j].emfx = 10.0*intensity*Math.cos(angle);
									materials[i][j].emfy = 10.0*intensity*Math.sin(-angle);
									materials[i][j].conductivity = 1.0;
								}
							} else if (brush == brush_adddielectric) {
								if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.m_vacuum) {
									materials[i][j].reset();
									materials[i][j].type = Material.m_dielectric;
									materials[i][j].epsr = 1+intensity;
								}
							} else if (brush == brush_addparamagnet) {
								if (mouseButton == MouseEvent.BUTTON1 && materials[i][j].type == Material.m_vacuum) {
									materials[i][j].reset();
									materials[i][j].type = Material.m_paramagnet;
									materials[i][j].mur = 1+intensity;
								}
							}
							
							if (mouseButton == MouseEvent.BUTTON3) {
								materials[i][j].reset();
							}
							
							updateMaterial(0, i, j);
							updateMaterial(0, i, j+1);
							updateMaterial(1, i, j);
							updateMaterial(1, i+1, j);
							mu_rz[i][j] = materials[i][j].mur;
						}
					}
				}
			}

		}
		mxp = mx;
		myp = my;
		//u[mouseX/RenderCanvas.INTERVAL][mouseY/RenderCanvas.INTERVAL] = 10.0;//10.0*Math.sin(time*20.0);
			//10.0*Math.sin(time*20.0);
	}
	
	public void updateMaterial(int orientation, int i, int j) {
		if (orientation == 0) {
			int i1 = i;
			int j1 = j-1;
			int i2 = i;
			int j2 = j;
			
			if (j1 < 0)
				j1 = 0;
			if (j2 == ny-1)
				j2 = ny-2;
			
			conductivityx[i][j] = Math.max(materials[i1][j1].conductivity, materials[i2][j2].conductivity);
			emfx[i][j] = Math.max(materials[i1][j1].emfx, materials[i2][j2].emfx);
			eps_rx[i][j] = Math.max(materials[i1][j1].epsr, materials[i2][j2].epsr);
		} else if (orientation == 1) {
			int i1 = i;
			int j1 = j;
			int i2 = i-1;
			int j2 = j;
			
			if (i1 < 0)
				i1 = 0;
			if (i2 == nx-1)
				i2 = nx-2;
			
			conductivityy[i][j] = Math.max(materials[i1][j1].conductivity, materials[i2][j2].conductivity);
			emfy[i][j] = Math.max(materials[i1][j1].emfy, materials[i2][j2].emfy);
			eps_ry[i][j] = Math.max(materials[i1][j1].epsr, materials[i2][j2].epsr);
		}
	}

	public void updateMainFields() {

		t6.start();
		
		stepnumber++;

		/*Interior*/
		
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
		//double sqrt2 = Math.sqrt(2);
		//bz[0][0] = ((sqrt2*dt/ds)*(-bz[1][1] - bzp[1][1] + bzp[0][0]) + bzp[0][0]+ bzp[0][1]+ bzp[1][0]+ bzp[1][1] - bz[0][1] - bz[1][0] - bz[1][1])/(1-sqrt2*dt/ds);
		bz[0][0] = 0.5*(bz[1][0] + bz[0][1]);
		bz[nx-2][0] = 0.5*(bz[nx-3][0] + bz[nx-2][1]);
		bz[0][ny-2] = 0.5*(bz[0][ny-3] + bz[1][ny-2]);
		bz[nx-2][ny-2] = 0.5*(bz[nx-3][ny-2] + bz[nx-2][ny-3]);

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

		/*E field*/
		
		
		phase += 10.0*Math.pow(10.0, opts.gui_emf_freq.getValue()/15.0) * dt;
		double emffactor = 1.0;
		if (opts.gui_emftype.getSelectedIndex() == 1) {
			emffactor = Math.cos(phase);
		}
		if (!opts.gui_emf.isSelected())
			emffactor = 0.0;

		for (int i = 0; i < nx-1; i++)
		{
			for (int j = 1; j < ny-1; j++)
			{
				jx[i][j] = (-D*(rho[i+1][j] - rho[i][j])/ds + emfx[i][j]*emffactor + 0.5*Lambda*ex[i][j])*conductivityx[i][j];
				ex[i][j] = (ex[i][j] + (csq*(bz[i][j]/mu_rz[i][j]-bz[i][j-1]/mu_rz[i][j-1])/ds - jx[i][j]/eps0)*(dt/eps_rx[i][j]))
						/(1+0.5*dt*Lambda/(eps0*eps_rx[i][j])*conductivityx[i][j]);
				jx[i][j] += 0.5*ex[i][j]*Lambda*conductivityx[i][j];
			}
		}

		for (int i = 1; i < nx-1; i++)
		{
			for (int j = 0; j < ny-1; j++)
			{
				jy[i][j] = (-D*(rho[i][j+1] - rho[i][j])/ds + emfy[i][j]*emffactor + 0.5*Lambda*ey[i][j])*conductivityy[i][j];
				ey[i][j] = (ey[i][j] + (-csq*(bz[i][j]/mu_rz[i][j]-bz[i-1][j]/mu_rz[i-1][j])/ds - jy[i][j]/eps0)*(dt/eps_ry[i][j]))
						/(1+0.5*dt*Lambda/(eps0*eps_ry[i][j])*conductivityy[i][j]);
				jy[i][j] += 0.5*ey[i][j]*Lambda*conductivityy[i][j];
			}
		}

		if (stepnumber%40 == 0) {
			correctEfield();
		}

		time += dt;
		advanceframe = false;

		t6.stop();
	}

	public void calcAuxillaryFields() {
		if (opts.gui_view.getSelectedIndex() == 4) {
			for (int i = 1; i < nx-2; i++)
			{
				for (int j = 1; j < ny-2; j++)
				{
					hz[i][j] = bz[i][j]/mu_rz[i][j];
				}
			}
		}

		else if (opts.gui_view.getSelectedIndex() == 5) {
			for (int i = 1; i < nx-2; i++)
			{
				for (int j = 1; j < ny-2; j++)
				{
					u[i][j] = (ex[i][j]*ex[i][j]*eps_rx[i][j]
							+ ex[i][j+1]*ex[i][j+1]*eps_rx[i][j+1]
									+ ey[i][j]*ey[i][j]*eps_ry[i][j]
											+ey[i+1][j]*ey[i+1][j]*eps_ry[i+1][j])/8
							+ bz[i][j]*bz[i][j]/(2*mu_rz[i][j]);
				}
			}
		}

		if (opts.gui_view_vec.getSelectedIndex() == 2) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					dx[i][j] = eps_rx[i][j]*ex[i][j];
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					dy[i][j] = eps_ry[i][j]*ey[i][j];
				}
			}
		}

		else if (opts.gui_view_vec.getSelectedIndex() == 3) {
			for (int i = 0; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					sy[i][j] = -ex[i][j]*0.5*(bz[i][j]/mu_rz[i][j] + bz[i][j-1]/mu_rz[i][j-1]);
				}
			}

			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 0; j < ny-1; j++)
				{
					sx[i][j] = ey[i][j]*0.5*(bz[i][j]/mu_rz[i][j] + bz[i-1][j]/mu_rz[i-1][j]);
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
			downscale(eps_rx, epsxMG[fineness], nx-2, ny-1, nxcgint, nycgint+1, ds/2.0, 0, ds, gridsize/2.0, 0, gridsize);
			downscale(eps_ry, epsyMG[fineness], nx-1, ny-2, nxcgint+1, nycgint, 0, ds/2.0, ds, 0, gridsize/2.0, gridsize);
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
				deltarho[i][j] = eps0*((ex[i][j]*eps_rx[i][j]-ex[i-1][j]*eps_rx[i-1][j] + ey[i][j]*eps_ry[i][j]-ey[i][j-1]*eps_ry[i][j-1])/ds) - rho[i][j];
			}
		}
		
		double err = 0;
		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				err += Math.abs(((eps_rx[i-1][j]+eps_rx[i][j]+eps_ry[i][j-1]+eps_ry[i][j])*poisson1[i][j]
						- (poisson1[i-1][j]*eps_rx[i-1][j] + poisson1[i+1][j]*eps_rx[i][j] + poisson1[i][j-1]*eps_ry[i][j-1] + poisson1[i][j+1]*eps_ry[i][j]))
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
				epsxMG[maxfineness+1][i][j] = eps_rx[i][j];
				epsyMG[maxfineness+1][i][j] = eps_ry[i][j];
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
				err += Math.abs(((eps_rx[i-1][j]+eps_rx[i][j]+eps_ry[i][j-1]+eps_ry[i][j])*poisson1[i][j]
						- (poisson1[i-1][j]*eps_rx[i-1][j] + poisson1[i+1][j]*eps_rx[i][j] + poisson1[i][j-1]*eps_ry[i][j-1] + poisson1[i][j+1]*eps_ry[i][j]))
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
		double scalingconstant = 0.5*Math.pow(10.0, opts.gui_brightness.getValue()/10.0);
		int scalarview = opts.gui_view.getSelectedIndex();
		int vectorview = opts.gui_view_vec.getSelectedIndex();
		setalphaBG(0);
		setalphaFG(1);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				int x = i;
				int y = j;

				if (scalarview == 0) {
					if (i < nx-1 && j < ny-1) {
						double exg = 0.5*(ex[i][j]+ex[i][j+1])*scalingconstant;
						double eyg = 0.5*(ey[i][j]+ey[i+1][j])*scalingconstant;
						setcol(exg, 0, eyg, 30);
					}
				}
				else if (scalarview == 1) {
					double bzg = this.bz[i][j]*c*scalingconstant;
					setcol(clamp(bzg, 0, 999), 0, clamp(-bzg, 0, 999), 30);
				}
				else if (scalarview == 2) {
					x -= 0.5;
					y -= 0.5;
					double rhog = rho[i][j]*5*ds*scalingconstant;
					setcol(clamp(rhog, 0, 999), 0, clamp(-rhog, 0, 999), 30);
				}
				else if (scalarview == 3) {
					if (i < nx-1 && j < ny-1) {
						double jxg = 0.5*(jx[i][j]+jx[i][j+1])*scalingconstant;
						double jyg = 0.5*(jy[i][j]+jy[i+1][j])*scalingconstant;
						setcol(jxg, 0, jyg, 30);
					}
				}
				else if (scalarview == 4) {
					double hzg = this.hz[i][j]*c*scalingconstant;
					setcol(clamp(hzg, 0, 999), 0, clamp(-hzg, 0, 999), 30);
				}
				else if (scalarview == 5) {
					double ug = this.u[i][j]*scalingconstant;
					setcol(0, ug, 0, 30);
				}
				
				fillRectangle((int)(x*scalefactor), (int)(y*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
			}
		}

		setalphaBG(1);
		setalphaFG(0.5);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (materials[i][j].conductivity > 0.5) {
					setColor(180, 180, 180);
					fillRectangle((int)(i*scalefactor), (int)(j*scalefactor), scalefactor, scalefactor);
				}
				if (materials[i][j].epsr > 1) {
					setColor(90, 50, 0);
					fillRectangle((int)(i*scalefactor), (int)(j*scalefactor), scalefactor, scalefactor);
				}
				if (materials[i][j].mur > 1) {
					setColor(77, 120, 100);
					fillRectangle((int)(i*scalefactor), (int)(j*scalefactor), scalefactor, scalefactor);
				}
				if (materials[i][j].emfx*materials[i][j].emfx + materials[i][j].emfy*materials[i][j].emfy > 0) {
					setColor(180, 0, 180);
					fillRectangle((int)(i*scalefactor), (int)(j*scalefactor), scalefactor, scalefactor);
				}
			}
		}

		setalphaBG(1);
		setalphaFG(0.3);
		
		double arrowlength = scalefactor;

		double vectorscalingconstant = 5.0*Math.pow(10.0, opts.gui_brightness_1.getValue()/10.0);
		for (int i = 0; i < 50; i++) {
			for (int j = 0; j < 50; j++) {
				double x = (nx-1)*(i+0.5)/50;
				double y = (ny-1)*(j+0.5)/50;
				Vector ctr = new Vector(x, y);
				Vector arrow = null;
				
				if (vectorview == 0)
					arrow = new Vector(getex(x, y), getey(x, y));
				else if (vectorview == 1)
					arrow = new Vector(getjx(x, y)/Lambda, getjy(x, y)/Lambda);
				else if (vectorview == 2)
					arrow = new Vector(getdx(x, y), getdy(x, y));
				else if (vectorview == 3)
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

		
		g.setColor(Color.WHITE);

		((Graphics2D)g).setRenderingHint(
		        RenderingHints.KEY_TEXT_ANTIALIASING,
		        RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		double mx = (mouseX/scalefactor);
		double my = (mouseY/scalefactor);
		int mi = (int)Math.floor(mx);
		int mj = (int)Math.floor(my);
		Material mat = materials[mi][mj];
		
		int vspacing = 12;
		int voffset = 1;
		int hoffset = 5;
		double eg = length(getex(mx,my), getey(mx,my));
		double bg = getbz(mx, my);
		double rhog = getrho(mx, my);
		double jg = length(getjx(mx,my), getjy(mx,my));
		g.drawString("E: " + getSI(eg, "V/m"), hoffset, voffset + 1*vspacing);
		g.drawString("B: " + getSI(bg, "T"), hoffset, voffset + 2*vspacing);
		g.drawString("\u03c1: " + getSI(rhog, "C/m^3"), hoffset, voffset + 3*vspacing);
		g.drawString("J: " + getSI(jg, "A/m^2"), hoffset, voffset + 4*vspacing);
		g.drawString("\u03c3: " + getSI(mat.conductivity, "S"), hoffset, voffset + 5*vspacing);
		g.drawString("emf: " + getSI(length(mat.emfx, mat.emfy), "V/m"), hoffset, voffset + 6*vspacing);
		g.drawString("\u03b5r: " + getSI(mat.epsr, ""), hoffset, voffset + 7*vspacing);
		g.drawString("\u03bcr: " + getSI(mat.mur, ""), hoffset, voffset + 8*vspacing);
		
		double brushsize = 1.5*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
		int r = (int)(scalefactor*brushsize/ds);
		int brushshape = opts.gui_brush_1.getSelectedIndex();
		if (brushshape == 0) {
		g.setColor(new Color(50, 50, 50));
		g.drawOval(mouseX - r+1, mouseY - r+1, 2*r-2, 2*r-2);
		g.setColor(new Color(200, 200, 200));
		g.drawOval(mouseX - r, mouseY - r, 2*r, 2*r);
		} else {
			g.setColor(new Color(50, 50, 50));
			g.drawRect(mouseX - r+1, mouseY - r+1, 2*r-2, 2*r-2);
			g.setColor(new Color(200, 200, 200));
			g.drawRect(mouseX - r, mouseY - r, 2*r, 2*r);
		}
		t5.stop();
		//g.setColor(.)
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
	
	public boolean load()
	{

		JFileChooser fd = new JFileChooser();
		fd.setFileFilter(new FileFilter(){
			public boolean accept(File f) {
				if (f.isDirectory()) {
					return true;
				}
				return true;
				//if (f.getName().endsWith(".txt")) return true;
				//return false;
			}
			@Override
			public String getDescription() {
				return ".txt";
			}
		});
		fd.setVisible(true);
		fd.showOpenDialog(opts);
		File infile = fd.getSelectedFile();
		if (infile == null) return false;
		try {
			FileInputStream fstr = new FileInputStream(infile);
			ObjectInputStream ostr = new ObjectInputStream(fstr);

			time = ostr.readDouble();
			phase = ostr.readDouble();
			opts.gui_paused.setSelected(ostr.readBoolean());
			opts.gui_view.setSelectedIndex(ostr.readInt());
			opts.gui_view_vec.setSelectedIndex(ostr.readInt());
			opts.gui_emftype.setSelectedIndex(ostr.readInt());
			opts.gui_simspeed.setValue(ostr.readInt());
			opts.gui_simspeed_2.setValue(ostr.readInt());
			opts.gui_brightness.setValue(ostr.readInt());
			opts.gui_brightness_1.setValue(ostr.readInt());
			opts.gui_emf.setSelected(ostr.readBoolean());
			opts.gui_emf_freq.setValue(ostr.readInt());
			opts.gui_brush.setSelectedIndex(ostr.readInt());
			opts.gui_brush_1.setSelectedIndex(ostr.readInt());
			opts.gui_brushsize.setValue(ostr.readInt());
			opts.gui_brushintensity2.setValue(ostr.readInt());
			opts.gui_direction.setValue(ostr.readInt());

			(ex) = (double[][]) ostr.readObject();
			(ey) = (double[][]) ostr.readObject();
			(bz) = (double[][]) ostr.readObject();
			(rho) = (double[][]) ostr.readObject();
			(jx) = (double[][]) ostr.readObject();
			(jy) = (double[][]) ostr.readObject();
			(bzp) = (double[][]) ostr.readObject();
			(materials) = (Material[][]) ostr.readObject();
			(conductivityx) = (double[][]) ostr.readObject();
			(conductivityy) = (double[][]) ostr.readObject();
			(emfx) = (double[][]) ostr.readObject();
			(emfy) = (double[][]) ostr.readObject();
			(eps_rx) = (double[][]) ostr.readObject();
			(eps_ry) = (double[][]) ostr.readObject();
			(mu_rz) = (double[][]) ostr.readObject();
			
			prescaleDielectric();
			
			ostr.close();
			fstr.close();
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException e) {
			return false;
		} catch (ClassNotFoundException e) {
			return false;
		}
		return true;
	}

	public boolean save()
	{
		JFileChooser fd = new JFileChooser();
		fd.setFileFilter(new FileFilter(){
			public boolean accept(File f) {
				if (f.isDirectory()) {
					return true;
				}
				if (f.getName().endsWith(".txt")) return true;
				return false;
			}
			@Override
			public String getDescription() {
				return ".txt";
			}
		});
		fd.showSaveDialog(opts);
		File outfile = fd.getSelectedFile();
		if (outfile == null) return false;
		try {
			FileOutputStream fstr = new FileOutputStream(outfile);
			ObjectOutputStream ostr = new ObjectOutputStream(fstr);
			
			ostr.writeDouble(time);
			ostr.writeDouble(phase);
			
			ostr.writeBoolean(opts.gui_paused.isSelected());
			ostr.writeInt(opts.gui_view.getSelectedIndex());
			ostr.writeInt(opts.gui_view_vec.getSelectedIndex());
			ostr.writeInt(opts.gui_emftype.getSelectedIndex());
			ostr.writeInt(opts.gui_simspeed.getValue());
			ostr.writeInt(opts.gui_simspeed_2.getValue());
			ostr.writeInt(opts.gui_brightness.getValue());
			ostr.writeInt(opts.gui_brightness_1.getValue());
			ostr.writeBoolean(opts.gui_emf.isSelected());
			ostr.writeInt(opts.gui_emf_freq.getValue());
			ostr.writeInt(opts.gui_brush.getSelectedIndex());
			ostr.writeInt(opts.gui_brush_1.getSelectedIndex());
			ostr.writeInt(opts.gui_brushsize.getValue());
			ostr.writeInt(opts.gui_brushintensity2.getValue());
			ostr.writeInt(opts.gui_direction.getValue());
			
			ostr.writeObject(ex);
			ostr.writeObject(ey);
			ostr.writeObject(bz);
			ostr.writeObject(rho);
			ostr.writeObject(jx);
			ostr.writeObject(jy);
			ostr.writeObject(bzp);
			ostr.writeObject(materials);
			ostr.writeObject(conductivityx);
			ostr.writeObject(conductivityy);
			ostr.writeObject(emfx);
			ostr.writeObject(emfy);
			ostr.writeObject(eps_rx);
			ostr.writeObject(eps_ry);
			ostr.writeObject(mu_rz);
			
			
			ostr.close();
			fstr.close();
		} catch (FileNotFoundException e) {
			return false;
		} catch (IOException e) {
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
		if (e.isShiftDown()) {
			linemode = true;
			mxp = (mouseX/scalefactor)*ds;
			myp = (mouseY/scalefactor)*ds;
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
				if (deltax > 10) {
					lineorientation = 1;
				} else if (deltax < -10) {
					lineorientation = 1;
				} else if (deltay > 10) {
					lineorientation = 2;
				} else if (deltay < -10) {
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
			measure = !measure;
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
    
	public Action[] keys = {key_p, key_f, key_1, key_2, key_3, key_4, key_5, key_sh1, key_sh2, key_sh3, key_c, key_r, key_m, key_eq, key_minus, key_tab};
	public KeyStroke[] keystrokes = {
			KeyStroke.getKeyStroke(KeyEvent.VK_P, 0),
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
			save();
		else if (e.getSource() == opts.gui_open)
			load();
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
			System.out.println(name + ", dt: " + String.format("%.2f", diff/1000000.0) + " ms");
		}
	}

	void stop(String msg) {
		if (outputenabled) {
			long tend = System.nanoTime();
			long diff = tend - tstart;
			System.out.println(name + ", dt: " + String.format("%.2f", diff/1000000.0) + " ms");
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
	double conductivity = 0.0;
	double emfx = 0.0;
	double emfy = 0.0;
	double epsr = 1.0;
	double mur = 1.0;

	static int m_vacuum = 0;
	static int m_conductor = 1;
	static int m_dielectric = 2;
	static int m_battery = 3;
	static int m_paramagnet = 4;
	
	int type = 0;
	
	public void reset() {
		type = 0;
		conductivity = 0;
		emfx = 0;
		emfy = 0;
		epsr = 1.0;
		mur = 1.0;
	}
}
