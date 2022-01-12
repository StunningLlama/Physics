/*
 * Brandon Li 
 * Microwave simulation for 2210
 * 11/14/2021
 */

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
import java.util.ArrayList;
import java.util.Random;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.KeyStroke;

public class Electrodynamics implements Runnable, MouseListener, MouseMotionListener, MouseWheelListener, ActionListener {
	
	//Units: cm, ns
	
	double[][] display;
	
	double[][] ex;
	double[][] ey;
	double[][] bz;
	double[][] rho;
	double[][] jx;
	double[][] jy;

	double[][] bzp;

	double[][] deltarho;
	double[][] deltarhoMG;
	double[][] poisson1;
	double[][] poisson2;
	
	double[][] conductivityx;
	double[][] conductivityy;
	double[][] emfx;
	double[][] emfy;

	int[][] indices;
	int[] xindices;
	int[] yindices;
	//AbstractMatrix delsq;
	//AbstractVector b;
	//AbstractVector poissonsol;
	
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
	
	double eps0 = 1/30.0;
	double mu0 = 1/30.0;
	double c = 1/Math.sqrt(eps0*mu0);
	double csq = 1/(eps0*mu0);
	double D = 2.0;
	double Lambda = 2.0;
	double q = 1.0;
	double sigma = c*(dt/ds);
	double sigmasq = sigma*sigma;
	double tauminus = 1-c*dt/ds;
	double tauplus = 1+c*dt/ds;
	double maxdt = 0.95*ds/(Math.sqrt(2)*c);
	
	double obstacleradius = 1.4;
	
	double time = 0.0;
	int iterationmultiplier = 1;
	
	//boolean paused = true;
	boolean advanceframe = false;
	boolean clear = false;
	boolean reset = false;
	boolean measure = false;
	//int brightness = 0;
	//int view = 1;
	RenderCanvas r;
	MainWindow opts;
	Image screen;
	double scalefactor = ds*15.0;

	int imgwidth = 0;
	int imgheight = 0;
	
	Timer t1 = new Timer("Poisson 1", false);
	Timer t2 = new Timer("Poisson 2", false);
	Timer t3 = new Timer("Poisson 3 intermediate", false);
	Timer t4 = new Timer("Poisson 3 total", false);
	Timer t5 = new Timer("Graphics", false);
	
	boolean experimentended = false;
	double experimentduration = 99999.0;
	
	int brush_addcharge = 0;
	int brush_addinsulator = 1;
	int brush_addemf = 2;
	//int brush = brush_addcharge;
	boolean drawingcharge = false;
	//double brushsize = 0.5;
	int lastsimspeed = 0;

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
		
		conductivityx = new double[nx][ny];
		conductivityy = new double[nx][ny];
		emfx = new double[nx][ny];
		emfy = new double[nx][ny];
	
		resetfields(true);
		r= new RenderCanvas(this);

		opts.setVisible(true);
		imgwidth = (int)Math.ceil(scalefactor*nx);
		imgheight = (int)Math.ceil(scalefactor*ny);
		screen = opts.createImage(imgwidth, imgheight);
		
		
		opts.add(r, BorderLayout.CENTER);
		r.addMouseListener(this);
		r.addMouseMotionListener(this);
		r.addMouseWheelListener(this);
		r.setPreferredSize(new Dimension(imgwidth-1, imgheight));
		opts.gui_reset.addActionListener(this);
		opts.gui_resetall.addActionListener(this);
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
		int frameindex = 0;
		
		while (true) {
			frameindex++;
			try {
				Thread.sleep(17);
			} catch (InterruptedException e) {}
			for (int i = 0; i < iterationmultiplier ; i++) {
				this.Physics();
				this.Physics();
			}
			if (frameindex%2==0)
			{
				r.repaint();
			}
		}
	}
	
	public void resetfields(boolean resetall) {
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				double x = i*ds;
				double y = j*ds;
				
				
				display[i][j]=0.0;
				ex [i][j] = 0.0;
				ey [i][j] = 0.0;

				bz [i][j] = 0.0;
				bzp [i][j] = 0.0;

				//rho[i][j] = (1+Math.signum(200.0 - (x-25)*(x-25) - (y-25)*(y-25)));
				//rho[i][j] = 100.0*(1+Math.signum(1.0 - (x-25)*(x-25) - (y-25)*(y-25)));

				jx [i][j] = 0.0;
				jy [i][j] = 0.0;


				poisson1 [i][j] = 0.0;
				poisson2 [i][j] = 0.0;
				
				if (resetall) {
					rho [i][j] = 0.0;
					conductivityx[i][j] = 0.0;
					conductivityy[i][j] = 0.0;
					emfx[i][j] = 0.0;
					emfy[i][j] = 0.0;
				}
			}
		}


		for (int i = 1; i < nx-2; i++)
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
		}

		correctEfield();
	}
	
	public double mxp = 0.0;
	public double myp = 0.0;
	
	public void handleMouseInput() {
		int brush = opts.gui_brush.getSelectedIndex();
		if (mouseIsPressed && brush == brush_addcharge) {
			if (!drawingcharge) {
				drawingcharge = true;
				opts.gui_paused.setSelected(true);
			}
		} else {
			if (drawingcharge) {
				drawingcharge = false;
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
					if (brush == brush_addcharge || brush == brush_addinsulator) {
						 cx = i*ds;
						 cy = j*ds;
					} else if (brush == brush_addemf) {
						 cx = (i+0.5)*ds;
						 cy = (j+0.5)*ds;
					}
					
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
					double r = Math.sqrt(p.dot(p));
					if (r <= brushsize) {
						if (brush == brush_addinsulator) {
							double newconductivity = 0.0;
							if (mouseButton == MouseEvent.BUTTON3)
								newconductivity = 1;
							conductivityx[i][j] = newconductivity;
							conductivityy[i][j] = newconductivity;
							conductivityx[i-1][j] = newconductivity;
							conductivityy[i][j-1] = newconductivity;
						} else if (brush == brush_addcharge) {
							double drho = 0.0;
							if (mouseButton == MouseEvent.BUTTON1)
								drho = intensity;
							else if (mouseButton == MouseEvent.BUTTON3)
								drho = -intensity;
							rho[i][j] += drho*Math.exp(-3.0*r/brushsize);
						} else if (brush == brush_addemf) {
							if (mouseButton == MouseEvent.BUTTON1) {
								emfx[i][j] = 10.0*intensity*Math.cos(angle);
								emfx[i][j+1] = 10.0*intensity*Math.cos(angle);
								emfy[i][j] = 10.0*intensity*Math.sin(-angle);
								emfy[i+1][j] = 10.0*intensity*Math.sin(-angle);
							} else if (mouseButton == MouseEvent.BUTTON3) {
								emfx[i][j] = 0;
								emfx[i][j+1] = 0;
								emfy[i][j] = 0;
								emfy[i+1][j] = 0;
							}
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
	
	int framenumber = 0;
	public void Physics() {
		
		if (clear || reset) {
			resetfields(reset);
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
			
			if (!experimentended && time > experimentduration) {
				measure = true;
				opts.gui_paused.setSelected(true);
				experimentended = true;
			}
			
			framenumber++;
			
			/*Interior*/

			if (framenumber%2==0) {

				//for (int j = 1; j < ny-1; j++) {
				//	bz[0][j] = 2.0*exp[0][j]-expp[0][j]+sigma*(exp[0+1][j]-exp[0][j]-expp[0+1][j]+expp[0][j]) + (sigmasq/2.0)*(exp[0][j+1]-2*exp[0][j]+exp[0][j-1]);
				//}
				
				for (int i = 1; i < nx-2; i++)
				{
					for (int j = 1; j < ny-2; j++)
					{
							bz[i][j] = bz[i][j] + (-(ey[i+1][j] - ey[i][j]) + (ex[i][j+1] - ex[i][j]))*dt/ds;
							//jx[i][j] = -D*(rho[i+1][j] - rho[i][j])/ds + ex[i][j]*Lambda*0.5*Math.abs(rho[i][j] + rho[i+1][j])/q;
							//jy[i][j] = -D*(rho[i][j+1] - rho[i][j])/ds + ey[i][j]*Lambda*0.5*Math.abs(rho[i][j] + rho[i][j+1])/q;
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
				
			} else {


				
				for (int i = 0; i < nx-1; i++)
				{
					for (int j = 1; j < ny-1; j++)
					{
						jx[i][j] = (-D*(rho[i+1][j] - rho[i][j])/ds +  + emfx[i][j])*conductivityx[i][j];
						ex[i][j] = ex[i][j] + csq*(bz[i][j]-bz[i][j-1])*dt/ds - jx[i][j]/eps0*dt;
						
						jx[i][j] += ex[i][j]*Lambda;
						//ex[i][j] = (ex[i][j] + csq*(bz[i][j]-bz[i][j-1])*dt/ds)/(1+dt*Lambda/eps0*conductivityx[i][j]);
						//jx[i][j] = dt*Lambda/eps0*conductivityx[i][j]*ex[i][j];
					}
				}

				for (int i = 1; i < nx-1; i++)
				{
					for (int j = 0; j < ny-1; j++)
					{
						jy[i][j] = (-D*(rho[i][j+1] - rho[i][j])/ds + ey[i][j]*Lambda + emfy[i][j])*conductivityy[i][j];
						ey[i][j] = ey[i][j] + -csq*(bz[i][j]-bz[i-1][j])*dt/ds - jy[i][j]/eps0*dt;
						//ey[i][j] = (ey[i][j] + -csq*(bz[i][j]-bz[i-1][j])*dt/ds)/(1+dt*Lambda/eps0*conductivityx[i][j]);
						//jy[i][j] = dt*Lambda/eps0*conductivityy[i][j]*ey[i][j];
					}
				}
			}
			
			if (framenumber%40 == 0) {
				correctEfield();
			}
			
			 bz[100][100] = Math.sin(time*100.0);
			
			/* Boundary */
			time += dt/2.0;
			advanceframe = false;
		}
	}
	
	public void correctEfield() {
		boolean debug = false;
		t4.start();
		t3.disableOutput();
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
				deltarho[i][j] = eps0*((ex[i][j]-ex[i-1][j] + ey[i][j]-ey[i][j-1])/ds) - rho[i][j];
			}
		}
		

		double err = 0;
		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				err += Math.abs((4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(ds*ds) - deltarho[i][j]/eps0);
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
			for (int i = 1; i < nxcgint+1; i++) {
				for (int j = 1; j < nycgint+1; j++) {
					deltarhoMG[i][j] = 0;
					double x1 = (i - 0.5)*(nx/(double)(1 << fineness));
					double x2 = (i + 0.5)*(nx/(double)(1 << fineness));
					double y1 = (j - 0.5)*(nx/(double)(1 << fineness));
					double y2 = (j + 0.5)*(nx/(double)(1 << fineness));
					
					int x1f = (int)Math.floor(x1);
					int x2f = (int)Math.floor(x2);
					int y1f = (int)Math.floor(y1);
					int y2f = (int)Math.floor(y2);
					
					double totalarea = 0.0;
					
					for (int i2 = x1f; i2 <= x2f; i2++) {
						for (int j2 = y1f; j2 <= y2f; j2++) {
							double x1r = Math.max(x1, i2);
							double x2r = Math.min(x2, i2+1.0);
							double y1r = Math.max(y1, j2);
							double y2r = Math.min(y2, j2+1.0);
							if (x2r > x1r && y2r > y1r) {
								//System.out.println(i2 + ", " + j2);
								deltarhoMG[i][j] += deltarho[i2][j2]*(x2r-x1r)*(y2r-y1r);
								totalarea += (x2r-x1r)*(y2r-y1r);
							}
						}
					}
					deltarhoMG[i][j] /= totalarea;
				}
			}
			
			int poissonsteps = stepsarray[fineness];
			t3.stop("Charge computation");

			double alpha = (gridsize*gridsize)/eps0;
			double beta = 4;

			t3.start();
			for (int poissonit = 0; poissonit < poissonsteps; poissonit++) {
				for (int i = 1; i < nxcgint+1; i++) {
					for (int j = 1; j < nycgint+1; j++) {
						poisson2[i][j] = ((poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]) + deltarhoMG[i][j]*alpha)/beta;
					}
				}
				for (int i = 1; i < nxcgint+1; i++) {
					for (int j = 1; j < nycgint+1; j++) {
						poisson1[i][j] = ((poisson2[i-1][j] + poisson2[i+1][j] + poisson2[i][j-1] + poisson2[i][j+1]) + deltarhoMG[i][j]*alpha)/beta;
					}
				}

				if (debug) {
					for (int i = 1; i < nxcgint+1; i++) {
						for (int j = 1; j < nycgint+1; j++) {
							bz[i][j] = (4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(gridsize*gridsize) - deltarhoMG[i][j]/eps0;
						}
					}
					for (int i = 1; i < nx-1; i++)
					{
						for (int j = 1; j < ny-1; j++)
						{
							//ex[i][j] = (poisson1[i+1][j]-poisson1[i][j])/(ds*gridsize);
							//ey[i][j] = (poisson1[i][j+1]-poisson1[i][j])/(ds*gridsize);
						}
					}

					r.repaint();
					try {
						Thread.sleep(17);
					} catch (InterruptedException e) {}
				}
			}
			t3.stop("Iteration");


			//for (int i = nxcgint-1; i >= 0; i--)
			//{
			//	for (int j = nycgint-1; j >= 0; j--) {
			//		poisson1[i+1][j+1] = poisson1[i/2 + 1][j/2 + 1];
			//	}
			//}
			

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

		t3.start();
		for (int poissonit = 0; poissonit < poissonsteps; poissonit++) {
			for (int i = 1; i < nxint+1; i++) {
				for (int j = 1; j < nyint+1; j++) {
					poisson2[i][j] = ((poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]) + deltarho[i][j]*alpha)/beta;
				}
			}
			for (int i = 1; i < nxint+1; i++) {
				for (int j = 1; j < nyint+1; j++) {
					poisson1[i][j] = ((poisson2[i-1][j] + poisson2[i+1][j] + poisson2[i][j-1] + poisson2[i][j+1]) + deltarho[i][j]*alpha)/beta;
				}
			}

			if (debug)
			{
				for (int i = 1; i < nx-1; i++) {
					for (int j = 1; j < ny-1; j++) {
						bz[i][j] = (4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(ds*ds) - deltarho[i][j]/eps0;
					}
				}
				r.repaint();
				try {
					Thread.sleep(50);
				} catch (InterruptedException e) {}
			}
		}
		
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
				err += Math.abs((4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(ds*ds) - deltarho[i][j]/eps0);
				//bz[i][j] = (4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(ds*ds) - deltarho[i][j]/eps0;
			}
		}
		System.out.println("Poisson error: " + err);
		t4.stop();
		//t3.stop();
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
		} else if (yfloor < 0) {
			xfloor = 0;
			fx = 0.0;
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
	
	public void render() {
		t5.start();
		Graphics2D g = (Graphics2D) screen.getGraphics();
		g.setColor(new Color(50, 50, 50));
		g.fillRect(0, 0, imgwidth, imgheight);
		double scalingconstant = 20.0*Math.pow(10.0, opts.gui_brightness.getValue()/10.0);
		int scalarview = opts.gui_view.getSelectedIndex();
		int vectorview = opts.gui_view_vec.getSelectedIndex();
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				int ex = (int) Math.round(this.ex[i][j]*scalingconstant);
				int ey = (int) Math.round(this.ey[i][j]*scalingconstant);
				int bz = (int) Math.round(this.bz[i][j]*c*scalingconstant);
				int rho = (int) Math.round(this.rho[i][j]*5*ds*scalingconstant);
				int jx = (int) Math.round(this.jx[i][j]*scalingconstant/Lambda);
				int jy = (int) Math.round(this.jy[i][j]*scalingconstant/Lambda);
				//int phi = (int) Math.round(this.poisson1[i][j]*scalingconstant);
				Color c = Color.WHITE;
				if (scalarview == 0) {
					//c = getcol(ex, 0, ey, 30);
					c = getcol(clamp(ex, 0, 999), 0, clamp(-ex, 0, 999), 30);
				}
				else if (scalarview == 1) {
					//c = getcol(0, bz, 0, 30);
					c = getcol(clamp(bz, 0, 999), 0, clamp(-bz, 0, 999), 30);
				}
				else if (scalarview == 2) {
					c = getcol(clamp(rho, 0, 999), 0, clamp(-rho, 0, 999), 30);
				}
				else if (scalarview == 3) {
					c = getcol(jx, 0, jy, 30);
				}
				else if (scalarview == 4) {
					c = getcol(clamp(rho, 0, 999), bz, clamp(-rho, 0, 999), 30);
				}
				Color d = new Color(180, 180, 180);
				double frac = 0.5*(1.0 - conductivityx[i][j]);
				g.setColor(interpcolor(c, d, frac));
				g.fillRect((int)(i*scalefactor), (int)(j*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
				
				//g.setColor();
				//g.fillRect((int)(i*scalefactor), (int)(j*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
				
				//g.drawLine(i*7, (int)(pos[i]*10.0+250), (i+1)*7, (int)(pos[i+1]*10.0+250));
				//g.setColor(Color.BLUE);
				//g.fillRect(i*7, (int)(pos[i]*10.0+250), 2, 2);
			}
		}
		
		double arrowlength = scalefactor;

		double vectorscalingconstant = 20.0*Math.pow(10.0, opts.gui_brightness_1.getValue()/10.0);
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
				double Estrength = vectorscalingconstant*Math.sqrt(arrow.dot(arrow));
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
				g.setColor(getcol((int)(Estrength*5), (int)(Estrength*5), (int)(Estrength*5), 30));
				g.drawLine((int)(body1.x*scalefactor), (int)(body1.y*scalefactor), (int)(body2.x*scalefactor), (int)(body2.y*scalefactor));
				g.drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip1.x*scalefactor), (int)(tip1.y*scalefactor));
				g.drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip2.x*scalefactor), (int)(tip2.y*scalefactor));
				//tip1.add(ctr);
			}
		}

		double brushsize = 1.5*Math.pow(10.0, opts.gui_brushsize.getValue()/10.0);
		int r = (int)(scalefactor*brushsize/ds);
		g.setColor(new Color(50, 50, 50));
		g.drawOval(mouseX - r+1, mouseY - r+1, 2*r-2, 2*r-2);
		g.setColor(new Color(200, 200, 200));
		g.drawOval(mouseX - r, mouseY - r, 2*r, 2*r);
		t5.stop();
		//g.setColor(.)
	}

	public Color getcol(int r, int g, int b, int base) {
		return new Color(clamp(Math.abs(r) + base, 0, 255), clamp(Math.abs(g) + base, 0, 255), clamp(Math.abs(b) + base, 0, 255));
	}
	
	public Color interpcolor(Color a, Color b, double frac) {
		return new Color((int)(b.getRed() * frac + a.getRed() * (1-frac)), (int)(b.getGreen() * frac + a.getGreen() * (1-frac)), (int)(b.getBlue() * frac + a.getBlue() * (1-frac)));
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
	
	int mouseX = 0;
	int mouseY = 250;
	int mouseButton = 0;
	boolean mouseIsPressed = false;
	PointerInfo p = MouseInfo.getPointerInfo();

	@Override
	public void mouseClicked(MouseEvent arg0) {}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent arg0) {
		mouseIsPressed = true;
		mouseButton = arg0.getButton();
	}
	@Override
	public void mouseReleased(MouseEvent arg0) {
		mouseIsPressed = false;}

	@Override
	public void mouseDragged(MouseEvent arg0) {
		mouseX = arg0.getX();
		mouseY = arg0.getY();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		mouseX = arg0.getX();
		mouseY = arg0.getY();
	}
	
	private Action key_p = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_paused.setSelected(!opts.gui_paused.isSelected());
        }
    };
    private Action key_f = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			advanceframe = true;
        }
    };
    private Action key_1 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(0);
        }
    };
    private Action key_2 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(1);
        }
    };
    private Action key_3 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(2);
        }
    };
    private Action key_4 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(3);
        }
    };
    private Action key_5 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_view.setSelectedIndex(4);
        }
    };
    private Action key_sh1 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(1);
        }
    };
    private Action key_sh2 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(2);
        }
    };
    private Action key_sh3 = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brush.setSelectedIndex(3);
        }
    };
    private Action key_c = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			clear = true;
        }
    };
    private Action key_r = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			reset = true;
        }
    };
    private Action key_m = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			measure = true;
        }
    };
    private Action key_eq = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brightness.setValue(opts.gui_brightness.getValue()+1);
        }
    };
    private Action key_minus = new AbstractAction(null) {
        @Override
        public void actionPerformed(ActionEvent e) {
			opts.gui_brightness.setValue(opts.gui_brightness.getValue()-1);
        }
    };
    
	public Action[] keys = {key_p, key_f, key_1, key_2, key_3, key_4, key_5, key_sh1, key_sh2, key_sh3, key_c, key_r, key_m, key_eq, key_minus};
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
			KeyStroke.getKeyStroke(KeyEvent.VK_MINUS, 0)};
    
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
		if (e.getSource() == opts.gui_resetall)
			reset = true;
	}

}

class RenderCanvas extends JPanel {

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
		tstart = System.nanoTime();
	}

	void disableOutput() {
		outputenabled = false;
	}
	void enableOutput() {
		outputenabled = true;
	}

	void stop() {
		long tend = System.nanoTime();
		long diff = tend - tstart;
		if (outputenabled)
			System.out.println(name + ", dt: " + String.format("%.2f", diff/1000000.0) + " ms");
	}

	void stop(String msg) {
		long tend = System.nanoTime();
		long diff = tend - tstart;
		if (outputenabled) {
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
