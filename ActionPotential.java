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

public class ActionPotential implements Runnable, MouseListener, MouseMotionListener, MouseWheelListener, ActionListener {
	
	double[][] ex; //x component of E field
	double[][] ey; //y component of E field
	double[][] bz; //z component of B field
	double[][] rhotot; //Net charge density
	double[][] jxtot; //Net current density, x component
	double[][] jytot; //Net current density, y component
	double[][][] rho; //Charge density of each ion type
	double[][][] jx; //Current density of each ion type, x component
	double[][][] jy; //Current density of each ion type, y component
	double[][] bzp; //Temporary storage for B field

	/*Used in poisson equation solver*/
	double[][] deltarho;
	double[][] deltarhoMG;
	double[][] poisson1;
	double[][] poisson2;
	
	double[][][] mobilityfactorx;
	double[][][] mobilityfactory;
	double[][][] emfx;
	double[][][] emfy;

	int chargetypes = 3;
	int c_sodium = 0;
	int c_potassium = 1;
	int c_chloride = 2;
	double[] diffusivities; //Diffusion constant for each ion
	double[] conductivities; //Molar conductivity of each ion
	double[] ecchargedensity; //Extracellular charge densities
	double[] ctchargedensity; //Cytoplasm charge densities
	
	/*Membrane model*/
	double[] NaActivation;
	double[] NaInactivation;
	double[] KActivation;
	
	/*Domain parameters*/
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
	int interiorpoints = (nx-2)*(ny-2);
	int membraneycoord_s = ny/2;
	int membraneycoord_e = ny/2+5;

	double uL = 5e-9;
	double uT = 1e-9;
	double uC = 1e-19;
	double faradayconst = 96485;
	double uV = uL*uL/(uC*uT*uT);
	
	double eps0 = 8.854e-12 * (uL*uL*uL)/(uT*uT*uC*uC);
	double c = 200;
	double mu0 = 1/(c*c*eps0);
	
	double csq = 1/(eps0*mu0);
	double q = 1.0;
	double sigma = c*(dt/ds);
	double sigmasq = sigma*sigma;
	double tauminus = 1-c*dt/ds;
	double tauplus = 1+c*dt/ds;
	double maxdt = 0.95*ds/(Math.sqrt(2)*c); //Set max dt to just below theoretical max, dictated by the CFL condition
	
	double time = 0.0;
	int lastsimspeed = 0;
	int framenumber = 0;
	
	//boolean paused = true;
	boolean advanceframe = false;
	boolean clear = false;
	boolean reset = false;
	boolean measure = false;
	boolean experimentended = false;
	double experimentduration = 99999.0;
	
	/*Graphics*/
	RenderCanvas r;
	MainWindow opts;
	Image screen;
	double scalefactor = ds*15.0;
	int imgwidth = 0;
	int imgheight = 0;
	
	/*Timing*/
	Timer t1 = new Timer("Poisson 1", false);
	Timer t2 = new Timer("Poisson 2", false);
	Timer t3 = new Timer("Poisson 3 intermediate", false);
	Timer t4 = new Timer("Poisson 3 total", true);
	Timer t5 = new Timer("Graphics", true);
	
	/*Mouse controls*/
	int brush_addcharge = 0;
	int brush_addinsulator = 1;
	int brush_addemf = 2;
	boolean drawingcharge = false;
	boolean linemode = false;
	int lineorientation = 0;
	int mouseX = 0;
	int mouseY = 250;
	int mouseXstart = 0;
	int mouseYstart = 250;
	double mxp = 0.0;
	double myp = 0.0;
	int mouseButton = 0;
	boolean mouseIsPressed = false;
	
	/*Debug*/
	double membranevoltage = 0;
	boolean vgpumpactive = false;

	public static void main(String[] args) {
		ActionPotential w = new ActionPotential();
		w.run();
	}
	
	public ActionPotential() {
		System.out.println("dtmax = " + maxdt);
		opts = new MainWindow();

		ex = new double[nx][ny];
		ey = new double[nx][ny];
		bz = new double[nx][ny];
		rhotot = new double[nx][ny];
		jxtot = new double[nx][ny];
		jytot = new double[nx][ny];
		rho = new double[chargetypes][nx][ny];
		jx = new double[chargetypes][nx][ny];
		jy = new double[chargetypes][nx][ny];

		bzp = new double[nx][ny];
		

		deltarho = new double[nx][ny];
		deltarhoMG = new double[nx][ny];
		poisson1 = new double[nx][ny];
		poisson2 = new double[nx][ny];
		
		mobilityfactorx = new double[chargetypes][nx][ny];
		mobilityfactory = new double[chargetypes][nx][ny];
		emfx = new double[chargetypes][nx][ny];
		emfy = new double[chargetypes][nx][ny];

		diffusivities = new double[chargetypes];
		conductivities = new double[chargetypes];
		ecchargedensity = new double[chargetypes];
		ctchargedensity = new double[chargetypes];

		NaActivation = new double[nx];
		NaInactivation = new double[nx];
		KActivation = new double[nx];

		setPhysicalConstants();
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
			int iterationmultiplier = opts.gui_simspeed_2.getValue();
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
	
	public void setPhysicalConstants() {
		
		diffusivities[c_sodium] = 1.96e-9;
		diffusivities[c_potassium] = 1.33e-9;
		diffusivities[c_chloride] = 2.03e-9;
		
		conductivities[c_sodium] = 73.6;
		conductivities[c_potassium] = 50.0;
		conductivities[c_chloride] = 76.2;
		
		ecchargedensity[c_sodium] = 440;
		ecchargedensity[c_potassium] = 20;
		ecchargedensity[c_chloride] = -460;
		
		ctchargedensity[c_sodium] = 50;
		ctchargedensity[c_potassium] = 400;
		ctchargedensity[c_chloride] = -450;
		
		for (int k = 0; k < chargetypes; k++) {
			diffusivities[k] *= uT/(uL*uL);
			conductivities[k] /= 10000*faradayconst*uT*uC;
			ecchargedensity[k] *= faradayconst*uL*uL*uL/uC;
			ctchargedensity[k] *= faradayconst*uL*uL*uL/uC;
		}
	}
	
	public void checkCFL() {
		if (c > ds/(Math.sqrt(2)*dt)) {
			System.out.println("Error: Wave speed too high. Ratio = " + (c/(ds/(Math.sqrt(2*dt)))));
			System.exit(0);
		}
		for (int k = 0; k < chargetypes; k++) {
			if (dt > ds*ds/(4*diffusivities[k])) {
				System.out.println("Error: Diffusion coefficient too high for k = " + k + ", Ratio = " + dt/(ds*ds/(4*diffusivities[k])));
				System.exit(0);
			}
		}
	}
	
	public void resetfields(boolean resetall) {
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				
				ex [i][j] = 0.0;
				ey [i][j] = 0.0;

				bz [i][j] = 0.0;
				bzp [i][j] = 0.0;


				jxtot [i][j] = 0.0;
				jytot [i][j] = 0.0;
				
				for (int k = 0; k < chargetypes; k++) {
					jx[k][i][j] = 0;
					jy[k][i][j] = 0;
				}


				poisson1 [i][j] = 0.0;
				poisson2 [i][j] = 0.0;
				
				if (resetall) {
					//double x = i*ds;
					//double y = j*ds;
					//double dist = Math.sqrt((x-25)*(x-25) + (y-25)*(y-25));
					//rho[i][j] = 100.0*Math.exp(-dist*dist);
					
					rhotot [i][j] = 0.0;

					for (int k = 0; k < chargetypes; k++) {
						rho[k][i][j] = 0;
						mobilityfactorx[k][i][j] = 1.0;
						mobilityfactory[k][i][j] = 1.0;
						emfx[k][i][j] = 0.0;
						emfy[k][i][j] = 0.0;
					}
				}
			}
		}


		if (resetall) {
			for (int j = 1; j < ny-1; j++) {
				for (int k = 0; k < chargetypes; k++) {
					mobilityfactorx[k][0][j] = 0.0;
					mobilityfactorx[k][nx-2][j] = 0.0;
				}
			}
			for (int i = 1; i < nx-1; i++) {
				for (int k = 0; k < chargetypes; k++) {
					mobilityfactory[k][i][0] = 0.0;
					mobilityfactory[k][i][ny-2] = 0.0;
				}
			}

			for (int i = 0; i < nx; i++) {
				NaActivation[i] = 0;
				NaInactivation[i] = 0;
				KActivation[i] = 0;
			}
			

			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					//if (i == 0 || i == nx - 1 || j == 0 || j == ny-1) {
						for (int k = 0; k < chargetypes; k++) {
							double f = (j - membraneycoord_s)/(double)(membraneycoord_e - membraneycoord_s);
							f = clamp(f, 0, 1);
							rho[k][i][j] = ecchargedensity[k]*f+ctchargedensity[k]*(1-f);
						}
					//}
				}
			}
		}

		/*for (int j = 1; j < ny-1; j++) {
			rho[0][0][j] = 1.0;
			rho[0][nx-1][j] = 1.0;
		}
		for (int i = 1; i < nx-1; i++) {
			rho[0][i][0] = 1.0;
			rho[0][i][ny-1] = 1.0;
		}*/
		

		correctEfield();
	}
	
	public void updateMembrane() {
		for (int i = 1; i < nx-1; i++) {
			// Na-K pump
			//emfy[c_sodium][i][membraneycoord] = 1.0*0;
			//emfy[c_potassium][i][membraneycoord] = -0.67*0;
			
			//double voltage = (ey[i][membraneycoord]*ds)*uV*1000;
			double voltage = 0;
			for (int j = membraneycoord_s; j < membraneycoord_e; j++) {
				voltage += (ey[i][j]*ds)*uV*1000;
			}
			//System.out.println(voltage);
			membranevoltage += voltage;
			double gNaLow = 0.001;
			double gNaHigh = 2.0;
			double gKLow = 0.5;
			double gKHigh = 2.0;
			double gCl = 0.00;
			
			double activationFactor = sigmoid((voltage+15.0)/20.0);
			double activationFactorK = sigmoid((voltage-15.0)/20.0);
			double inactivationFactor = sigmoid((voltage+50.0)/10.0);

			double factor = 5.0;
			double ra = factor*activationFactor/0.3;
			double rb = factor*(1-activationFactor)/1.2;
			double rc = factor*inactivationFactor/2.5;
			double rd = factor*(1-inactivationFactor)/0.6;

			double re = factor*activationFactorK/0.5;
			double rf = factor*(1-activationFactorK)/1.5;
			
			double atmp = NaActivation[i];
			double iatmp = NaInactivation[i];
			double katmp = KActivation[i];
			
			NaActivation[i] += dt*(ra*(1-atmp-iatmp) - rb*atmp + rd*iatmp - rc*atmp);
			NaInactivation[i] += dt*(rc*atmp-rd*iatmp);
			KActivation[i] += dt*(re*(1-katmp) - rf*katmp);
			
			double NaFactor = clamp(4*Math.pow(NaActivation[i], 3), 0, 100);
			double KFactor = clamp(Math.pow(KActivation[i], 3), 0, 100);
			if (!vgpumpactive) {
				NaFactor = 0;
				KFactor = 0;
			}
			for (int j = membraneycoord_s; j < membraneycoord_e; j++) {
			mobilityfactory[c_sodium][i][j] = gNaLow + (gNaHigh-gNaLow)*NaFactor;
			mobilityfactory[c_potassium][i][j] = gKLow + (gKHigh-gKLow)*KFactor;
			mobilityfactory[c_chloride][i][j] = gCl;
			}
		}
		membranevoltage /= (nx-2);
	}
	
	public double sigmoid(double x) {
		return 1.0/(1.0+Math.exp(-4.0*x));
	}
	
	public void handleMouseInput() {
		int brush = opts.gui_brush.getSelectedIndex();
		int brushshape = opts.gui_brush_1.getSelectedIndex();
		int k = opts.gui_chargetype.getSelectedIndex();
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
			double intensity = 1.0*Math.pow(10.0, opts.gui_brushintensity2.getValue()/10.0);
			
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
					double r = 0;
					if (brushshape == 0)
						r = Math.sqrt(p.dot(p));
					else
						r = Math.max(Math.abs(p.x), Math.abs(p.y));
					if (r <= brushsize) {
						if (brush == brush_addinsulator) {
							double newconductivity = 0.0;
							if (mouseButton == MouseEvent.BUTTON3)
								newconductivity = 1;
							for (int k2 = 0; k2 < chargetypes; k2++) {
								mobilityfactorx[k2][i][j] = newconductivity;
								mobilityfactory[k2][i][j] = newconductivity;
								mobilityfactorx[k2][i-1][j] = newconductivity;
								mobilityfactory[k2][i][j-1] = newconductivity;
							}
						} else if (brush == brush_addcharge) {
							double drho = 0.0;
							if (mouseButton == MouseEvent.BUTTON1)
								drho = intensity*ecchargedensity[k]*0.01;
							else if (mouseButton == MouseEvent.BUTTON3)
								drho = -intensity*ecchargedensity[k]*0.01;
							rho[k][i][j] += drho*Math.exp(-3.0*r/brushsize);
							rhotot[i][j] += drho*Math.exp(-3.0*r/brushsize);
						} else if (brush == brush_addemf) {
							if (mouseButton == MouseEvent.BUTTON1) {
								emfx[k][i][j] = 10.0*intensity*Math.cos(angle);
								emfx[k][i][j+1] = 10.0*intensity*Math.cos(angle);
								emfy[k][i][j] = 10.0*intensity*Math.sin(-angle);
								emfy[k][i+1][j] = 10.0*intensity*Math.sin(-angle);
							} else if (mouseButton == MouseEvent.BUTTON3) {
								emfx[k][i][j] = 0;
								emfx[k][i][j+1] = 0;
								emfy[k][i][j] = 0;
								emfy[k][i+1][j] = 0;
							}
						}
					}
				}
			}

		}
		mxp = mx;
		myp = my;
	}
	
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
			checkCFL();
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
						rhotot[i][j] = 0;
						for (int k = 0; k < chargetypes; k++) {
							rho[k][i][j] = rho[k][i][j] - (jx[k][i][j]-jx[k][i-1][j] + jy[k][i][j]-jy[k][i][j-1])*dt/ds;
							rhotot[i][j] += rho[k][i][j];
						}
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
				
				updateMembrane();

				double[] tmpconductivity = new double[chargetypes];
				
				for (int i = 0; i < nx-1; i++)
				{
					for (int j = 1; j < ny-1; j++)
					{
						jxtot[i][j] = 0;
						double totalconductivity = 0;
						
						for (int k = 0; k < chargetypes; k++)
						{
							tmpconductivity[k] = 0.5*(rho[k][i+1][j] + rho[k][i][j])*conductivities[k]*mobilityfactorx[k][i][j];
							if (tmpconductivity[k] < 0)
								tmpconductivity[k] = 0;
							jx[k][i][j] = (-diffusivities[k]*(rho[k][i+1][j] - rho[k][i][j])/ds)*mobilityfactorx[k][i][j] + emfx[k][i][j] + 0.5*tmpconductivity[k]*ex[i][j];
							totalconductivity += tmpconductivity[k];
							jxtot[i][j] += jx[k][i][j];
						}
						
						ex[i][j] = (ex[i][j] + csq*(bz[i][j]-bz[i][j-1])*dt/ds - jxtot[i][j]/eps0*dt)/(1+0.5*dt*totalconductivity/eps0);

						for (int k = 0; k < chargetypes; k++)
						{
							jx[k][i][j] += 0.5*ex[i][j]*tmpconductivity[k];
						}
						
						jxtot[i][j] += 0.5*ex[i][j]*totalconductivity;
					}
				}

				for (int i = 1; i < nx-1; i++)
				{
					for (int j = 0; j < ny-1; j++)
					{
						jytot[i][j] = 0;
						double totalconductivity = 0;
						
						for (int k = 0; k < chargetypes; k++)
						{
							tmpconductivity[k] = 0.5*(rho[k][i][j+1] + rho[k][i][j])*conductivities[k]*mobilityfactory[k][i][j];
							if (tmpconductivity[k] < 0)
								tmpconductivity[k] = 0;
							jy[k][i][j] = (-diffusivities[k]*(rho[k][i][j+1] - rho[k][i][j])/ds)*mobilityfactory[k][i][j] + emfy[k][i][j] + 0.5*tmpconductivity[k]*ey[i][j];
							totalconductivity += tmpconductivity[k];
							jytot[i][j] += jy[k][i][j];
						}
						
						ey[i][j] = (ey[i][j] + -csq*(bz[i][j]-bz[i-1][j])*dt/ds - jytot[i][j]/eps0*dt)/(1+0.5*dt*totalconductivity/eps0);
						
						for (int k = 0; k < chargetypes; k++)
						{
							jy[k][i][j] += 0.5*ey[i][j]*tmpconductivity[k];
						}
						
						jytot[i][j] += 0.5*ey[i][j]*totalconductivity;
					}
				}
			}
			
			if (framenumber%40 == 0) {
				correctEfield();
			}
			
			// bz[100][100] = 10.0*Math.sin(time*100.0);
			
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
				deltarho[i][j] = eps0*((ex[i][j]-ex[i-1][j] + ey[i][j]-ey[i][j-1])/ds) - rhotot[i][j];
			}
		}
		
		double err = 0;
		for (int i = 1; i < nx-1; i++) {
			for (int j = 1; j < ny-1; j++) {
				err += Math.abs((4*poisson1[i][j] - (poisson1[i-1][j] + poisson1[i+1][j] + poisson1[i][j-1] + poisson1[i][j+1]))/(ds*ds) - deltarho[i][j]/eps0);
			}
		}
		System.out.println("Poisson init error: " + err);
		
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
		System.out.println("Poisson final error: " + err);
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
	
	public double getrhototal(double x, double y) {
		return bilinearinterp(rhotot, x, y);
	}

	public double getjx(int k, double x, double y) {
		return bilinearinterp(jx[k], x - 0.5, y);
	}
	
	public double getjy(int k, double x, double y) {
		return bilinearinterp(jy[k], x, y - 0.5);
	}
	
	public double getjxtotal(double x, double y) {
		return bilinearinterp(jxtot, x - 0.5, y);
	}
	
	public double getjytotal(double x, double y) {
		return bilinearinterp(jytot, x, y - 0.5);
	}
	
	public void render() {
		t5.start();
		Graphics2D g = (Graphics2D) screen.getGraphics();

		((Graphics2D)g).setRenderingHint(
		        RenderingHints.KEY_ANTIALIASING,
		        RenderingHints.VALUE_ANTIALIAS_ON);
		
		g.setColor(new Color(50, 50, 50));
		g.fillRect(0, 0, imgwidth, imgheight);
		double scalingconstant = 20.0*Math.pow(10.0, opts.gui_brightness.getValue()/10.0);
		int scalarview = opts.gui_view.getSelectedIndex();
		int vectorview = opts.gui_view_vec.getSelectedIndex();
		int k = opts.gui_chargetype.getSelectedIndex();
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				double x = i;
				double y = j;
				
				int ex = (int) Math.round(this.ex[i][j]*scalingconstant);
				int ey = (int) Math.round(this.ey[i][j]*scalingconstant);
				int bz = (int) Math.round(this.bz[i][j]*c*scalingconstant);
				int rho = (int) Math.round(this.rho[k][i][j]*5*ds*scalingconstant);
				int jx = (int) Math.round(this.jx[k][i][j]*scalingconstant);
				int jy = (int) Math.round(this.jy[k][i][j]*scalingconstant);
				//int phi = (int) Math.round(this.poisson1[i][j]*scalingconstant);
				Color c = Color.WHITE;
				if (scalarview == 0) {
					//c = getcol(ex, 0, ey, 30);
					c = getcol(ex, 0, ey, 30);
				}
				else if (scalarview == 1) {
					//c = getcol(0, bz, 0, 30);
					c = getcol(clamp(bz, 0, 999), 0, clamp(-bz, 0, 999), 30);
					x += 0.5;
					y += 0.5;
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
				double frac = 0.5*(1.0 - mobilityfactory[0][i][j]);
				//g.setColor(interpcolor(c, d, frac));
				g.setColor(c);
				g.fillRect((int)((x-0.5)*scalefactor), (int)((y-0.5)*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
				
				//g.setColor();
				//g.fillRect((int)(i*scalefactor), (int)(j*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
				
				//g.drawLine(i*7, (int)(pos[i]*10.0+250), (i+1)*7, (int)(pos[i+1]*10.0+250));
				//g.setColor(Color.BLUE);
				//g.fillRect(i*7, (int)(pos[i]*10.0+250), 2, 2);
			}
		}
		

		for (int i = 0; i < nx; i++) {
				double x = i;
				double y = membraneycoord_s;
				int c1 = (int) Math.round(this.NaActivation[i]*256);
				int c2 = (int) Math.round(this.NaInactivation[i]*256);
				int c3 = (int) Math.round(this.KActivation[i]*256);
				Color d = getcol(c1, c2, c3, 50);
				g.setColor(d);
				g.fillRect((int)((x-0.5)*scalefactor), (int)((y-0.5)*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
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
					arrow = new Vector(getex(x, y)*ds*uV*1000, getey(x, y)*ds*uV*1000);
				else if (vectorview == 1)
					arrow = new Vector(getjx(k, x, y), getjy(k, x, y));
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

		
		g.drawString("Voltage: " + membranevoltage, 20, 20);
		
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
			vgpumpactive = !vgpumpactive;
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
		if (e.getSource() == opts.gui_resetall)
			reset = true;
	}

}

class RenderCanvas extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7369516276529576171L;
	ActionPotential parent;
	@Override
	public void paintComponent(Graphics real) {
		parent.render();
		real.drawImage(parent.screen, 0, 0, parent.opts);
	}
	
	public RenderCanvas(ActionPotential w) {
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
