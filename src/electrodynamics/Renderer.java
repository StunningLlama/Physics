// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferInt;
import java.util.ArrayList;
import java.util.Random;
import java.util.TimerTask;
import java.util.concurrent.BrokenBarrierException;

import javax.swing.JPanel;

import electrodynamics.Controls.Brush;
import electrodynamics.plot.Plot;
import electrodynamics.util.DistributionSampler;
import electrodynamics.util.FastRandom;
import electrodynamics.util.Timer;
import electrodynamics.util.Utils;
import electrodynamics.util.Vector;

public class Renderer extends TimerTask {
	Simulation e;
	
	/* Graphics */
	
	BufferedImage img_back;
	BufferedImage img_front;
	public int[] imgData;
	public double[][] scalarfield;
	public double[][] gradscalarfield;
	public float[][] image_r;
	public float[][] image_g;
	public float[][] image_b;
	public float col_r = 0;
	public float col_g = 0;
	public float col_b = 0;
	public float alphaBG = 0;
	public float alphaFG = 0;

	public ArrayList<Text> texts = new ArrayList<>();
	public Font bigfont = new Font(Font.SANS_SERIF, Font.PLAIN, 15);
	public Font regularfont = new Font(Font.SANS_SERIF, Font.PLAIN, 12);
	public int scalefactor;
	public int imgwidth = 0;
	public int imgheight = 0;
	public int targetframerate = 60;
	public int frameduration = 1000/targetframerate;
	
	double t_prev = 0;
	double delta_t = 0;

	ArrayList<Dot> dots = new ArrayList<>();
	ArrayList<ChargeCarrierDot> ccdots = new ArrayList<>();
	int numdots = 2000;
	double rho_n_max = 0;
	double rho_p_max = 0;
	double C_prev = 0;

	int probetexttimer = 0;


	// Electron and hole dots
	
	DistributionSampler rho_p_dist = new DistributionSampler();
	DistributionSampler rho_n_dist = new DistributionSampler();
	DistributionSampler rho_G_dist = new DistributionSampler();
	FastRandom frand = new FastRandom();
	double cc_default_dot_density = 5e11;		// How many electron/hole dots to draw
	double tau = 1e-12;		// How long a dot stays on the screen
	
	// Performance profiling
	
	Timer FPStimer = new Timer("Graphics FPS", 10, true);
	Timer t5 = new Timer("Graphics", 20, true);
	
	public Renderer(Simulation e) {
		this.e = e;
	}
	
	public void setResolution(int resolution) {
		scalarfield = new double[e.nx][e.ny];
		gradscalarfield = new double[e.nx][e.ny];
		image_r = new float[e.nx][e.ny];
		image_g = new float[e.nx][e.ny];
		image_b = new float[e.nx][e.ny];
		scalefactor = (int)(768.0/e.ny);
		if (scalefactor < 1) scalefactor = 1;

		rho_p_dist.init(e.nx, e.ny);
		rho_n_dist.init(e.nx, e.ny);
		rho_G_dist.init(e.nx, e.ny);

		imgwidth = (int)Math.ceil(scalefactor*e.nx);
		imgheight = (int)Math.ceil(scalefactor*e.ny);

		e.canvas.setPreferredSize(new Dimension(imgwidth, imgheight));
		img_back = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
		img_front = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
		imgData = ((DataBufferInt)img_back.getRaster().getDataBuffer()).getData();
		
		reset();
	}
	
	public void reset() {
		t_prev = e.time;
		
		dots.clear();
		for (int i = 0; i < numdots; i++) {
			dots.add(new Dot(0, 0, 0, 0));
		}
		
		ccdots.clear();
		C_prev = 0;
	}
	
	
	public void setalphaBG(double alpha) {
		alphaBG = (float)alpha;
	}

	public void setalphaFG(double alpha) {
		alphaFG = (float)alpha;
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

	public void setPixel(int i, int j) {
		if (i < 0 || j < 0 || i >= e.nx || j >= e.ny)
			return;

		image_r[i][j] = image_r[i][j]*alphaBG + col_r*alphaFG;
		image_g[i][j] = image_g[i][j]*alphaBG + col_g*alphaFG;
		image_b[i][j] = image_b[i][j]*alphaBG + col_b*alphaFG;
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

	public void drawPixelRectangle(int x, int y, int w, int h) {
		for (int i = x; i < x+w; i++) {
			for (int j = y; j < y+h; j++) {
				setPixel(i, j);
			}
		}
	}

	public void stampPixelData() {
		e.t9.start();
		int scansize = e.nx*scalefactor;
		for (int x = 0; x < e.nx*scalefactor; x++) {
			for (int y = 0; y < e.ny*scalefactor; y++) {
				int i = x/scalefactor;
				int j = y/scalefactor;
				double scale = 1f/max(image_r[i][j], image_g[i][j], image_b[i][j], 1f);
				int rgb = clamp((int)(256*image_r[i][j]*scale), 0, 255) << 16
						| clamp((int)(256*image_g[i][j]*scale), 0, 255) << 8
						| clamp((int)(256*image_b[i][j]*scale), 0, 255);
				imgData[x + y*scansize] = rgb;
			}
		}
		e.t9.stop();
	}

	public void drawPixel(int x, int y) {
		drawPixel(x, y, col_r, col_g, col_b, alphaFG, alphaBG);
	}

	public void drawPixel(int x, int y, float col_r, float col_g, float col_b, float alphaFG, float alphaBG) {
		int scansize = scalefactor*e.nx;
		int rgb =  imgData[x+y*scansize];
		int r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
		int g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
		int b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
		if (r > 255)
			r = 255;
		if (g > 255)
			g = 255;
		if (b > 255)
			b = 255;
		imgData[x+y*scansize] = 255<<24 | r << 16 | g << 8 | b;
	}

	public void drawPixelWithContrast(int x, int y, float col_r, float col_g, float col_b, float alphaFG, float alphaBG) {
		int scansize = scalefactor*e.nx;
		int rgb =  imgData[x+y*scansize];
		if (((rgb>>16)&255) + ((rgb>>8)&255) + ((rgb)&255) > 512) {
			alphaBG = 1-alphaFG;
			alphaFG = 0;
		}

		int r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
		int g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
		int b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
		if (r > 255)
			r = 255;
		if (g > 255)
			g = 255;
		if (b > 255)
			b = 255;
		imgData[x+y*scansize] = 255<<24 | r << 16 | g << 8 | b;
	}

	public void drawRectangle(int x, int y, int w, int h) {
		drawRectangle(x, y, w, h, col_r, col_g, col_b, alphaFG, alphaBG);
	}

	public void drawRectangle(int x, int y, int w, int h, float col_r, float col_g, float col_b, float alphaFG, float alphaBG) {
		try {
			int scansize = scalefactor*e.nx;
			for (int i = x; i < x+w; i++) {
				for (int j = y; j < y+h; j++) {
					int rgb =  imgData[i+j*scansize];
					int r = ((int)(((rgb>>16)&255)*alphaBG + 255*col_r*alphaFG));
					int g = ((int)(((rgb>>8)&255)*alphaBG + 255*col_g*alphaFG));
					int b = ((int)(((rgb)&255)*alphaBG + 255*col_b*alphaFG));
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					imgData[i+j*scansize] = 255<<24 | r << 16 | g << 8 | b;
				}
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			return;
		}
	}

	public void drawLine(int x0, int y0, int x1, int y1, boolean draw_starting_point) {
		drawLine(x0, y0, x1, y1, draw_starting_point, col_r, col_g, col_b, alphaFG, alphaBG);
	}

	public void drawLine(int x0, int y0, int x1, int y1, boolean draw_starting_point, float col_r, float col_g, float col_b, float alphaFG, float alphaBG) {
		try {
			int scansize = scalefactor*e.nx;

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
		if ((val != val) || (val < min)) return min;
		if (val > max) return max;
		return val;
	}

	public double bump(double x, double w) {
		if (x <= 0 || x >= 1) return 0;
		if (x <= (1-w) && x >= w) return 1;

		double y = x/w;
		if (x > (1-w))
			y = (1.0-x)/w;

		return 3*y*y - 2*y*y*y;
	}

	@Override
	public void run() {
		e.rwLock.readLock().lock();
		try {
			t5.start();
			drawPixels();
			drawVectors();
			drawText();
			copyImage(img_back, img_front);
			e.canvas.repaint();
			t5.stop();

			FPStimer.stop();
			FPStimer.start();
		} catch (Exception e1) {
			e.displayErrorMessage(e1);
		}
		finally {
			e.rwLock.readLock().unlock();
		}
	}
	
	BufferedImage copyImage(BufferedImage src, BufferedImage dst) {
	    if (src.getType() != dst.getType() || 
	        src.getWidth() != dst.getWidth() || 
	        src.getHeight() != dst.getHeight()) {
	        throw new IllegalArgumentException("Images must be same size and type");
	    }

	    DataBuffer srcBuffer = src.getRaster().getDataBuffer();
	    DataBuffer dstBuffer = dst.getRaster().getDataBuffer();

	    if (srcBuffer instanceof DataBufferInt && dstBuffer instanceof DataBufferInt) {
	        int[] srcData = ((DataBufferInt) srcBuffer).getData();
	        int[] dstData = ((DataBufferInt) dstBuffer).getData();
	        System.arraycopy(srcData, 0, dstData, 0, srcData.length);
	    } else {
	        // Fallback if not int-packed
	        int[] temp = src.getRaster().getPixels(0, 0, src.getWidth(), src.getHeight(), (int[]) null);
	        dst.getRaster().setPixels(0, 0, dst.getWidth(), dst.getHeight(), temp);
	    }
	    return dst;
	}
	
	void drawPixels() {

		for (int i = 0; i < e.nx; i++) {
			for (int j = 0; j < e.ny; j++) {
				image_r[i][j] = 0;
				image_g[i][j] = 0;
				image_b[i][j] = 0;
			}
		}

		setalphaBG(0);
		setalphaFG(1);

		Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();

		if (e.opts.gui_elem_colors.isSelected()) {
			for (int i = 0; i < e.nx; i++) {
				for (int j = 0; j < e.ny; j++) {
					setColor(e.materials[i][j].type.color_r, e.materials[i][j].type.color_g, e.materials[i][j].type.color_b);

					setPixel(i, j);
				}
			}
		} else {
			for (int i = 0; i < e.nx; i++) {
				for (int j = 0; j < e.ny; j++) {
					setColor(e.materials[i][j].type.color_grayscale, e.materials[i][j].type.color_grayscale, e.materials[i][j].type.color_grayscale);
					setPixel(i, j);
				}
			}
		}

		if (e.controls.moving_selection) {
			if (e.opts.gui_elem_colors.isSelected()) {
				for (int i = 0; i < e.nx; i++) {
					for (int j = 0; j < e.ny; j++) {
						int si = i-e.controls.delta_mx_index;
						int sj = j-e.controls.delta_my_index;
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection[si][sj].type != MaterialType.VACUUM) {
							setColor(e.controls.selection[si][sj].type.color_r, e.controls.selection[si][sj].type.color_g, e.controls.selection[si][sj].type.color_b);
							setPixel(i, j);
						}
					}
				}
			} else {
				for (int i = 0; i < e.nx; i++) {
					for (int j = 0; j < e.ny; j++) {
						int si = i-e.controls.delta_mx_index;
						int sj = j-e.controls.delta_my_index;
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection[si][sj].type != MaterialType.VACUUM) {
							setColor(e.controls.selection[si][sj].type.color_grayscale, e.controls.selection[si][sj].type.color_grayscale, e.controls.selection[si][sj].type.color_grayscale);
							setPixel(i, j);
						}
					}
				}
			}
		}

		ScalarView scalarview = (ScalarView) e.opts.gui_view.getSelectedItem();

		if (scalarview != ScalarView.NONE) {
			setalphaBG(1.0);
			setalphaFG(1.0);
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					switch (scalarview) {
					case NONE:
						break;
					case B_FIELD:
						scalarfield[i][j] = e.parity*0.25*(e.Bz[i][j]+e.Bz[i-1][j]+e.Bz[i][j-1]+e.Bz[i-1][j-1]);
						break;
					case E_FIELD:
						scalarfield[i][j] = Utils.length(0.5*(e.Ex[i][j]+e.Ex[i][j+1]), 0.5*(e.Ey[i][j]+e.Ey[i+1][j]));
						break;
					case H_FIELD:
						scalarfield[i][j] = e.parity*0.25*(e.Hz[i][j]+e.Hz[i-1][j]+e.Hz[i][j-1]+e.Hz[i-1][j-1]);
						break;
					case CURRENT:
						scalarfield[i][j] = Utils.length(0.5*(e.Jx_free[i][j]+e.Jx_free[i-1][j]), 0.5*(e.Jy_free[i][j]+e.Jy_free[i][j-1]));
						break;
					case POTENTIAL:
						scalarfield[i][j] = e.phi[i][j];
						break;
					case CHARGE:
						scalarfield[i][j] = e.rho_free[i][j];
						break;
					case BACKGROUND_CHARGE:
						scalarfield[i][j] = e.rho_back[i][j];
						break;
					case ELECTRON_CHARGE:
						scalarfield[i][j] = e.rho_n[i][j];
						break;
					case HOLE_CHARGE:
						scalarfield[i][j] = e.rho_p[i][j];
						break;
					case COMBINED_CHARGE:
						scalarfield[i][j] = Double.NaN;
						break;
					case ENERGY:
						scalarfield[i][j] = e.u[i][j];
						break;
					case ENTROPY:
						scalarfield[i][j] = e.S[i][j];
						break;
					case HEAT:
						scalarfield[i][j] = e.Q[i][j];
						break;
					case ELECTRON_POTENTIAL:
						scalarfield[i][j] = e.conducting[i][j]*(e.F_n[i][j]/e.q_n+e.phi[i][j]-e.W_semi/e.eVtoJ);
						break;
					case HOLE_POTENTIAL:
						scalarfield[i][j] = e.conducting[i][j]*(e.F_p[i][j]/e.q_p+e.phi[i][j]-e.W_semi/e.eVtoJ);
						break;
					case DEBUG:
						scalarfield[i][j] = (e.visited[i][j]? 1:0);
						break;
					case GENERATION:
						scalarfield[i][j] = e.G[i][j];
						break;
					case RECOMBINATION:
						scalarfield[i][j] = e.R[i][j];
						break;
					case AVERAGE_POTENTIAL:
						scalarfield[i][j] = e.F[i][j];
						break;
					case LIGHT:
						scalarfield[i][j] = -e.materials[i][j].semiconducting*(e.G[i][j]-e.R[i][j]);
						break;
					}
				}
			}
		}

		for (int i = 1; i < e.nx-1; i++) {
			for (int j = 1; j < e.ny-1; j++) {
				gradscalarfield[i][j] = Utils.length(scalarfield[i+1][j]-scalarfield[i-1][j], scalarfield[i][j+1]-scalarfield[i][j-1])/(2*e.ds);
			}
		}
		ScalarView.ColorScheme colorscheme = ((ScalarView) e.opts.gui_view.getSelectedItem()).colorscheme;
		float scalingconstant = (float) (10.0*Math.pow(10.0, e.opts.gui_brightness.getValue()/10.0)/scalarview.scale);


		/*double scale_min = -1/scalingconstant;
	double scale_max = 1/scalingconstant;
	for (int i = e.nx-20; i < e.nx; i++) {
		double t = (i-(e.nx - 20))/(double)(e.nx - (e.nx-20));
		scalarfield[i][2] = scale_min + t*(scale_max - scale_min);
	}*/


		if (colorscheme == ScalarView.ColorScheme.RED_BLUE) {
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					setColorFloat((float) scalarfield[i][j]*scalingconstant, 0, -(float) scalarfield[i][j]*scalingconstant);
					setPixel(i, j);
				}
			}
		} else if (colorscheme == ScalarView.ColorScheme.CYAN_YELLOW) {
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					setColorFloat((float) scalarfield[i][j]*scalingconstant, Math.abs((float) scalarfield[i][j]*scalingconstant), -(float) scalarfield[i][j]*scalingconstant);
					setPixel(i, j);
				}
			}
		} else if (colorscheme == ScalarView.ColorScheme.GREEN) {
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					setColorFloat(0, Math.abs((float) scalarfield[i][j]*scalingconstant), 0);
					setPixel(i, j);
				}
			}
		} else if (colorscheme == ScalarView.ColorScheme.WHITE) {
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					setColorFloat((float) scalarfield[i][j]*scalingconstant, (float) scalarfield[i][j]*scalingconstant, (float) scalarfield[i][j]*scalingconstant);
					setPixel(i, j);
				}
			}
		} else if (colorscheme == ScalarView.ColorScheme.OTHER) {
			if (scalarview == ScalarView.COMBINED_CHARGE) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float rc = (float)(clamp(0.2*Math.log(e.rho_p[i][j]*scalingconstant), 0, 1));
						float bc = (float)(clamp(0.2*Math.log(-e.rho_n[i][j]*scalingconstant), 0, 1));
						float gc = Math.min(rc, bc);
						double it = Math.max(rc, bc);
						setalphaBG(1-0.5*gc);
						setalphaFG(0.6*it);
						setColorFloat(rc, gc, bc);
						setPixel(i, j);
					}
				}
			}  else if (scalarview == ScalarView.CHARGE) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float rc = Math.min(Math.max((float) scalarfield[i][j]*scalingconstant, 0), 1);
						float bc = Math.min(Math.max(-(float) scalarfield[i][j]*scalingconstant, 0), 1);
						float gc = Math.min(rc, bc);

						setalphaBG(1-0.5*gc);
						setColorFloat(rc, gc, bc);
						setPixel(i, j);
					}
				}
			}
		}

		boolean highlight = (e.opts.gui_brush_highlight.isSelected() && Brush.isBrushShapeImportant(brush));
		for (int i = 0; i < e.nx; i++) {
			for (int j = 0; j < e.ny; j++) {
				if (e.materials[i][j].type == MaterialType.EMF && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
					setalphaBG(0.25);
					setalphaFG(0.75);

					int offset = 10*(2*((i+j)%2)-1);
					if (e.controls.selected_EMF[i][j])
						offset = 60*(2*((i+j)%2)-1)-50;

					if ((e.materials[i+1][j].type != MaterialType.EMF
					|| e.materials[i-1][j].type != MaterialType.EMF
					|| e.materials[i][j+1].type != MaterialType.EMF
					| e.materials[i][j-1].type != MaterialType.EMF))
					{
						offset = -30;
					}

					int delta_r = MaterialType.EMF.color_r+offset;
					int delta_g = MaterialType.EMF.color_g+offset;
					int delta_b = MaterialType.EMF.color_b+offset;
					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				} else if (e.materials[i][j].type == MaterialType.AC_EMF && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
					setalphaBG(0.25);
					setalphaFG(0.75);

					int offset = 10*(2*((i+j)%2)-1);
					if (e.controls.selected_EMF[i][j])
						offset = 60*(2*((i+j)%2)-1)-50;

					if ((e.materials[i+1][j].type != MaterialType.AC_EMF
					|| e.materials[i-1][j].type != MaterialType.AC_EMF
					|| e.materials[i][j+1].type != MaterialType.AC_EMF
					| e.materials[i][j-1].type != MaterialType.AC_EMF))
					{
						offset = -30;
					}

					int delta_r = MaterialType.AC_EMF.color_r+offset;
					int delta_g = MaterialType.AC_EMF.color_g+offset;
					int delta_b = MaterialType.AC_EMF.color_b+offset;
					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				} else if (e.materials[i][j].type == MaterialType.SWITCH && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
					setalphaBG(0.25);
					setalphaFG(0.75);

					int offset = 10*(2*((i+j)%2)-1);

					if ((e.materials[i+1][j].type != MaterialType.SWITCH
					|| e.materials[i-1][j].type != MaterialType.SWITCH
					|| e.materials[i][j+1].type != MaterialType.SWITCH
					| e.materials[i][j-1].type != MaterialType.SWITCH))
					{
						offset = -30;
					}

					if (e.materials[i][j].activated == 0)
						offset -= 200;

					int delta_r = MaterialType.SWITCH.color_r+offset;
					int delta_g = MaterialType.SWITCH.color_g+offset;
					int delta_b = MaterialType.SWITCH.color_b+offset;
					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				} else if (!(e.materials[i][j].activated == 1)) {
					setalphaBG(0.25);
					setalphaFG(0.75);
					setColor(0, 0, 0);
					setPixel(i, j);
				}

				if (e.controls.selected[i][j] || (highlight && e.controls.under_brush[i][j]))
				{
					setalphaBG(0.75);
					setalphaFG(0.25);

					int s = e.controls.selected[i][j]? 1:0;
					int h = (highlight && e.controls.under_brush[i][j])? 1:0;
					int delta_r = 256*s + 256*h;
					int delta_g = 100*s + 256*h;
					int delta_b = 256*s + 256*h;

					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				}

				if (e.L[i][j] > 0)
				{
					setalphaBG(0.5);
					setalphaFG(0.75);
					setColorFloat(1, 1, 1);
					setPixel(i, j);
				}
			}
		}


		if (e.opts.gui_interface.isSelected())
		{
			setalphaBG(1);
			setalphaFG(1);
			setColorFloat(0.7f, 0.7f, 0.7f);

			if (Brush.drawLine(brush) && e.controls.mouse_pressed) {
				drawPixelLine(e.controls.mx_start_index, e.controls.my_start_index, e.controls.mx_index, e.controls.my_index);
			}


			setalphaFG(1.0);
			setColorFloat(0.5f, 1.0f, 1.0f);

			for (CurrentProbe p: e.currentprobes) {
				drawPixelRectangle(p.x1-1, p.y1-1, 3, 3);
				drawPixelRectangle(p.x2-1, p.y2-1, 3, 3);
			}

			for (VoltageProbe p: e.voltageprobes) {
				drawPixelRectangle(p.x-1, p.y-1, 3, 3);
			}

			if (e.ground != null) {
				drawPixelRectangle(e.ground.x-1, e.ground.y-1, 3, 3);
			}

			setalphaFG(0.3);
			setColorFloat(0.5f, 1.0f, 1.0f);

			for (CurrentProbe p: e.currentprobes) {
				drawPixelLine(p.x1, p.y1, p.x2, p.y2);
			}

			setalphaFG(1.0);
			setColorFloat(1.0f, 1.0f, 1.0f);

			for (Plot p: e.plots) {
				if (p.frame.isVisible()) {
					drawPixelRectangle((int)p.x1-1, (int)p.y1-1, 3, 3);
					drawPixelRectangle((int)p.x2-1, (int)p.y2-1, 3, 3);
				}
			}

			setalphaFG(1.0);
			setColorFloat(1.0f, 1.0f, 1.0f);

			for (Plot p: e.plots) {
				if (p != null && p.frame.isVisible()) {
					drawPixelLine((int)p.x1, (int)p.y1, (int)p.x2, (int)p.y2);
				}
			}

			setalphaFG(0.8);
			setColorFloat(1.0f, 1.0f, 1.0f);

			if (e.controls.texting) {
				drawPixelLine(e.controls.text_x, e.controls.text_y, e.controls.text_x, e.controls.text_y+7);
			}
		}

		stampPixelData();

	}
	
	void drawVectors() {
		delta_t = e.time - t_prev;
		t_prev = e.time;

		if ((VectorMode) e.opts.gui_view_vec_mode.getSelectedItem() == VectorMode.SPECIES && !e.opts.gui_paused.isSelected()) {

			double C = cc_default_dot_density*Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/20.0);

			rho_n_dist.prepare(e.rho_n);
			rho_p_dist.prepare(e.rho_p);
			rho_G_dist.prepare(e.G);

			double A = e.ds*e.ds;
			double N_G = rho_G_dist.getTotalAmount()*A*e.e_charge*C*delta_t;
			double N_n = rho_n_dist.getTotalAmount()*A*C*delta_t/tau;
			double N_p = rho_p_dist.getTotalAmount()*A*C*delta_t/tau;

			double N_n_excess = Math.max(C-C_prev, 0)*rho_n_dist.getTotalAmount()*A;
			double N_p_excess = Math.max(C-C_prev, 0)*rho_p_dist.getTotalAmount()*A;

			double P_deficit = Math.max(-(C-C_prev)/C_prev, 0);

			for (int i = ccdots.size() - 1; i >= 0; i--) {
				ChargeCarrierDot d = ccdots.get(i);
				d.time -= delta_t;
				if (d.time < 0)
					ccdots.remove(i);
				else {
					double R_tmp = bilinearinterp(e.R, d.x, d.y)*e.e_charge;
					if (d.species == Species.HOLE) {
						if (R_tmp > 0 && frand.next() < R_tmp/bilinearinterp(e.rho_p, d.x, d.y)*delta_t) {
							ccdots.remove(i);
							continue;
						}
					} else {
						if (R_tmp > 0 && frand.next() < -R_tmp/bilinearinterp(e.rho_n, d.x, d.y)*delta_t) {
							ccdots.remove(i);
							continue;
						}
					}

					if (P_deficit > 0 && frand.next() < P_deficit) {
						ccdots.remove(i);
						continue;
					}
				}
			}

			rho_n_dist.generateSamples(N_n, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau, Species.ELECTRON, frand.next())); });
			rho_p_dist.generateSamples(N_p, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau, Species.HOLE, frand.next())); });
			rho_G_dist.generateSamples(N_G, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), Species.ELECTRON, frand.next())); });
			rho_G_dist.generateSamples(N_G, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), Species.HOLE, frand.next())); });
			rho_n_dist.generateSamples(N_n_excess, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), Species.ELECTRON, frand.next())); });
			rho_p_dist.generateSamples(N_p_excess, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), Species.HOLE, frand.next())); });

			C_prev = C;
		}

		if ((VectorView) e.opts.gui_view_vec.getSelectedItem() != VectorView.NONE) {
			setalphaBG(1.0);
			try {
				if (e.graphics_threads.size() == e.n_threads) {
					e.graphics_start_barrier.await();
					e.graphics_end_barrier.await();
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}
	
	void drawText() {
		clearStrings();
		
		Graphics2D g = (Graphics2D) img_back.getGraphics();
		g.setRenderingHint(
		RenderingHints.KEY_TEXT_ANTIALIASING,
		RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		if (e.opts.gui_interface.isSelected())
		{
			for (VoltageProbe p: e.voltageprobes) {
				if (e.ground != null)
					drawStringWithBackgroundAndBorder("V = " + Utils.getSI(p.potential - e.ground.potential, "V", 1e-6), p.x*scalefactor-5, p.y*scalefactor - 12, g);
				else
					drawStringWithBackgroundAndBorder("V = " + Utils.getSI(p.potential, "V", 1e-6), p.x*scalefactor-5, p.y*scalefactor - 12, g);
			}

			if (e.ground != null)
				drawStringWithBackgroundAndBorder("Ground = " + Utils.getSI(e.ground.potential - e.ground.potential, "V", 1e-6), e.ground.x*scalefactor-5, e.ground.y*scalefactor - 12, g);

			for (CurrentProbe p: e.currentprobes) {
				double xa = 0.5*(p.x1+p.x2)*scalefactor;
				double ya = 0.5*(p.y1+p.y2)*scalefactor;

				double dx = p.x2 - p.x1;
				double dy = p.y2 - p.y1;
				double len = Utils.length(dx, dy);
				dx = dx/len;
				dy = dy/len;
				if (Math.abs(dx) > Math.abs(dy))
				{
					dx = -Math.abs(dx);
				} else {
					dy = -2*Math.abs(dy);
				}

				drawStringWithBackgroundAndBorder("I = " + Utils.getSI(p.current*e.depth, "A", 1e-9), (int)(xa-8*dy)-5, (int)(ya+12*dx)+5, g);
			}

			drawStringBackgrounds(g);
			drawStrings(g);
			startNewStringLayer();

			{
				double mx_t = e.controls.mx_index;
				double my_t = e.controls.my_index;
				//double mx_t = (mouseX/(double)scalefactor);
				//double my_t = (mouseY/(double)scalefactor);
				int mi = e.controls.mx_index;
				int mj = e.controls.my_index;

				if (mi < 0)
					mi = 0;
				if (mi >= e.nx)
					mi = e.nx-1;
				if (mj < 0)
					mj = 0;
				if (mj >= e.ny)
					mj = e.ny-1;

				Material mat = e.materials[mi][mj];

				int vspacing = 12;
				int voffset = 1 + e.controls.my;
				int hoffset = 5 + e.controls.mx+15;

				if (e.opts.gui_tooltip.isSelected()) {
					if (voffset + 22*vspacing > e.ny*scalefactor) {
						voffset = voffset - ((voffset + 22*vspacing) - e.ny*scalefactor);
					}
					if (hoffset + 120 > e.ny*scalefactor) {
						hoffset = hoffset - ((hoffset + 120) - e.ny*scalefactor);
					}
				}

				String name = "Material: " + mat.type.name + (mat.modified? " (Modified)" : "");

				this.drawBigStringWithBackground(name, hoffset, voffset + 1*vspacing, g);
				if (e.opts.gui_tooltip.isSelected()) {
					voffset = voffset+3;
					int line = 2;
					drawTwoColumnString("E" , 							Utils.getSI(Utils.length(bilinearinterp(e.Ex, mx_t-0.5,my_t), bilinearinterp(e.Ey, mx_t,my_t-0.5)), "V/m", 1e-6), hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("B" , 							Utils.getSI(e.parity*bilinearinterp(e.Bz, mx_t-0.5, my_t-0.5), "T", 1e-9),	hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03d5" , 						Utils.getSI(bilinearinterp(e.phi,mx_t, my_t), "V", 1e-6),					hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u2130" , 						Utils.getSI(mat.emf, "V/m", 1e-6),										hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03b5/\u03b5\u2080" , 		Utils.getSI(mat.eps_r, ""),										hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03bc/\u03bc\u2080" , 		Utils.getSI(mat.mu_r, ""),											hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u2099" , 				Utils.getSI(bilinearinterp(e.rho_n,mx_t, my_t), "C/m^3", 1e-9),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u209A" , 				Utils.getSI(bilinearinterp(e.rho_p,mx_t, my_t), "C/m^3", 1e-9),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u2080",					Utils.getSI(bilinearinterp(e.rho_back,mx_t, my_t), "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1" , 						Utils.getSI(bilinearinterp(e.rho_free,mx_t, my_t), "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J\u2099" ,						Utils.getSI(Utils.length(bilinearinterp(e.Jx_n,mx_t-0.5,my_t), bilinearinterp(e.Jy_n,mx_t,my_t-0.5)), "A/m^2", 1),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J\u209A" , 					Utils.getSI(Utils.length(bilinearinterp(e.Jx_p,mx_t-0.5,my_t), bilinearinterp(e.Jy_p,mx_t,my_t-0.5)), "A/m^2", 1),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J" , 							Utils.getSI(Utils.length(bilinearinterp(e.Jx_free,mx_t-0.5,my_t), bilinearinterp(e.Jy_free,mx_t,my_t-0.5)), "A/m^2", 1),	hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F\u2099" , 					Utils.getSI(bilinearinterp(e.F_n,mx_t, my_t)/e.q_n+bilinearinterp(e.phi,mx_t, my_t)-e.W_semi/e.eVtoJ, "V", 1e-9),	hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F\u209a" , 					Utils.getSI(bilinearinterp(e.F_p,mx_t, my_t)/e.q_p+bilinearinterp(e.phi,mx_t, my_t)-e.W_semi/e.eVtoJ, "V", 1e-9),	hoffset, voffset + line*vspacing, g); line++;
					//drawTwoColumnString("E\u2099" , 					Utils.getSI(bilinearinterp(E0_n,mx_t, my_t)/eVtoJ, "eV"),	hoffset, voffset + line*vspacing, g); line++;
					//drawTwoColumnString("E\u209a" , 					Utils.getSI(bilinearinterp(E0_p,mx_t, my_t)/eVtoJ, "eV"),	hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F" , 							Utils.getSI(bilinearinterp(e.F,mx_t, my_t), "V", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("CMF\u2099" ,					Utils.getSI(Utils.length(bilinearinterp(e.cmfx_n,mx_t-0.5,my_t), bilinearinterp(e.cmfy_n,mx_t,my_t-0.5))/e.q_n, "V/m", 1e-6),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("CMF\u209A" , 					Utils.getSI(Utils.length(bilinearinterp(e.cmfx_p,mx_t-0.5,my_t), bilinearinterp(e.cmfy_n,mx_t,my_t-0.5))/e.q_p, "V/m", 1e-6),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("x" , 							Utils.getSI(mx_t*e.ds, "m"),								hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("y" , 							Utils.getSI(e.ds*e.ny-(my_t+1)*e.ds, "m"),								hoffset, voffset + line*vspacing, g); line++;
				}
			}

			drawStringBackgrounds(g);
			drawStrings(g);
			startNewStringLayer();

			int vspacing = 13;
			int voffset = 3;
			int hoffset = 5;
			int line = 1;
			drawStringWithBackground("Time: " + Utils.getSI(e.time, "s"), hoffset, voffset + line*vspacing, g); line++;
			drawStringWithBackground("Steps/s: " + Utils.getSI(e.opts.gui_simspeed_2.getValue()/e.simFPStimer.getAverageTime(), ""), hoffset, voffset + line*vspacing, g); line++;
			if (e.opts.gui_paused.isSelected())
			{
				drawStringWithBackground("Paused", hoffset, voffset + line*vspacing, g); line++;
			}
			if (e.sign_violation) {
				drawStringWithBackground("Warning: Numerical instability detected. Please decrease timestep.", hoffset, voffset + line*vspacing, g); line++;
			}
			if (e.controls.debugging) {
				long total = Runtime.getRuntime().totalMemory();
				long used  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
				drawStringWithBackground("Used memory " + Utils.getSI(used, "B"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground("Total memory " + Utils.getSI(total, "B"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(e.t4.getName() + " " + Utils.getSI(e.t4.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(t5.getName() + " " + Utils.getSI(t5.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(e.t6.getName() + " " + Utils.getSI(e.t6.getAverageTime()*e.opts.gui_simspeed_2.getValue(), "s"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(e.t7.getName() + " " + Utils.getSI(e.t7.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(e.t8.getName() + " " + Utils.getSI(e.t8.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(FPStimer.getName() + " " + Utils.getSI(1/FPStimer.getAverageTime(), "Hz"), hoffset, voffset + line*vspacing, g); line++;
				drawStringWithBackground(e.simFPStimer.getName() + " " + Utils.getSI(1/e.simFPStimer.getAverageTime(), "Hz"), hoffset, voffset + line*vspacing, g); line++;
			}

			if (probetexttimer > 0) {
				drawStringWithBackground("Data saved to " + e.datafilename, hoffset, voffset + line*vspacing, g); line++;
				probetexttimer--;
			}

			drawStringBackgrounds(g);
			drawStrings(g);

			Brush brush = (Brush) e.opts.gui_brush.getSelectedItem();
			if (!e.opts.gui_brush_highlight.isSelected()) {
				int r = (int)(scalefactor*e.controls.brushsize/e.ds);
				int brushshape = e.opts.gui_brush_1.getSelectedIndex();
				if (Brush.isBrushShapeImportant(brush))
					if (brushshape == 0) {
						g.setColor(new Color(50, 50, 50));
						g.drawOval(e.controls.mx - r, e.controls.my - r, 2*r, 2*r);
						g.setColor(new Color(200, 200, 200));
						g.drawOval(e.controls.mx - r-1, e.controls.my - r-1, 2*r, 2*r);
					} else {
						g.setColor(new Color(50, 50, 50));
						g.drawRect(e.controls.mx - r, e.controls.my - r, 2*r-2, 2*r-2);
						g.setColor(new Color(200, 200, 200));
						g.drawRect(e.controls.mx - r-1, e.controls.my - r-1, 2*r, 2*r);
					}
			}
		}
	}

	class GraphicsThread extends Thread {

		int n_thread;
		int n_threads;
		Random rand = new Random();

		public GraphicsThread(int n, int n_threads) {
			n_thread = n;
			this.n_threads = n_threads;
			System.out.println("Graphics thread " + n_thread);
		}

		int lower(int n_max) {
			return Math.min((n_thread*n_max)/n_threads, n_max);
		}

		int upper(int n_max) {
			return Math.min(((n_thread+1)*n_max)/n_threads, n_max);
		}


		@Override
		public void run() {
			try {
				while (true) {
					e.graphics_start_barrier.await();

					double arrowlength = 10.0/scalefactor;

					double vectorscalingconstant = 0;

					VectorMode vector_display_mode = (VectorMode)e.opts.gui_view_vec_mode.getSelectedItem();

					int density = 75;

					double randomness = 0;

					if (vector_display_mode == VectorMode.ARROWS) {
						randomness = 0.5;
					} else if (vector_display_mode == VectorMode.LINES) {
						randomness = 0.75;
					}

					double[][] vf_x = null;
					double[][] vf_y = null;
					boolean isCurrent = false;

					switch ((VectorView) e.opts.gui_view_vec.getSelectedItem()) {
					case NONE:
						break;
					case D_FIELD:
						vf_x = e.Dx;
						vf_y = e.Dy;
						break;
					case E_FIELD:
						vf_x = e.Ex;
						vf_y = e.Ey;
						break;
					case ELECTRON_CURRENT:
						vf_x = e.Jx_n;
						vf_y = e.Jy_n;
						isCurrent = true;
						break;
					case HOLE_CURRENT:
						vf_x = e.Jx_p;
						vf_y = e.Jy_p;
						isCurrent = true;
						break;
					case TOTAL_CURRENT:
						vf_x = e.Jx_free;
						vf_y = e.Jy_free;
						isCurrent = true;
						break;
					case POYNTING:
						vf_x = e.Sx;
						vf_y = e.Sy;
						break;
					case EMF:
						vf_x = e.emfx;
						vf_y = e.emfy;
						break;
					}

					if (vf_x != null)
					{
						Vector ctr = new Vector(0,0);
						Vector arrow = new Vector(0,0);
						Vector tip1 = new Vector(0,0);
						Vector tip2 = new Vector(0,0);
						Vector body1 = new Vector(0,0);
						Vector body2 = new Vector(0,0);

						if (vector_display_mode == VectorMode.LINES) {
							vectorscalingconstant = 0.01*Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/5.0)/((VectorView) e.opts.gui_view_vec.getSelectedItem()).scale;
							rand.setSeed(n_thread);
							density = 75;
							for (int i = lower(density); i < upper(density); i++) {
								for (int j = 0; j < density; j++) {

									//double x = (e.nx-1)*(i+0.5)/50;
									//double y = (e.ny-1)*(j+0.5)/50;
									double x = e.nx*(i+randomness*(rand.nextFloat()-0.5))/density;
									double y = e.ny*(j+randomness*(rand.nextFloat()-0.5))/density;
									for (int sign = -1; sign <= 1; sign += 2) {

										double prevx = x;
										double prevy = y;
										double dx = 0;
										double dy = 0;


										int steps = 30;
										for (int k = 0; k < steps; k++) {
											dx = bilinearinterp(vf_x, prevx-0.5, prevy);
											dy = bilinearinterp(vf_y, prevx, prevy-0.5);

											double fieldmagnitude = Math.sqrt(dx*dx+dy*dy);
											double alphaFG = bump(0.5*(1.0-k/(double)(steps-1)), 0.5)*Math.min(1, vectorscalingconstant*fieldmagnitude);

											if (fieldmagnitude != 0) {
												dx /= fieldmagnitude;
												dy /= fieldmagnitude;
											}

											double nextx = prevx + dx*arrowlength*0.25*sign;
											double nexty = prevy + dy*arrowlength*0.25*sign;

											drawLine((int)((prevx+0.5)*scalefactor), (int)((prevy+0.5)*scalefactor), (int)((nextx+0.5)*scalefactor), (int)((nexty+0.5)*scalefactor), k == 0 && sign == 1,
												1f, 1f, 1f, (float) alphaFG, 1f);

											prevx = nextx;
											prevy = nexty;
										}
									}
								}

							}
						} else if (vector_display_mode == VectorMode.ARROWS) {
							vectorscalingconstant = 0.01*Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/5.0)/((VectorView) e.opts.gui_view_vec.getSelectedItem()).scale;
							rand.setSeed(n_thread);
							for (int i = lower(density); i < upper(density); i++) {
								for (int j = 0; j < density; j++) {

									//double x = (e.nx-1)*(i+0.5)/50;
									//double y = (e.ny-1)*(j+0.5)/50;
									double x = e.nx*(i+randomness*(rand.nextFloat()-0.5))/density;
									double y = e.ny*(j+randomness*(rand.nextFloat()-0.5))/density;
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
									double alphaFG = (0.1*Math.sqrt(fieldmagnitude));
									drawLine((int)(body1.x*scalefactor), (int)(body1.y*scalefactor), (int)(body2.x*scalefactor), (int)(body2.y*scalefactor), true,
										(float)fieldmagnitude, (float)fieldmagnitude, (float)fieldmagnitude, (float)alphaFG, 1f);
									drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip1.x*scalefactor), (int)(tip1.y*scalefactor), false,
										(float)fieldmagnitude, (float)fieldmagnitude, (float)fieldmagnitude, (float)alphaFG, 1f);
									drawLine((int)(body2.x*scalefactor), (int)(body2.y*scalefactor), (int)(tip2.x*scalefactor), (int)(tip2.y*scalefactor), false,
										(float)fieldmagnitude, (float)fieldmagnitude, (float)fieldmagnitude, (float)alphaFG, 1f);
								}
							}
						} else if (vector_display_mode == VectorMode.DOTS) {
							boolean paused = e.opts.gui_paused.isSelected();

							vectorscalingconstant = Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/10.0)/((VectorView) e.opts.gui_view_vec.getSelectedItem()).scale;
							int lower = lower(dots.size());
							int upper = upper(dots.size());
							for (int i = lower; i < upper; i++) {
								Dot d = dots.get(i);
								if (d.time <= 0 || d.x < 0 || d.y < 0 || d.x >= e.nx || d.y >= e.ny) {
									d.x = e.nx*rand.nextDouble();
									d.y = e.ny*rand.nextDouble();
									d.lifespan = 100*(1+rand.nextDouble());
									d.time = d.lifespan;
									if (isCurrent && bilinearinterp(e.conducting, d.x, d.y) == 0) {
										d.lifespan = 0;
										d.time = 0;
									}
								}
							}

							setColorFloat(1, 1, 1);
							for (int i = lower; i < upper; i++) {
								Dot d = dots.get(i);
								if (d.time > 0) {
									double dx = 0;
									double dy = 0;
									int steps = 10;

									if (!paused) {
										for (int k = 0; k < steps; k++) {
											dx = bilinearinterp(vf_x,d.x-0.5, d.y)*5e-7*vectorscalingconstant/steps;
											dy = bilinearinterp(vf_y,d.x, d.y-0.5)*5e-7*vectorscalingconstant/steps;
											double maxspeed = 0.5;
											double factor = Math.min(1, maxspeed/Math.sqrt(dx*dx+dy*dy));
											d.x += dx*factor;
											d.y += dy*factor;
										}

										double p = d.time/d.lifespan;
										d.brightness = Math.min(1, Math.max(0.1, 10*Math.sqrt(dx*dx+dy*dy)))*bump(p, 1/3.0);
										d.time -= 1;
									}

									double alphaFG = d.brightness;
									drawRectangle((int)((d.x+0.5)*scalefactor)-1, (int)((d.y+0.5)*scalefactor)-1, 3, 3,
										1f, 1f, 1f, (float)alphaFG, 1f);
								}
							}
						} else if (vector_display_mode == VectorMode.CONTOUR) {
							float scalingconstant = (float) (10.0*Math.pow(10.0, e.opts.gui_brightness.getValue()/10.0)/((ScalarView)e.opts.gui_view.getSelectedItem()).scale);
							double spacing = 0.2/scalingconstant;
							double contourwidth = 1e-7;

							for (int i = lower(scalefactor*e.nx); i < upper(scalefactor*e.nx); i++) {
								for (int j = 0; j < scalefactor*e.ny; j++) {
									double u = bilinearinterp(scalarfield, (double)i/scalefactor - 0.5, (double)j/scalefactor - 0.5)/spacing;
									double v = bilinearinterp(gradscalarfield, (double)i/scalefactor - 0.5, (double)j/scalefactor - 0.5)/spacing;
									double f = ((((u%1)+1.5)%1)/Math.abs(v))/contourwidth;
									if (f < 1) {
										drawPixel(i, j, 1f, 1f, 1f, (float)(2*Math.min(f, 1-f)), 1f);
									}
								}
							}
						} else if (vector_display_mode == VectorMode.SPECIES) {
							boolean paused = e.opts.gui_paused.isSelected();

							vectorscalingconstant = Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/10.0)/((VectorView) e.opts.gui_view_vec.getSelectedItem()).scale;
							int lower = lower(ccdots.size());
							int upper = upper(ccdots.size());

							int steps = 10;
							double dt_dot = delta_t/steps;

							for (int i = lower; i < upper; i++) {
								ChargeCarrierDot d = ccdots.get(i);

								if (d.time > 0) {
									if (!paused) {
										if (d.species == Species.ELECTRON) {
											for (int k = 0; k < steps; k++) {
												double s = dt_dot/(e.ds*bilinearinterp(e.rho_n, d.x, d.y));
												d.x += s*bilinearinterp(e.Jx_n,d.x-0.5, d.y);
												d.y += s*bilinearinterp(e.Jy_n,d.x, d.y-0.5);
											}
										} else if (d.species == Species.HOLE) {
											for (int k = 0; k < steps; k++) {
												double s = dt_dot/(e.ds*bilinearinterp(e.rho_p, d.x, d.y));
												d.x += s*bilinearinterp(e.Jx_p,d.x-0.5, d.y);
												d.y += s*bilinearinterp(e.Jy_p,d.x, d.y-0.5);
											}
										}

										double p = d.time/d.lifespan;
										d.brightness = 1.0*bump(p, 1/3.0);
									}

									double alphaFG = d.brightness;
									if (d.species == Species.ELECTRON/* && -d.random_id*bilinearinterp(e.rho_n, d.x, d.y) < 1*/)
										drawRectangle((int)((d.x+0.5)*scalefactor)-1, (int)((d.y+0.5)*scalefactor)-1, 3, 3,
											0.25f, 0.25f, 1f, (float)alphaFG, 1-(float)alphaFG);
									else if (d.species == Species.HOLE/* && d.random_id*bilinearinterp(e.rho_p, d.x, d.y) < 1*/)
										drawRectangle((int)((d.x+0.5)*scalefactor)-1, (int)((d.y+0.5)*scalefactor)-1, 3, 3,
											1f, 0.25f, 0.25f, (float)alphaFG, 1-(float)alphaFG);

								}
							}
						}
					}
					e.graphics_end_barrier.await();
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}

	public void clearStrings() {
		texts.clear();
	}

	public void drawStringBackgrounds(Graphics g) {
		if (!e.opts.gui_text_bg.isSelected() || !e.opts.gui_interface.isSelected())
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
		if (!e.opts.gui_interface.isSelected())
			return;

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
	
	public double bilinearinterp(double[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-2 && Math.abs(y-Math.round(y)) < 1e-2) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= e.nx) i = e.nx - 1;
			if (j >= e.ny) j = e.ny - 1;
			return array[i][j];
		}

		if (xfloor < 0) {
			xfloor = 0;
			fx = 0.0;
		} else if (xfloor >= e.nx - 1) {
			xfloor = e.nx - 2;
			fx = 1.0;
		}
		if (yfloor < 0) {
			yfloor = 0;
			fy = 0.0;
		} else if (yfloor >= e.ny - 1) {
			yfloor = e.ny - 2;
			fy = 1.0;
		}
		double va = array[xfloor][yfloor]*(1.0-fx) + array[xfloor+1][yfloor]*fx;
		double vb = array[xfloor][yfloor+1]*(1.0-fx) + array[xfloor+1][yfloor+1]*fx;

		return va*(1.0-fy) + vb*fy;
	}

	public double bilinearinterp(int[][] array, double x, double y) {
		int xfloor = (int)Math.floor(x);
		int yfloor = (int)Math.floor(y);
		double fx = x - xfloor;
		double fy = y - yfloor;
		if (Math.abs(x-Math.round(x)) < 1e-2 || Math.abs(y-Math.round(y)) < 1e-2) {
			int i = (int)Math.round(x);
			int j = (int)Math.round(y);
			if (i < 0) i = 0;
			if (j < 0) j = 0;
			if (i >= e.nx) i = e.nx - 1;
			if (j >= e.ny) j = e.ny - 1;
			return array[i][j];
		}

		if (xfloor < 0) {
			xfloor = 0;
			fx = 0.0;
		} else if (xfloor >= e.nx - 1) {
			xfloor = e.nx - 2;
			fx = 1.0;
		}
		if (yfloor < 0) {
			yfloor = 0;
			fy = 0.0;
		} else if (yfloor >= e.ny - 1) {
			yfloor = e.ny - 2;
			fy = 1.0;
		}
		double va = array[xfloor][yfloor]*(1.0-fx) + array[xfloor+1][yfloor]*fx;
		double vb = array[xfloor][yfloor+1]*(1.0-fx) + array[xfloor+1][yfloor+1]*fx;

		return va*(1.0-fy) + vb*fy;
	}
	
	public class Text {
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

	public enum ScalarView {
		NONE("No scalar overlay",															"",				ColorScheme.OTHER,			1),
		E_FIELD("View E field magnitude",													"V/m",			ColorScheme.GREEN,			1),
		B_FIELD("View B field",																"T",			ColorScheme.CYAN_YELLOW,	1e-5),
		CHARGE("View \u03c1: Net charge density",											"C/m^3",		ColorScheme.OTHER,			1e5),
		CURRENT("View J: Total current magnitude",											"A/m^2",		ColorScheme.GREEN,			1e7),
		H_FIELD("View H field",																"A/m",			ColorScheme.CYAN_YELLOW,	1e-5/1.257e-6),
		POTENTIAL("View \u03d5: Electric scalar potential", 								"V",			ColorScheme.RED_BLUE,		1),
		ENERGY("View u: Electromagnetic energy density",									"J/m^3",		ColorScheme.GREEN,			1),
		ELECTRON_CHARGE("View \u03c1\u2099: Electron charge density",						"C/m^3",		ColorScheme.RED_BLUE,		1),
		HOLE_CHARGE("View \u03c1\u209A: Hole charge density",								"C/m^3",		ColorScheme.RED_BLUE,		1),
		COMBINED_CHARGE("View: Combined electron+hole charge density",						"log[C/m^3]",	ColorScheme.OTHER,			1),
		BACKGROUND_CHARGE("View \u03c1\u2080: Background charge density",					"C/m^3",		ColorScheme.RED_BLUE,		1),
		HEAT("View Q: Heat dissipation",													"W/m^3",		ColorScheme.RED_BLUE,		1e12),
		ENTROPY("View s: Entropy generation (Free energy dissipation)",						"J/(m^3 s)",	ColorScheme.RED_BLUE,		1e12),
		ELECTRON_POTENTIAL("View F\u2099: Electron chemical potential (quasi Fermi level)",	"V",			ColorScheme.RED_BLUE, 		1),
		HOLE_POTENTIAL("View F\u209A: Hole chemical potential (quasi Fermi level)",			"V",			ColorScheme.RED_BLUE, 		1),
		AVERAGE_POTENTIAL("View F: Average electrochemical potential",						"V",			ColorScheme.RED_BLUE,		1),
		GENERATION("View G: Carrier generation rate",										"1/(m^3 s)",	ColorScheme.GREEN,			1e31),
		RECOMBINATION("View R: Carrier recombination rate",									"1/(m^3 s)",	ColorScheme.GREEN,			1e31),
		LIGHT("View: Emitted light",														"",				ColorScheme.WHITE,			1e30),
		DEBUG("Debug",																		"",				ColorScheme.RED_BLUE,		1);
	
		enum ColorScheme {
			RED_BLUE, CYAN_YELLOW, GREEN, WHITE, OTHER;
		}
	
		public String name;
		public String unit;
		public ColorScheme colorscheme;
		public double scale; //Typical order of magnitude of the quantity
	
		ScalarView(String name, String unit, ColorScheme colorScheme, double scale)
		{
			this.name = name;
			this.unit = unit;
			this.colorscheme = colorScheme;
			this.scale = scale;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}

	public enum VectorView {
		NONE("No vector overlay",							1),
		E_FIELD("View E field",								1),
		D_FIELD("View D field",								8.85e-12),
		ELECTRON_CURRENT("View J\u2099: Electron current",	1),
		HOLE_CURRENT("View J\u209A: Hole current",			1),
		TOTAL_CURRENT("View J: Total current",				1),
		EMF("View \u2130: External electromotive force",	1),
		POYNTING("View S: Poynting vector",					1);
	
		public String name;
		public double scale;
	
		VectorView(String name, double scale)
		{
			this.name = name;
			this.scale = scale;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}

	public enum VectorMode {
		ARROWS("Show vectors"),
		LINES("Show lines"),
		DOTS("Show dots"),
		CONTOUR("Show scalar contours"),
		SPECIES("Show charge carriers");
	
		String name;
		VectorMode(String name)
		{
			this.name = name;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}
	
	class Dot {
		double x = 0;
		double y = 0;
		double lifespan = 0;
		double time = 0;
		double brightness = 0;
		
		public Dot(double x, double y, double lifespan, double time) {
			this.x = x;
			this.y = y;
			this.lifespan = lifespan;
			this.time = time;
		}
	}

	class ChargeCarrierDot extends Dot {
		Species species;
		double random_id;
		
		public ChargeCarrierDot(double x, double y, double lifespan, double time, Species species, double random_id) {
			super(x, y, lifespan, time);
			this.species = species;
			this.random_id = random_id;
		}
	}
}


class RenderCanvas extends JPanel {

	private static final long serialVersionUID = 7369516276529576171L;

	Simulation parent;
	@Override
	public void paintComponent(Graphics real) {
		//synchronized(parent.renderer) {
			real.drawImage(parent.renderer.img_front, 0, 0, parent.opts);
		//}
	}

	public RenderCanvas(Simulation w) {
		parent = w;
	}
}
