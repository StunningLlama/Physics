// Copyright (c) Brandon Li 2025-2026
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
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import javax.swing.JPanel;

import electrodynamics.Controls.Brush;
import electrodynamics.Simulation.Species;
import electrodynamics.plot.Plot;
import electrodynamics.plot.ProbePlot;
import electrodynamics.probe.ChargeProbe;
import electrodynamics.probe.CurrentProbe;
import electrodynamics.probe.FluxProbe;
import electrodynamics.probe.Probe;
import electrodynamics.probe.VoltageProbe;
import electrodynamics.util.DistributionSampler;
import electrodynamics.util.FastRandom;
import electrodynamics.util.PeriodicTask;
import electrodynamics.util.Timer;
import electrodynamics.util.Utils;
import electrodynamics.util.Vector;

public class Renderer extends PeriodicTask {
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
	public int scalefactor;
	public double scalefactor_real;
	public int imgwidth = 0;
	public int imgheight = 0;
	public int canvas_size = 768;
	public double targetframerate = 60;
	public double frameduration = 1000/targetframerate;
	
	double t_prev = 0;
	double delta_t = 0;

	ArrayList<Dot> dots = new ArrayList<>();
	ArrayList<ChargeCarrierDot> ccdots = new ArrayList<>();
	int numdots = 2000;
	double rho_n_max = 0;
	double rho_p_max = 0;
	double C_prev = 0;

	int probetexttimer = 0;

	/* Multithreading */
	
	CyclicBarrier graphics_start_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	CyclicBarrier graphics_mid_barrier = new CyclicBarrier(SemiSim.n_threads);
	CyclicBarrier graphics_end_barrier = new CyclicBarrier(SemiSim.n_threads + 1);

	/* Electron and hole dots */
	
	DistributionSampler rho_p_dist = new DistributionSampler();
	DistributionSampler rho_n_dist = new DistributionSampler();
	DistributionSampler rho_G_dist = new DistributionSampler();
	FastRandom frand = new FastRandom();
	double cc_default_dot_density = 5e11;		// How many electron/hole dots to draw
	double tau = 1e-12;		// How long a dot stays on the screen

	public int new_canvas_size;
	
	/* Performance profiling */
	
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

		rho_p_dist.init(e.nx, e.ny);
		rho_n_dist.init(e.nx, e.ny);
		rho_G_dist.init(e.nx, e.ny);
		
		setCanvasSize(canvas_size);
		
		resetChargeDots();
	}
	
	public void setCanvasSize(int canvas_size) {
		this.canvas_size = canvas_size;
		
		scalefactor = (int)((double)canvas_size/e.ny);
		scalefactor_real = (double)canvas_size/e.ny;
		if (scalefactor < 1) scalefactor = 1;

		imgwidth = (int)Math.ceil(scalefactor*e.nx);
		imgheight = (int)Math.ceil(scalefactor*e.ny);

		e.canvas.setPreferredSize(new Dimension(canvas_size, canvas_size));
		img_back = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
		img_front = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
		imgData = ((DataBufferInt)img_back.getRaster().getDataBuffer()).getData();
	}
	
	public void resetChargeDots() {
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

		int scansize = scalefactor*e.nx;
		int lines = scalefactor*e.ny;

		if (x >= 0 && y >= 0 && x+w < scansize && y+h < lines) {

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
			drawPixels();
			drawOverlay();
			drawText();
			copyImage(img_back, img_front);
			e.canvas.repaint();

			FPStimer.stop();
			FPStimer.start();
		} catch (Exception e1) {
			SemiSim.displayErrorMessage(e1);
		}
		finally {
			e.rwLock.readLock().unlock();
		}

        SemiSim.instance.threadPool.schedule(this, nextDelay(frameduration), TimeUnit.MILLISECONDS);
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
		boolean showborders = e.opts.menu_borders.isSelected();

		if (e.opts.menu_elem_colors.isSelected()) {
			for (int i = 0; i < e.nx; i++) {
				for (int j = 0; j < e.ny; j++) {
					setColor(e.materials[i][j].type.color_r, e.materials[i][j].type.color_g, e.materials[i][j].type.color_b);
					
					double alpha = 1.0;
					if (showborders && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
						if (e.materials[i+1][j].type != e.materials[i][j].type || e.materials[i][j+1].type != e.materials[i][j].type)
							alpha -= 0.1;
						if (e.materials[i-1][j].type != e.materials[i][j].type || e.materials[i][j-1].type != e.materials[i][j].type)
							alpha += 0.1;
					}

					setalphaFG(alpha);
					
					setPixel(i, j);
				}
			}
		} else {
			for (int i = 0; i < e.nx; i++) {
				for (int j = 0; j < e.ny; j++) {
					setColor(e.materials[i][j].type.color_grayscale, e.materials[i][j].type.color_grayscale, e.materials[i][j].type.color_grayscale);

					double alpha = 1.0;
					if (showborders && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
						if (e.materials[i+1][j].type != e.materials[i][j].type || e.materials[i][j+1].type != e.materials[i][j].type)
							alpha -= 0.1;
						if (e.materials[i-1][j].type != e.materials[i][j].type || e.materials[i][j-1].type != e.materials[i][j].type)
							alpha += 0.1;
					}

					setalphaFG(alpha);

					setPixel(i, j);
				}
			}
		}

		if (e.controls.moving_selection) {
			if (e.opts.menu_elem_colors.isSelected()) {
				for (int i = 0; i < e.nx; i++) {
					for (int j = 0; j < e.ny; j++) {
						int si = i-e.controls.delta_mx;
						int sj = j-e.controls.delta_my;
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection[si][sj].m.type != MaterialType.VACUUM) {
							setColor(e.controls.selection[si][sj].m.type.color_r, e.controls.selection[si][sj].m.type.color_g, e.controls.selection[si][sj].m.type.color_b);
							setPixel(i, j);
						}
					}
				}
			} else {
				for (int i = 0; i < e.nx; i++) {
					for (int j = 0; j < e.ny; j++) {
						int si = i-e.controls.delta_mx;
						int sj = j-e.controls.delta_my;
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection[si][sj].m.type != MaterialType.VACUUM) {
							setColor(e.controls.selection[si][sj].m.type.color_grayscale, e.controls.selection[si][sj].m.type.color_grayscale, e.controls.selection[si][sj].m.type.color_grayscale);
							setPixel(i, j);
						}
					}
				}
			}
		}

		ScalarView scalarview = e.controls.scalarview.getOption();
		ScalarMode scalarmode = e.controls.scalarmode.getOption();

		if (scalarview != ScalarView.NONE && scalarmode != ScalarMode.NONE) {
			double offset = 0;
			
			if (e.ground != null) {
				if (scalarview == ScalarView.ELECTRON_POTENTIAL) {
					offset = e.conducting[e.ground.x][e.ground.y]*(e.F_n[e.ground.x][e.ground.y]/e.q_n+e.phi[e.ground.x][e.ground.y]-e.W_semi/e.eVtoJ);
				} else if (scalarview == ScalarView.HOLE_POTENTIAL) {
					offset = e.conducting[e.ground.x][e.ground.y]*(e.F_p[e.ground.x][e.ground.y]/e.q_p+e.phi[e.ground.x][e.ground.y]-e.W_semi/e.eVtoJ);
				} else if (scalarview == ScalarView.AVERAGE_POTENTIAL) {
					offset = e.F[e.ground.x][e.ground.y];
				}
				
				if (!Double.isFinite(offset))
					offset = 0;
			}

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
						scalarfield[i][j] = e.conducting[i][j]*(e.F_n[i][j]/e.q_n+e.phi[i][j]-e.W_semi/e.eVtoJ) - offset;
						break;
					case HOLE_POTENTIAL:
						scalarfield[i][j] = e.conducting[i][j]*(e.F_p[i][j]/e.q_p+e.phi[i][j]-e.W_semi/e.eVtoJ) - offset;
						break;
					case DEBUG:
						scalarfield[i][j] = e.debug[i][j];
						break;
					case GENERATION:
						scalarfield[i][j] = e.G[i][j];
						break;
					case RECOMBINATION:
						scalarfield[i][j] = e.R[i][j];
						break;
					case AVERAGE_POTENTIAL:
						scalarfield[i][j] = e.F[i][j] - offset;
						break;
					case LIGHT:
						scalarfield[i][j] = -e.materials[i][j].semiconducting*(e.G[i][j]-e.R[i][j]);
						break;
					}
				}
			}
		}

		if (scalarmode == ScalarMode.CONTOUR_COLORS || scalarmode == ScalarMode.CONTOUR) {
			for (int i = 1; i < e.nx-1; i++) {
				for (int j = 1; j < e.ny-1; j++) {
					double gx = scalarfield[i+1][j]-scalarfield[i-1][j];
					double gy = scalarfield[i][j+1]-scalarfield[i][j-1];
					gradscalarfield[i][j] = Utils.length((Double.isFinite(gx))? gx : 0, (Double.isFinite(gy))? gy : 0)/(2*e.ds);
				}
			}
		}


		if (scalarmode == ScalarMode.COLORS || scalarmode == ScalarMode.CONTOUR_COLORS) {
			ScalarView.ColorScheme colorscheme = e.controls.scalarview.getOption().colorscheme;
			float scalingconstant = (float) (10.0*Math.pow(10.0, e.opts.gui_brightness.getValue()/10.0)/scalarview.scale);

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
		}

		boolean showbrush = Brush.isBrushShapeImportant(brush);
		boolean highlight = e.opts.gui_brush_highlight.isSelected();
		for (int i = 0; i < e.nx; i++) {
			for (int j = 0; j < e.ny; j++) {
				MaterialType type = e.materials[i][j].type;
				if (type.isInteractable() && i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1) {
					int ci = 0;
					int cj = 0;
					int shading = 0;
					
					if (type == MaterialType.SWITCH) {
						shading = 2*((i/2+j/2)%2)-1;
					} else {
						int dir = (((int)Math.round(2*e.materials[i][j].emf_direction/Math.PI) % 4) + 4)%4;
						if (dir == 0) {
							ci = i%3;
							cj = (i*2/3+j)%4;
						} else if (dir == 1) {
							ci = j%3;
							cj = (j*2/3+i)%4;
						} else if (dir == 2) {
							ci = (e.nx-i)%3;
							cj = ((e.nx-i)*2/3+j)%4;
						} else if (dir ==  3) {
							ci = (e.nx-j)%3;
							cj = ((e.nx-j)*2/3+i)%4;
						}
						shading = (ci == 0 && cj != 3) || (ci == 1 && cj == 1)? -1 : 1;
					}
					
					//int shading2 = 2*((i+j)%2)-1;

					setalphaBG(0.25);
					setalphaFG(0.75);

					int offset = 10*shading;
					if (e.controls.selected_EMF[i][j])
						offset = 50*shading-30;
					
					if (e.materials[i+1][j].type != type || e.materials[i-1][j].type != type || e.materials[i][j+1].type != type || e.materials[i][j-1].type != type)
					{
						offset = -30;
						if (e.materials[i][j].activated == 0) {
							offset -= 100;
						}
					} else if (e.materials[i][j].activated == 0) {
						offset -= 200;
					}


					int delta_r = type.color_r+offset;
					int delta_g = type.color_g+offset;
					int delta_b = type.color_b+offset;
					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				}

				if (e.controls.selected[i][j] || (showbrush && highlight && e.controls.under_brush[i][j]))
				{
					setalphaBG(0.75);
					setalphaFG(0.25);

					int s = e.controls.selected[i][j]? 1:0;
					int h = (showbrush && highlight && e.controls.under_brush[i][j])? 1:0;
					int delta_r = 256*s + 256*h;
					int delta_g = 100*s + 256*h;
					int delta_b = 256*s + 256*h;

					setColor(delta_r, delta_g, delta_b);

					setPixel(i, j);
				}
				
				if (showbrush && !highlight && e.controls.under_brush[i][j]) {
					if (i > 0 && j > 0 && i < e.nx-1 && j < e.ny-1 && !(e.controls.under_brush[i-1][j] && e.controls.under_brush[i+1][j] && e.controls.under_brush[i][j-1] && e.controls.under_brush[i][j+1]))
					{
						setalphaBG(0.25);
						setalphaFG(0.75);

						if ((image_r[i][j] + image_g[i][j] + image_b[i][j])/3.0 < 0.75)
							setColor(256, 256, 256);
						else
							setColor(50, 50, 50);
						setPixel(i, j);
					}
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


		setalphaBG(1);
		setalphaFG(1);
		setColorFloat(0.7f, 0.7f, 0.7f);

		if (Brush.drawLine(brush) && e.controls.mouse_pressed_prev) {
			drawPixelLine(e.controls.mx_start, e.controls.my_start, e.controls.mx, e.controls.my);
		}

		if (brush == Brush.ZOOM && e.controls.mouse_pressed_prev) {
			int x1 = e.controls.mx_start;
			int y1 = e.controls.my_start;
			int x2 = e.controls.mx;
			int y2 = e.controls.my;
			drawPixelLine(x1, y1, x2, y1);
			drawPixelLine(x2, y1, x2, y2);
			drawPixelLine(x2, y2, x1, y2);
			drawPixelLine(x1, y2, x1, y1);
		}


		if (e.opts.menu_interface.isSelected())
		{
			if (e.opts.menu_probes.isSelected())
			{
				for (Probe p0 : e.probes) {
					if (p0 instanceof VoltageProbe) {
						VoltageProbe p = (VoltageProbe) p0;
						setalphaFG(1.0);
						setColorFloat(0.5f, 1.0f, 1.0f);
						drawPixelRectangle(p.x-1, p.y-1, 3, 3);
						
						setalphaFG(0.1);
						setColorFloat(1.0f, 1.0f, 1.0f);
						drawPixelLine(p.x, p.y, p.labelcoord.x, p.labelcoord.y);
					}
					else if (p0 instanceof CurrentProbe) {
						CurrentProbe p = (CurrentProbe) p0;
						setalphaFG(1.0);
						setColorFloat(0.5f, 1.0f, 1.0f);
						drawPixelRectangle(p.x1-1, p.y1-1, 3, 3);
						drawPixelRectangle(p.x2-1, p.y2-1, 3, 3);

						setalphaFG(0.3);
						setColorFloat(0.5f, 1.0f, 1.0f);
						drawPixelLine(p.x1, p.y1, p.x2, p.y2);
						double dx_perp = (p.y2-p.y1);
						double dy_perp = -(p.x2-p.x1);
						double norm = Utils.length(dx_perp, dy_perp);
						if (norm > 0) {
							dx_perp = dx_perp/norm;
							dy_perp = dy_perp/norm;

							int dist = 3;

							setPixel((int)Math.round(p.x1 + dist*dx_perp), (int)Math.round(p.y1+dist*dy_perp));
							setPixel((int)Math.round(p.x2 + dist*dx_perp), (int)Math.round(p.y2+dist*dy_perp));
						}
						
						setalphaFG(0.1);
						setColorFloat(1.0f, 1.0f, 1.0f);
						drawPixelLine((p.x1 + p.x2)/2, (p.y1+p.y2)/2, p.labelcoord.x, p.labelcoord.y);
					}
					else if (p0 instanceof ChargeProbe) {
						ChargeProbe p = (ChargeProbe) p0;
						setalphaFG(1.0);
						setColorFloat(0.5f, 1.0f, 1.0f);

						setalphaFG(0.3);
						setColorFloat(0.5f, 1.0f, 1.0f);
						drawPixelLine(p.x1, p.y1, p.x2, p.y1);
						drawPixelLine(p.x2, p.y1, p.x2, p.y2);
						drawPixelLine(p.x2, p.y2, p.x1, p.y2);
						drawPixelLine(p.x1, p.y2, p.x1, p.y1);
						
						setalphaFG(0.1);
						setColorFloat(1.0f, 1.0f, 1.0f);
						drawPixelLine((p.x1 + p.x2)/2, (p.y1+p.y2)/2, p.labelcoord.x, p.labelcoord.y);
					} else if (p0 instanceof FluxProbe) {
						FluxProbe p = (FluxProbe) p0;
						setalphaFG(1.0);
						setColorFloat(0.5f, 1.0f, 1.0f);

						setalphaFG(0.3);
						setColorFloat(0.5f, 1.0f, 1.0f);
						drawPixelLine(p.x1, p.y1, p.x2, p.y1);
						drawPixelLine(p.x2, p.y1, p.x2, p.y2);
						drawPixelLine(p.x2, p.y2, p.x1, p.y2);
						drawPixelLine(p.x1, p.y2, p.x1, p.y1);
						
						setalphaFG(0.1);
						setColorFloat(1.0f, 1.0f, 1.0f);
						drawPixelLine((p.x1 + p.x2)/2, (p.y1+p.y2)/2, p.labelcoord.x, p.labelcoord.y);
					}
				}
			}

			setalphaFG(1.0);
			setColorFloat(1.0f, 1.0f, 1.0f);

			for (Plot p: e.plots) {
				if (p.frame.isVisible() && !(p instanceof ProbePlot)) {
					drawPixelRectangle((int)p.x1-1, (int)p.y1-1, 3, 3);
					drawPixelRectangle((int)p.x2-1, (int)p.y2-1, 3, 3);
				}
			}

			setalphaFG(1.0);
			setColorFloat(1.0f, 1.0f, 1.0f);

			for (Plot p: e.plots) {
				if (p != null && p.frame.isVisible() && !(p instanceof ProbePlot)) {
					drawPixelLine((int)p.x1, (int)p.y1, (int)p.x2, (int)p.y2);
				}
			}

		}

		setalphaFG(0.8);
		setColorFloat(1.0f, 1.0f, 1.0f);
		
		if (e.controls.texting) {
			drawPixelLine(e.controls.text_x, e.controls.text_y, e.controls.text_x, e.controls.text_y+7);
		}

		stampPixelData();

	}
	
	void drawOverlay() {
		
		boolean show_carriers = e.opts.gui_carriers.isSelected();
		delta_t = e.time - t_prev;
		t_prev = e.time;

		t5.start();
		if (show_carriers && (!e.opts.gui_paused.isSelected() || delta_t > 0)) {

			double C = cc_default_dot_density*Math.pow(10.0, e.opts.gui_carrier_density.getValue()/20.0);

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
					double R_tmp = Utils.bilinearinterp(e.R, d.x, d.y, e.nx, e.ny)*e.e_charge;
					if (d.species == Species.HOLE) {
						if (R_tmp > 0 && frand.next() < R_tmp/Utils.bilinearinterp(e.rho_p, d.x, d.y, e.nx, e.ny)*delta_t) {
							ccdots.remove(i);
							continue;
						}
					} else {
						if (R_tmp > 0 && frand.next() < -R_tmp/Utils.bilinearinterp(e.rho_n, d.x, d.y, e.nx, e.ny)*delta_t) {
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
		t5.stop();

		setalphaBG(1.0);
		try {
			if (SemiSim.instance.graphics_threads.size() == SemiSim.n_threads) {
				graphics_start_barrier.await();
				graphics_end_barrier.await();
			}
		} catch (InterruptedException | BrokenBarrierException e) {
			e.printStackTrace();
		}
	}
	
	void drawText() {
		clearStrings();

		Graphics2D g = (Graphics2D) img_back.getGraphics();
		g.setRenderingHint(
		RenderingHints.KEY_TEXT_ANTIALIASING,
		RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		if (e.opts.menu_probes.isSelected()) {
			int index = 0;
			for (Probe p: e.probes) {
				if (p == e.ground)
					drawMonospacedString("Ground = " + Utils.getSI_fixedsigfigs(e.ground.potential - e.ground.potential, "V", 1e-6), e.ground.labelcoord.x*scalefactor, e.ground.labelcoord.y*scalefactor, g);
				else if (p instanceof VoltageProbe)
					drawMonospacedString("V" + e.getProbeName(index) + " = " + Utils.getSI_fixedsigfigs(((VoltageProbe)p).potential, "V", 1e-6), p.labelcoord.x*scalefactor, p.labelcoord.y*scalefactor, g);
				else if (p instanceof CurrentProbe)
					drawMonospacedString("I" + e.getProbeName(index) + " = " + Utils.getSI_fixedsigfigs(((CurrentProbe)p).current, "A", 1e-9), p.labelcoord.x*scalefactor, p.labelcoord.y*scalefactor, g);
				else if (p instanceof ChargeProbe)
					drawMonospacedString("Q" + e.getProbeName(index) + " = " + Utils.getSI_fixedsigfigs(((ChargeProbe)p).charge, "C"), p.labelcoord.x*scalefactor, p.labelcoord.y*scalefactor, g);
				else if (p instanceof FluxProbe)
					drawMonospacedString("Φ" + e.getProbeName(index) + " = " + Utils.getSI_fixedsigfigs(((FluxProbe)p).flux, "Wb"), p.labelcoord.x*scalefactor, p.labelcoord.y*scalefactor, g);
				index++;
			}

			drawStrings(g);
			startNewStringLayer();
		}

		{
			int mx = e.controls.mx;
			int my = e.controls.my;

			if (mx < 0)
				mx = 0;
			if (mx >= e.nx)
				mx = e.nx-1;
			if (my < 0)
				my = 0;
			if (my >= e.ny)
				my = e.ny-1;

			Material mat = e.materials[mx][my];

			int vspacing = 12;
			int voffset = 1 + (int)(e.controls.my_screen*scalefactor/scalefactor_real);
			int hoffset = 5 + (int)(e.controls.mx_screen*scalefactor/scalefactor_real)+15;
			//int voffset = 1 + (int)(e.controls.my*scalefactor);
			//int hoffset = 5 + (int)(e.controls.mx*scalefactor)+15;

			if (!e.controls.zoomed) {
				if (e.opts.menu_tooltip.isSelected()) {
					if (voffset + 22*vspacing > e.ny*scalefactor) {
						voffset = voffset - ((voffset + 22*vspacing) - e.ny*scalefactor);
					}
					if (hoffset + 120 > e.ny*scalefactor) {
						hoffset = hoffset - ((hoffset + 120) - e.ny*scalefactor);
					}
				}

				if (e.opts.menu_materialname.isSelected()) {
					String name = "Material: " + mat.type.name + (mat.modified? " (Modified)" : "");
					this.drawBigString(name, hoffset, voffset + 1*vspacing, g);
				}
				
				if (e.opts.menu_tooltip.isSelected()) {
					voffset = voffset+3;
					int line = 2;
					drawTwoColumnString("E" , 							Utils.getSI(Utils.bilinearinterp_length(e.Ex, e.Ey, mx, my, e.nx, e.ny), "V/m", 1e-6), 				hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("B" , 							Utils.getSI(e.parity*Utils.bilinearinterp(e.Bz, mx-0.5, my-0.5, e.nx, e.ny), "T", 1e-9),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03d5" , 						Utils.getSI(e.phi[mx][my], "V", 1e-6),				hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u2130" , 						Utils.getSI(mat.emf, "V/m", 1e-6),					hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03b5/\u03b5\u2080" , 		Utils.getSI(mat.eps_r, ""),							hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03bc/\u03bc\u2080" , 		Utils.getSI(mat.mu_r, ""),							hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u2099" , 				Utils.getSI(e.rho_n[mx][my], "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u209A" , 				Utils.getSI(e.rho_p[mx][my], "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1\u2080",					Utils.getSI(e.rho_back[mx][my], "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("\u03c1" , 						Utils.getSI(e.rho_free[mx][my], "C/m^3", 1e-9),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J\u2099" ,						Utils.getSI(Utils.bilinearinterp_length(e.Jx_n, e.Jy_n, mx, my, e.nx, e.ny), "A/m^2", 1),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J\u209A" , 					Utils.getSI(Utils.bilinearinterp_length(e.Jx_p, e.Jy_p, mx, my, e.nx, e.ny), "A/m^2", 1),			hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("J" , 							Utils.getSI(Utils.bilinearinterp_length(e.Jx_free, e.Jy_free, mx, my, e.nx, e.ny), "A/m^2", 1),		hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F\u2099" , 					Utils.getSI(e.F_n[mx][my]/e.q_n + e.phi[mx][my] - e.W_semi/e.eVtoJ, "V", 1e-9),						hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F\u209a" , 					Utils.getSI(e.F_p[mx][my]/e.q_p + e.phi[mx][my] - e.W_semi/e.eVtoJ, "V", 1e-9),						hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("F" , 							Utils.getSI(e.F[mx][my], "V", 1e-9),				hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("x" , 							Utils.getSI(mx*e.ds, "m"),							hoffset, voffset + line*vspacing, g); line++;
					drawTwoColumnString("y" , 							Utils.getSI(e.ds*e.ny-(my+1)*e.ds, "m"),			hoffset, voffset + line*vspacing, g); line++;
				}
			}
			
			drawStrings(g);
			startNewStringLayer();
		}

		int vspacing = 13;
		int voffset = 3;
		int hoffset = 5;
		int line = 1;
		
		if (e.opts.menu_time.isSelected()) {
			drawString("Time: " + Utils.getSI(e.time, "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString("Steps/s: " + Utils.getSI(e.opts.gui_simspeed_2.getValue()/e.simFPStimer.getAverageTime(), ""), hoffset, voffset + line*vspacing, g); line++;

			String sv_a = e.controls.scalarmode.getOption() != ScalarMode.NONE? (e.controls.scalarmode.getOption().shorthand + ": " + e.controls.scalarview.getOption().shorthand) : "";
			String sv_b = e.controls.vectormode.getOption() != VectorMode.NONE? (e.controls.vectormode.getOption().shorthand + ": " + e.controls.vectorview.getOption().shorthand) : "";
			String sv = Arrays.asList(sv_a, sv_b).stream().filter(s -> s != null && !s.isEmpty()).collect(Collectors.joining(" / "));
			drawString(sv, hoffset, voffset + line*vspacing, g); line++;
			if (e.opts.gui_paused.isSelected())
			{
				drawString("Paused", hoffset, voffset + line*vspacing, g); line++;
			}
			if (e.sign_violation) {
				drawString("Warning: Numerical instability detected. Please decrease timestep.", hoffset, voffset + line*vspacing, g); line++;
			}
		}
		
		if (e.controls.debugging) {
			long total = Runtime.getRuntime().totalMemory();
			long used  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			drawString("Used memory " + Utils.getSI(used, "B"), hoffset, voffset + line*vspacing, g); line++;
			drawString("Total memory " + Utils.getSI(total, "B"), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t4.getName() + " " + Utils.getSI(e.t4.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString(t5.getName() + " " + Utils.getSI(t5.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t6.getName() + " " + Utils.getSI(e.t6.getAverageTime()*e.opts.gui_simspeed_2.getValue(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t7.getName() + " " + Utils.getSI(e.t7.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t8.getName() + " " + Utils.getSI(e.t8.getAverageTime(), "s"), hoffset, voffset + line*vspacing, g); line++;
			drawString(FPStimer.getName() + " " + Utils.getSI(1/FPStimer.getAverageTime(), "Hz"), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.simFPStimer.getName() + " " + Utils.getSI(1/e.simFPStimer.getAverageTime(), "Hz"), hoffset, voffset + line*vspacing, g); line++;
		}
		if (e.numerical_overflow) {
			line++;
			drawError("Error: Numerical overflow detected. Please reset simulation.", hoffset, voffset + line*vspacing, g); line += 2;
		}

		if (probetexttimer > 0) {
			drawString("Data saved to " + e.datafilename, hoffset, voffset + line*vspacing, g); line++;
			probetexttimer--;
		}

		drawStrings(g);
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
					graphics_start_barrier.await();

					VectorMode vector_display_mode = e.controls.vectormode.getOption();
					ScalarMode scalar_display_mode = e.controls.scalarmode.getOption();
					ScalarView scalar_view = e.controls.scalarview.getOption();
					VectorView vector_view = e.controls.vectorview.getOption();

					if (vector_display_mode != VectorMode.NONE && vector_view != VectorView.NONE) {
						double arrowlength = 10.0/scalefactor;

						double vectorscalingconstant = 0;


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

						switch (vector_view) {
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

						Vector ctr = new Vector(0,0);
						Vector arrow = new Vector(0,0);
						Vector tip1 = new Vector(0,0);
						Vector tip2 = new Vector(0,0);
						Vector body1 = new Vector(0,0);
						Vector body2 = new Vector(0,0);

						if (vector_display_mode == VectorMode.LINES) {
							vectorscalingconstant = 0.01*Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/5.0)/e.controls.vectorview.getOption().scale;
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
											dx = Utils.bilinearinterp(vf_x, prevx-0.5, prevy, e.nx, e.ny);
											dy = Utils.bilinearinterp(vf_y, prevx, prevy-0.5, e.nx, e.ny);

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
							vectorscalingconstant = 0.01*Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/5.0)/e.controls.vectorview.getOption().scale;
							rand.setSeed(n_thread);
							for (int i = lower(density); i < upper(density); i++) {
								for (int j = 0; j < density; j++) {

									//double x = (e.nx-1)*(i+0.5)/50;
									//double y = (e.ny-1)*(j+0.5)/50;
									double x = e.nx*(i+randomness*(rand.nextFloat()-0.5))/density;
									double y = e.ny*(j+randomness*(rand.nextFloat()-0.5))/density;
									ctr.x = x+0.5;
									ctr.y = y+0.5;

									arrow.x = Utils.bilinearinterp(vf_x,x-0.5, y, e.nx, e.ny);
									arrow.y = Utils.bilinearinterp(vf_y,x, y-0.5, e.nx, e.ny);

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

							vectorscalingconstant = Math.pow(10.0, e.opts.gui_brightness_vec.getValue()/10.0)/e.controls.vectorview.getOption().scale;
							int lower = lower(dots.size());
							int upper = upper(dots.size());
							for (int i = lower; i < upper; i++) {
								Dot d = dots.get(i);
								if (d.time <= 0 || d.x < 0 || d.y < 0 || d.x >= e.nx || d.y >= e.ny) {
									d.x = e.nx*rand.nextDouble();
									d.y = e.ny*rand.nextDouble();
									d.lifespan = 100*(1+rand.nextDouble());
									d.time = d.lifespan;
									if (isCurrent && Utils.bilinearinterp(e.conducting, d.x, d.y, e.nx, e.ny) == 0) {
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
											dx = Utils.bilinearinterp(vf_x,d.x-0.5, d.y, e.nx, e.ny)*5e-7*vectorscalingconstant/steps;
											dy = Utils.bilinearinterp(vf_y,d.x, d.y-0.5, e.nx, e.ny)*5e-7*vectorscalingconstant/steps;
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
						}
					}
					
					graphics_mid_barrier.await();

					if (scalar_view != ScalarView.NONE && (scalar_display_mode == ScalarMode.CONTOUR_COLORS || scalar_display_mode == ScalarMode.CONTOUR)) {
						float scalingconstant = (float) (10.0*Math.pow(10.0, e.opts.gui_brightness.getValue()/10.0)/e.controls.scalarview.getOption().scale);
						double spacing = 0.2/scalingconstant;
						double contourwidth = 1e-7;

						for (int i = lower(scalefactor*e.nx); i < upper(scalefactor*e.nx); i++) {
							for (int j = 0; j < scalefactor*e.ny; j++) {
								double u = Utils.bilinearinterp(scalarfield, (double)i/scalefactor - 0.5, (double)j/scalefactor - 0.5, e.nx, e.ny)/spacing;
								double v = Utils.bilinearinterp(gradscalarfield, (double)i/scalefactor - 0.5, (double)j/scalefactor - 0.5, e.nx, e.ny)/spacing;
								double f = ((((u%1)+1.5)%1)/Math.abs(v))/contourwidth;
								if (f < 1) {
									drawPixel(i, j, 1f, 1f, 1f, (float)(2*Math.min(f, 1-f)), 1f);
								}
							}
						}
					}

					graphics_mid_barrier.await();

					boolean show_carriers = e.opts.gui_carriers.isSelected();
					if (show_carriers) {

						boolean fast = e.opts.menu_carriers_metal.isSelected();
						
						int lower = lower(ccdots.size());
						int upper = upper(ccdots.size());

						int steps = fast? 2 : 10;
						double dt_dot = delta_t/steps;

						for (int i = lower; i < upper; i++) {
							ChargeCarrierDot d = ccdots.get(i);

							if (d.time > 0) {
								
								boolean dorender = !fast || Utils.bilinearinterp(e.semiconducting, d.x, d.y, e.nx, e.ny) > 0;
								
								if (delta_t > 0) {
									if (d.species == Species.ELECTRON) {
										for (int k = 0; k < steps; k++) {
											double s = dt_dot/(e.ds*Utils.bilinearinterp(e.rho_n, d.x, d.y, e.nx, e.ny));
											d.x += s*Utils.bilinearinterp(e.Jx_n,d.x-0.5, d.y, e.nx, e.ny);
											d.y += s*Utils.bilinearinterp(e.Jy_n,d.x, d.y-0.5, e.nx, e.ny);
										}
									} else if (d.species == Species.HOLE) {
										for (int k = 0; k < steps; k++) {
											double s = dt_dot/(e.ds*Utils.bilinearinterp(e.rho_p, d.x, d.y, e.nx, e.ny));
											d.x += s*Utils.bilinearinterp(e.Jx_p,d.x-0.5, d.y, e.nx, e.ny);
											d.y += s*Utils.bilinearinterp(e.Jy_p,d.x, d.y-0.5, e.nx, e.ny);
										}
									}

									double p = d.time/d.lifespan;
									d.brightness = 1.0*bump(p, 1/3.0);
								}

								double alphaFG = d.brightness;
								if (dorender && d.species == Species.ELECTRON/* && -d.random_id*bilinearinterp(e.rho_n, d.x, d.y) < 1*/)
									drawRectangle((int)((d.x+0.5)*scalefactor)-1, (int)((d.y+0.5)*scalefactor)-1, 3, 3,
									0.25f, 0.25f, 1f, (float)alphaFG, 1-(float)alphaFG);
								else if (dorender && d.species == Species.HOLE/* && d.random_id*bilinearinterp(e.rho_p, d.x, d.y) < 1*/)
									drawRectangle((int)((d.x+0.5)*scalefactor)-1, (int)((d.y+0.5)*scalefactor)-1, 3, 3,
									1f, 0.25f, 0.25f, (float)alphaFG, 1-(float)alphaFG);

							}
						}
					}
					graphics_end_barrier.await();
				}
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}
		}
	}

	public void clearStrings() {
		texts.clear();
	}

	public void drawStrings(Graphics g) {

		boolean dodraw = e.opts.menu_text_bg.isSelected() && e.opts.menu_interface.isSelected();

		for (Text text : texts) {
			text.computeDimensions(g);
		}
		
		for (Text text : texts) {
			if (text.hasBackground && (dodraw || text.isError)) {
				g.setFont(text.getFont());
				int x = text.x-3;
				int y = text.y-text.height+6;

				g.setColor(Color.GRAY);
				g.fillRect(x-2, y-2, text.width+4, text.height+4);
			}
		}

		for (Text text : texts) {
			if (text.hasBackground && (dodraw || text.isError)) {
				g.setFont(text.getFont());

				int x = text.x-3;
				int y = text.y-text.height+6;

				g.setColor(Color.BLACK);
				if (text.isError)
					g.setColor(Color.RED);
				
				g.fillRect(x, y, text.width, text.height);
			}
		}
		
		dodraw = e.opts.menu_interface.isSelected();

		for (Text text : texts) {
			if (dodraw || text.isError) {
				g.setFont(text.getFont());

				g.setColor(Color.DARK_GRAY);
				g.drawString(text.text, text.x+1, text.y+1);
				g.setColor(Color.WHITE);
				g.drawString(text.text, text.x, text.y);
			}
		}
	}

	public void startNewStringLayer() {
		texts.clear();
	}

	public void drawString(String str1, int x, int y, Graphics g) {
		Text t = new Text(str1, x, y);
		texts.add(t);
	}
	
	public void drawMonospacedString(String str1, int x, int y, Graphics g) {
		Text t = new Text(str1, x, y);
		t.monospaced = true;
		texts.add(t);
	}
	
	public void drawError(String str1, int x, int y, Graphics g) {
		Text t = new Text(str1, x, y);
		t.isError = true;
		texts.add(t);
	}

	public void drawBigString(String str1, int x, int y, Graphics g) {
		Text t = new Text(str1, x, y);
		t.big = true;
		texts.add(t);
	}

	public void drawTwoColumnString(String str1, String str2, int x, int y, Graphics g) {
		drawString(String.format("%-10s", str1), x, y, g);
		drawString(str2, x+40, y, g);
		texts.get(texts.size()-1).minwidth = 80;
	}
	
	public static class Text {

		public static Font bigfont = new Font(Font.SANS_SERIF, Font.PLAIN, 15);
		public static Font regularfont = new Font(Font.SANS_SERIF, Font.PLAIN, 12);
		public static Font monospacefont = getMonospacedFont();
		
		String text;
		int x;
		int y;
		int minwidth;
		int width = 0;
		int height = 0;
		
		boolean big = false;
		boolean hasBackground = true;
		boolean monospaced = false;
		boolean isError = false;

		public Text(String text, int x, int y) {
			this.text = text;
			this.x = x;
			this.y = y;
			minwidth = 0;
		}
		
		public Font getFont() {
			if (big)
				return bigfont;
			else if (monospaced)
				return monospacefont;
			else
				return regularfont;
		}
		
		public void computeDimensions(Graphics g) {
			g.setFont(getFont());
			width = Math.max(minwidth, g.getFontMetrics().stringWidth(text)+8);
			height = g.getFontMetrics().getHeight()+4;
		}
		
		public static Font getMonospacedFont() {
			Font f = Font.decode("Consolas-PLAIN-12");
			if (f.getFamily() == "Dialog")
				f = Font.decode("Andale Mono-PLAIN-12");
			if (f.getFamily() == "Dialog")
				f = new Font(Font.MONOSPACED, Font.PLAIN, 12);
			return f;
		}
	}

	public enum ScalarView {
		NONE("No scalar overlay",															"None",		"",				ColorScheme.OTHER,			1),
		E_FIELD("View: E field magnitude",													"E",		"V/m",			ColorScheme.GREEN,			1e5),
		B_FIELD("View: B field",															"B",		"T",			ColorScheme.CYAN_YELLOW,	1e-5),
		CHARGE("View \u03c1: Net charge density",											"\u03c1",	"C/m^3",		ColorScheme.OTHER,			1e1),
		CURRENT("View J: Total current magnitude",											"J",		"A/m^2",		ColorScheme.GREEN,			1e7),
		H_FIELD("View H field",																"H",		"A/m",			ColorScheme.CYAN_YELLOW,	1e-5/1.257e-6),
		POTENTIAL("View \u03d5: Electric scalar potential", 								"\u03d5",	"V",			ColorScheme.RED_BLUE,		1),
		ENERGY("View u: Electromagnetic energy density",									"u",		"J/m^3",		ColorScheme.GREEN,			1),
		ELECTRON_CHARGE("View \u03c1\u2099: Electron charge density",						"\u03c1\u2099",	"C/m^3",	ColorScheme.RED_BLUE,		1e1),
		HOLE_CHARGE("View \u03c1\u209A: Hole charge density",								"\u03c1\u209A",	"C/m^3",	ColorScheme.RED_BLUE,		1e1),
		COMBINED_CHARGE("View: Combined electron+hole charge density",						"\u03c1\u2099 and \u03c1\u209A",	"log[C/m^3]",	ColorScheme.OTHER,			1),
		BACKGROUND_CHARGE("View \u03c1\u2080: Static charge density (doping)",				"\u03c1\u2080",	"C/m^3",	ColorScheme.RED_BLUE,		1e1),
		HEAT("View Q: Heat dissipation",													"q",		"W/m^3",		ColorScheme.RED_BLUE,		1e12),
		ENTROPY("View s: Entropy generation (Free energy dissipation)",						"s",		"J/(m^3 s)",	ColorScheme.RED_BLUE,		1e12),
		ELECTRON_POTENTIAL("View F\u2099: Electron chemical potential (quasi Fermi level)",	"F\u2099",	"V",			ColorScheme.RED_BLUE, 		1),
		HOLE_POTENTIAL("View F\u209A: Hole chemical potential (quasi Fermi level)",			"F\u209A",	"V",			ColorScheme.RED_BLUE, 		1),
		AVERAGE_POTENTIAL("View F: Average electrochemical potential",						"F",		"V",			ColorScheme.RED_BLUE,		1),
		GENERATION("View G: Carrier generation rate",										"G",		"1/(m^3 s)",	ColorScheme.GREEN,			1e31),
		RECOMBINATION("View R: Carrier recombination rate",									"R",		"1/(m^3 s)",	ColorScheme.GREEN,			1e31),
		LIGHT("View: Emitted light",														"Light",	"",				ColorScheme.WHITE,			1e30),
		DEBUG("Debug",																		"Debug",	"",				ColorScheme.RED_BLUE,		1e-3);
	
		enum ColorScheme {
			RED_BLUE, CYAN_YELLOW, GREEN, WHITE, OTHER;
		}
	
		public String name;
		public String unit;
		public String shorthand;
		public ColorScheme colorscheme;
		public double scale; //Typical order of magnitude of the quantity
	
		ScalarView(String name, String shorthand, String unit, ColorScheme colorScheme, double scale)
		{
			this.name = name;
			this.shorthand = shorthand;
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
		NONE("No vector overlay",							"None",		1),
		E_FIELD("View E: Electric field",					"E",		10),
		D_FIELD("View D: Displacement field",				"D",		10*8.85e-12),
		ELECTRON_CURRENT("View J\u2099: Electron current",	"J\u2099",	1),
		HOLE_CURRENT("View J\u209A: Hole current",			"J\u209A",	1),
		TOTAL_CURRENT("View J: Total current",				"J",		1),
		EMF("View \u2130: External electromotive force",	"\u2130",	1),
		POYNTING("View S: Poynting vector",					"S",		1);
	
		public String name;
		public String shorthand;
		public double scale;
	
		VectorView(String name, String shorthand, double scale)
		{
			this.name = name;
			this.shorthand = shorthand;
			this.scale = scale;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}

	public enum ScalarMode {
		NONE("Turn off scalar overlay", ""),
		COLORS("Show colors", "color"),
		CONTOUR("Show contours", "contours"),
		CONTOUR_COLORS("Show contours and colors", "color+contours");
	
		String name;
		public String shorthand;
		ScalarMode(String name, String shorthand)
		{
			this.name = name;
			this.shorthand = shorthand;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}
	
	public enum VectorMode {
		NONE("Turn off vector overlay", ""),
		ARROWS("Show vectors", "vectors"),
		LINES("Show lines", "lines"),
		DOTS("Show moving dots", "dots");
	
		String name;
		public String shorthand;
		VectorMode(String name, String shorthand)
		{
			this.name = name;
			this.shorthand = shorthand;
		}
	
		@Override
		public String toString() {
			return name;
		}
	}
	
	public class Dot {
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

	public class ChargeCarrierDot extends Dot {
		Species species;
		double random_id;
		
		public ChargeCarrierDot(double x, double y, double lifespan, double time, Species species, double random_id) {
			super(x, y, lifespan, time);
			this.species = species;
			this.random_id = random_id;
		}
	}
	
	public class RenderCanvas extends JPanel {

		private static final long serialVersionUID = 7369516276529576171L;
		public int zoom_bound_x = 0;
		public int zoom_bound_y = 0;
		public int offset_x = 0;
		public int offset_y = 0;

		Simulation e;
		@Override
		public void paintComponent(Graphics real) {
			Graphics2D g = ((Graphics2D)real);
			g.setBackground(Color.BLACK);
			g.clearRect(0, 0, g.getClipBounds().width, g.getClipBounds().height);
			
			if (!e.controls.zoomed) {
				boolean need_to_rescale = (e.renderer.scalefactor*e.ny != e.renderer.canvas_size);
				if (need_to_rescale)
					g.drawImage(e.renderer.img_front, 0, 0, e.renderer.canvas_size, e.renderer.canvas_size, e.opts);
				else
					g.drawImage(e.renderer.img_front, 0, 0, e.opts);
			} else {
				int smin = Math.min(e.controls.zoom_i2 - e.controls.zoom_i1 + 1, e.controls.zoom_j2 - e.controls.zoom_j1 + 1);
				int smax = Math.max(e.controls.zoom_i2 - e.controls.zoom_i1 + 1, e.controls.zoom_j2 - e.controls.zoom_j1 + 1);
				int dim2 = (int) e.renderer.canvas_size*smin/smax;
				
				if (e.controls.zoom_i2 - e.controls.zoom_i1 < e.controls.zoom_j2 - e.controls.zoom_j1) {
					zoom_bound_x = dim2 - 1;
					zoom_bound_y = e.renderer.canvas_size - 1;
					offset_x = (e.renderer.canvas_size - dim2)/2;
					offset_y = 0;
				} else {
					zoom_bound_x = e.renderer.canvas_size - 1;
					zoom_bound_y = dim2 - 1;
					offset_x = 0;
					offset_y = (e.renderer.canvas_size - dim2)/2;
				}
				g.drawImage(e.renderer.img_front, offset_x, offset_y, zoom_bound_x+1+offset_x, zoom_bound_y+1+offset_y, 
					e.controls.zoom_i1*e.renderer.scalefactor, e.controls.zoom_j1*e.renderer.scalefactor, (e.controls.zoom_i2+1)*e.renderer.scalefactor, (e.controls.zoom_j2+1)*e.renderer.scalefactor, e.opts);
			}
		}

		public RenderCanvas(Simulation w) {
			e = w;
		}
	}
}
