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
import electrodynamics.plot.Plot;
import electrodynamics.probe.Probe;
import electrodynamics.units.Quantity;
import electrodynamics.util.DistributionSampler;
import electrodynamics.util.FastList;
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
	public int[] depth_buf;
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
	public double targetframerate = 60;
	public double frameduration = 1000/targetframerate;
	
	double t_prev = 0;
	double delta_t = 0;

	ArrayList<Dot> dots = new ArrayList<>();
	FastList<ChargeCarrierDot> ccdots = new FastList<>();
	int numdots = 2000;
	double rho_n_max = 0;
	double rho_p_max = 0;
	double C_prev = 0;
	public int carrier_diffusion_warning_timer = 0;

	int probetexttimer = 0;
	boolean show_carriers;
	VectorMode vector_display_mode;
	ScalarMode scalar_display_mode;
	ScalarView scalar_view;
	VectorView vector_view;
	
	float scalingconstant;
	float scalar_offset;

	/* Multithreading */
	
	CyclicBarrier graphics_start_barrier = new CyclicBarrier(SemiSim.n_threads + 1);
	CyclicBarrier graphics_mid_barrier = new CyclicBarrier(SemiSim.n_threads);
	CyclicBarrier graphics_end_barrier = new CyclicBarrier(SemiSim.n_threads + 1);

	/* Electron and hole dots */
	
	DistributionSampler rho_p_dist = new DistributionSampler();
	DistributionSampler rho_n_dist = new DistributionSampler();
	DistributionSampler rho_G_dist = new DistributionSampler();
	FastRandom frand = new FastRandom();
	double cc_default_dot_density;	// How many electron/hole dots to draw, in (C/m)^-1
	double tau;						// How long a dot stays on the screen
	double tau_events;				// How long a flash stays
	
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
		
		setCanvasSize();
		
		resetChargeDots();
	}
	
	public void setCanvasSize() {
		int min_canvas_size = Math.min(e.canvas.getWidth(), e.canvas.getHeight());
		
		scalefactor = (int)((double)min_canvas_size/e.ny);
		scalefactor_real = (double)min_canvas_size/e.ny;
		if (scalefactor < 1) scalefactor = 1;

		int imgwidth_new = (int)Math.ceil(scalefactor*e.nx);	
		int imgheight_new = (int)Math.ceil(scalefactor*e.ny);
		
		if (imgwidth_new != imgwidth || imgheight_new != imgheight) {
			imgwidth = imgwidth_new;
			imgheight = imgheight_new;
			img_back = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
			img_front = (BufferedImage) e.opts.createImage(imgwidth, imgheight);
			imgData = ((DataBufferInt)img_back.getRaster().getDataBuffer()).getData();
			depth_buf = new int[imgData.length];
		}
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
	
	// Wavelength to RGB
	double[] lambda_r = {0.23822, 0.24973, 0.26128, 0.26553, 0.26476, 0.2536, 0.23768, 0.22023, 0.20389, 0.18888, 0.169, 0.12949, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.19944, 0.4376, 0.58059, 0.68286, 0.76516, 0.82701, 0.88226, 0.95151, 1.0112, 1.0442, 1.0729, 1.1021, 1.1058, 1.086, 1.0531, 1.0138, 0.95956, 0.89531, 0.83074, 0.76591, 0.69323};
	double[] lambda_g = {0.19971, 0.17932, 0.14919, 0.12037, 0.093617, 0.09458, 0.10426, 0.11275, 0.11462, 0.11727, 0.12455, 0.14908, 0.1933, 0.25787, 0.32737, 0.38708, 0.44274, 0.49603, 0.54306, 0.59225, 0.64252, 0.69595, 0.74892, 0.79799, 0.84018, 0.86705, 0.88432, 0.89322, 0.89768, 0.89173, 0.87677, 0.86048, 0.83734, 0.81313, 0.78437, 0.75141, 0.70945, 0.65984, 0.60119, 0.54323, 0.47406, 0.38514, 0.29591, 0.21448, 0.12878, 0, 0, 0, 0, 0, 0.064417};
	double[] lambda_b = {0.36184, 0.45077, 0.54787, 0.62678, 0.69213, 0.72851, 0.7625, 0.8006, 0.84243, 0.87363, 0.90106, 0.91545, 0.92739, 0.92995, 0.90334, 0.85833, 0.79228, 0.6955, 0.58576, 0.4834, 0.39609, 0.31148, 0.22847, 0.14139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.062112, 0.10835, 0.13623, 0.15637, 0.17185, 0.18456, 0.19388, 0.20082, 0.20595, 0.20985};

	public void setColorWavelength(double lambda, double intensity) {
		col_r = (float)(interp(lambda_r, (lambda-400)/5)*intensity);
		col_g = (float)(interp(lambda_g, (lambda-400)/5)*intensity);
		col_b = (float)(interp(lambda_b, (lambda-400)/5)*intensity);
	}
	
	public double interp(double[] y, double x) {
		int i = (int) x;
		double f = x - i;
		if (i < 0 || x < 0)
			return y[0];
		if (i >= y.length-1)
			return y[y.length-1];
		return (1-f)*y[i] + f*y[i+1];
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
		Arrays.fill(depth_buf, 0);
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
	
	public void drawRectangle(int x, int y, int w, int h, float col_r, float col_g, float col_b, float alphaFG, float alphaBG, int depth) {

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
					if (depth >= depth_buf[i+j*scansize]) {
						imgData[i+j*scansize] = 255<<24 | r << 16 | g << 8 | b;
						depth_buf[i+j*scansize] = depth;
					}
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
			t5.start();
			drawPixels();
			drawOverlay();
			//drawText();
			copyImage(img_back, img_front);
			e.canvas.repaint();
			t5.stop();

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
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection.mat[si][sj].m.type != MaterialType.VACUUM) {
							setColor(e.controls.selection.mat[si][sj].m.type.color_r, e.controls.selection.mat[si][sj].m.type.color_g, e.controls.selection.mat[si][sj].m.type.color_b);
							setPixel(i, j);
						}
					}
				}
			} else {
				for (int i = 0; i < e.nx; i++) {
					for (int j = 0; j < e.ny; j++) {
						int si = i-e.controls.delta_mx;
						int sj = j-e.controls.delta_my;
						if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && e.controls.selection.mat[si][sj].m.type != MaterialType.VACUUM) {
							setColor(e.controls.selection.mat[si][sj].m.type.color_grayscale, e.controls.selection.mat[si][sj].m.type.color_grayscale, e.controls.selection.mat[si][sj].m.type.color_grayscale);
							setPixel(i, j);
						}
					}
				}
			}
		}

		ScalarView scalarview = e.controls.scalarview.getOption();
		ScalarMode scalarmode = e.controls.scalarmode.getOption();

		scalingconstant = (float) (10.0*Math.pow(10.0, e.opts.gui_brightness.getValue()/10.0)/scalarview.scale);
		scalar_offset = 0;

		if (scalarview != ScalarView.NONE && scalarmode != ScalarMode.NONE) {
			/*double offset = 0;
			
			if (e.hasGround() && e.prefs.chkbox_potential.isSelected()) {
				Ground ground = e.getGround();
				//TODO
				if (scalarview == ScalarView.ELECTRON_POTENTIAL) {
					offset = e.conducting[ground.x][ground.y]*(e.mu_n[ground.x][ground.y]/e.q_n-e.global_voltage_offset);
				} else if (scalarview == ScalarView.HOLE_POTENTIAL) {
					offset = e.conducting[ground.x][ground.y]*(e.mu_p[ground.x][ground.y]/e.q_p-e.global_voltage_offset);
				} else if (scalarview == ScalarView.AVERAGE_POTENTIAL) {
					offset = e.V_avg[ground.x][ground.y];
				}
				
				if (!Double.isFinite(offset))
					offset = 0;
			}*/
			

			switch (scalarview) {
			case ELECTRON_POTENTIAL:
			case HOLE_POTENTIAL:
				scalar_offset = -(float) e.global_voltage_offset;
				break;
			default:
				break;
			}

			setalphaBG(1.0);
			setalphaFG(1.0);
			
			e.computeScalarField(scalarfield, 0, 0, scalarview);

			ScalarView.ColorScheme colorscheme = scalarview.colorscheme;
			
			if (e.opts.menu_colormap.isSelected())
			{
				int i1 = (int) (e.nx*0.85);
				int i2 = (int) (e.nx*0.92);

				int j1 = (int) (e.ny*0.82);
				int j2 = (int) (e.ny*0.92);
				

				for (int i = i1; i <= i2; i++) {
					for (int j = j1; j <= j2; j++) {
						double t = 0;
						if (colorscheme == ScalarView.ColorScheme.RED_BLUE || colorscheme == ScalarView.ColorScheme.CYAN_YELLOW || scalarview == ScalarView.CHARGE) {
							t = 2*((j2-j)/(double)(j2-j1)-0.5);
							scalarfield[i][j] = t/scalingconstant + scalar_offset;
						} else if (colorscheme == ScalarView.ColorScheme.WHITE || colorscheme == ScalarView.ColorScheme.GREEN) {
							t = (j2-j)/(double)(j2-j1);
							scalarfield[i][j] = t/scalingconstant + scalar_offset;
						} else if (scalarview == ScalarView.LIGHT) {
							t = (j2-j)/(double)(j2-j1);
							scalarfield[i][j] = -(400+(650-400)*t);
						}
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

			if (colorscheme == ScalarView.ColorScheme.RED_BLUE) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float v = (float) (scalarfield[i][j]-scalar_offset)*scalingconstant;
						setColorFloat(v, 0, -v);
						setPixel(i, j);
					}
				}
			} else if (colorscheme == ScalarView.ColorScheme.CYAN_YELLOW) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float v = (float) (scalarfield[i][j]-scalar_offset)*scalingconstant;
						setColorFloat(v, Math.abs(v), -v);
						setPixel(i, j);
					}
				}
			} else if (colorscheme == ScalarView.ColorScheme.GREEN) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float v = (float) (scalarfield[i][j]-scalar_offset)*scalingconstant;
						setColorFloat(0, Math.abs(v), 0);
						setPixel(i, j);
					}
				}
			} else if (colorscheme == ScalarView.ColorScheme.WHITE) {
				for (int i = 1; i < e.nx-1; i++) {
					for (int j = 1; j < e.ny-1; j++) {
						float v = (float) (scalarfield[i][j]-scalar_offset)*scalingconstant;
						setColorFloat(v, v, v);
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
				} else if (scalarview == ScalarView.LIGHT) {
					for (int i = 1; i < e.nx-1; i++) {
						for (int j = 1; j < e.ny-1; j++) {
							if (scalarfield[i][j] < 0) {
								setColorWavelength(-scalarfield[i][j], 1);
							} else {
								double hc = 1.986e-25;
								double Emax = (e.materials[i][j].Ec - e.materials[i][j].Ev) + 0.5/e.beta; // Peak emission energy
								double lambda_nm = 1e9*hc/Emax;
								setColorWavelength(lambda_nm, scalarfield[i][j]*scalingconstant);
							}
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

		if (brush == Brush.ZOOM && e.controls.mouse_pressed_prev && !e.controls.shift_down) {
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
					p0.draw(this);
				}
				
				if (e.controls.moving_selection) {
					for (Probe p0 : e.controls.selection.probes) {
						p0.translate(e.controls.delta_mx, e.controls.delta_my);
						p0.draw(this);
						p0.translate(-e.controls.delta_mx, -e.controls.delta_my);
					}
				}
			}

			for (Plot p: e.plots) {
				if (p.frame.isVisible()) {
					if (p.path != null)
						p.path.draw(this);
				}
			}
			
			if (e.controls.plotpath != null)
				e.controls.plotpath.draw(this);
		}

		setalphaFG(0.8);
		setColorFloat(1.0f, 1.0f, 1.0f);
		
		if (e.controls.texting) {
			drawPixelLine(e.controls.text_x, e.controls.text_y, e.controls.text_x, e.controls.text_y+7);
		}

		stampPixelData();

	}
	
	void drawOverlay() {
		delta_t = e.time - t_prev;
		t_prev = e.time;

		show_carriers = e.opts.gui_carriers.isSelected();
		vector_display_mode = e.controls.vectormode.getOption();
		scalar_display_mode = e.controls.scalarmode.getOption();
		scalar_view = e.controls.scalarview.getOption();
		vector_view = e.controls.vectorview.getOption();
		
		try {
			if (SemiSim.instance.graphics_threads.size() == SemiSim.n_threads) {
				graphics_start_barrier.await();
				graphics_end_barrier.await();
			}
		} catch (InterruptedException | BrokenBarrierException e) {
			e.printStackTrace();
		}
	}
	
	void drawText(Graphics2D g) {
		clearStrings();

		//Graphics2D g = (Graphics2D) img_back.getGraphics();
		g.setRenderingHint(
		RenderingHints.KEY_TEXT_ANTIALIASING,
		RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		if (e.opts.menu_probes.isSelected()) {
			for (Probe p: e.probes) {
				drawMonospacedStringSimCoords(p.getText(e.units), p.labelcoord.x, p.labelcoord.y, g);
			}

			drawStrings(g);
			startNewStringLayer();
		}
		
		ScalarView scalarview = e.controls.scalarview.getOption();
		ScalarMode scalarmode = e.controls.scalarmode.getOption();
		if (scalarview != ScalarView.NONE && scalarmode != ScalarMode.NONE) {
			ScalarView.ColorScheme colorscheme = scalarview.colorscheme;

			if (e.opts.menu_colormap.isSelected())
			{
				int i1 = (int) (e.nx*0.85);
				int i2 = (int) (e.nx*0.92);

				int j1 = (int) (e.ny*0.82);
				int j2 = (int) (e.ny*0.95);

				double tmin = 0;
				double tmax = 1;
				
				if (!(colorscheme == ScalarView.ColorScheme.GREEN || colorscheme == ScalarView.ColorScheme.WHITE))
				{
					tmin = 2*(tmin-0.5);
					tmax = 2*(tmax-0.5);
				}
				if (scalarview == ScalarView.LIGHT) {
					drawMonospacedStringSimCoords(e.units.toString(400e-9, Quantity.LENGTH), (i1+i2)/2, j2, g);
					texts.get(texts.size()-1).isHorizontalCentered = true;
					drawMonospacedStringSimCoords(e.units.toString(650e-9, Quantity.LENGTH), (i1+i2)/2, j1, g);
					texts.get(texts.size()-1).isHorizontalCentered = true;
				} else {
					drawMonospacedStringSimCoords(e.units.toString(tmin/scalingconstant+scalar_offset, scalarview.unit), (i1+i2)/2, j2, g);
					texts.get(texts.size()-1).isHorizontalCentered = true;
					drawMonospacedStringSimCoords(e.units.toString(tmax/scalingconstant+scalar_offset, scalarview.unit), (i1+i2)/2, j1, g);
					texts.get(texts.size()-1).isHorizontalCentered = true;
				}
				texts.get(texts.size()-1).isBottomJustified = true;
			}
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
			int voffset = 1 + (int)(e.controls.my_screen);
			int hoffset = 5 + (int)(e.controls.mx_screen)+15;
			//int voffset = 1 + (int)(e.controls.my*scalefactor);
			//int hoffset = 5 + (int)(e.controls.mx*scalefactor)+15;

			if (e.opts.menu_tooltip.isSelected()) {
				if (voffset + 22*vspacing > e.canvas.getHeight()) {
					voffset = voffset - ((voffset + 22*vspacing) - e.canvas.getHeight());
				}
				if (hoffset + 120 > e.canvas.getWidth()) {
					hoffset = hoffset - ((hoffset + 120) - e.canvas.getWidth());
				}
			}

			if (e.opts.menu_materialname.isSelected()) {
				String name = "Material: " + mat.getDisplayName();
				this.drawBigString(name, hoffset, voffset + 1*vspacing, g);
			}

			if (e.opts.menu_tooltip.isSelected()) {
				voffset = voffset+3;
				int line = 2;
				drawTwoColumnString("E" , 							e.units.toString(Utils.bilinearinterp_length(e.Ex, e.Ey, mx, my, e.nx, e.ny), Quantity.ELECTRIC_FIELD, 1e-6), 				hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("B" , 							e.units.toString(e.parity*Utils.bilinearinterp(e.Bz, mx-0.5, my-0.5, e.nx, e.ny), Quantity.MAGNETIC_FLUX_DENSITY, 1e-9),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03d5" , 						e.units.toString(e.phi[mx][my], Quantity.ELECTRIC_POTENTIAL, 1e-6),				hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u2130" , 						e.units.toString(mat.emf, Quantity.ELECTRIC_FIELD, 1e-6),					hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03b5/\u03b5\u2080" , 		e.units.toString(mat.eps_r, Quantity.DIMENSIONLESS),							hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03bc/\u03bc\u2080" , 		e.units.toString(mat.mu_r, Quantity.DIMENSIONLESS),							hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u2099" , 				e.units.toString(e.rho_n[mx][my], Quantity.CHARGE_DENSITY, 1e-9),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u209A" , 				e.units.toString(e.rho_p[mx][my], Quantity.CHARGE_DENSITY, 1e-9),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1\u2080",					e.units.toString(e.rho_back[mx][my], Quantity.CHARGE_DENSITY, 1e-9),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("\u03c1" , 						e.units.toString(e.rho_free[mx][my], Quantity.CHARGE_DENSITY, 1e-9),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J\u2099" ,						e.units.toString(Utils.bilinearinterp_length(e.Jx_n, e.Jy_n, mx, my, e.nx, e.ny), Quantity.CURRENT_DENSITY, 1),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J\u209A" , 					e.units.toString(Utils.bilinearinterp_length(e.Jx_p, e.Jy_p, mx, my, e.nx, e.ny), Quantity.CURRENT_DENSITY, 1),			hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("J" , 							e.units.toString(Utils.bilinearinterp_length(e.Jx_free, e.Jy_free, mx, my, e.nx, e.ny), Quantity.CURRENT_DENSITY, 1),		hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("F\u2099" , 					e.units.toString(-e.mu_n[mx][my]/e.q_n, Quantity.ELECTRIC_POTENTIAL, 1e-9),						hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("F\u209a" , 					e.units.toString(-e.mu_p[mx][my]/e.q_p, Quantity.ELECTRIC_POTENTIAL, 1e-9),						hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("V" , 							e.units.toString(e.V_avg[mx][my], Quantity.ELECTRIC_POTENTIAL, 1e-9),				hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("x" , 							e.units.toString(mx*e.ds, Quantity.LENGTH),							hoffset, voffset + line*vspacing, g); line++;
				drawTwoColumnString("y" , 							e.units.toString(e.ds*e.ny-(my+1)*e.ds, Quantity.LENGTH),			hoffset, voffset + line*vspacing, g); line++;
			}

			drawStrings(g);
			startNewStringLayer();
		}

		int vspacing = 13;
		int voffset = 3;
		int hoffset = 5;
		int line = 1;
		
		if (e.opts.menu_time.isSelected()) {
			drawString("Time: " + e.units.toString(e.time, Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString("Steps/s: " + e.units.toString(e.opts.gui_simspeed_2.getValue()/e.simFPStimer.getAverageTime(), Quantity.DIMENSIONLESS), hoffset, voffset + line*vspacing, g); line++;

			String sv_a = e.controls.scalarmode.getOption() != ScalarMode.NONE? (e.controls.scalarmode.getOption().shorthand + ": " + e.controls.scalarview.getOption().shorthand) : "";
			String sv_b = e.controls.vectormode.getOption() != VectorMode.NONE? (e.controls.vectormode.getOption().shorthand + ": " + e.controls.vectorview.getOption().shorthand) : "";
			String sv = Arrays.asList(sv_a, sv_b).stream().filter(s -> s != null && !s.isEmpty()).collect(Collectors.joining(" / "));
			drawString(sv, hoffset, voffset + line*vspacing, g); line++;
			if (e.opts.gui_paused.isSelected())
			{
				drawString("Paused", hoffset, voffset + line*vspacing, g); line++;
			}
			if (e.sign_violation_timer > 0) {
				e.sign_violation_timer--;
				drawString("Warning: Numerical instability detected. Please decrease timestep or increase junction smoothing.", hoffset, voffset + line*vspacing, g); line++;
			}
		}
		
		if (e.controls.debugging) {
			long total = Runtime.getRuntime().totalMemory();
			long used  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			drawString("Used memory " + e.units.toString(used, Quantity.INFORMATION), hoffset, voffset + line*vspacing, g); line++;
			drawString("Total memory " + e.units.toString(total, Quantity.INFORMATION), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t4.getName() + " " + e.units.toString(e.t4.getAverageTime(), Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString(t5.getName() + " " + e.units.toString(t5.getAverageTime(), Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t6.getName() + " " + e.units.toString(e.t6.getAverageTime()*e.opts.gui_simspeed_2.getValue(), Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t7.getName() + " " + e.units.toString(e.t7.getAverageTime(), Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.t8.getName() + " " + e.units.toString(e.t8.getAverageTime(), Quantity.TIME), hoffset, voffset + line*vspacing, g); line++;
			drawString(FPStimer.getName() + " " + e.units.toString(1/FPStimer.getAverageTime(), Quantity.FREQUENCY), hoffset, voffset + line*vspacing, g); line++;
			drawString(e.simFPStimer.getName() + " " + e.units.toString(1/e.simFPStimer.getAverageTime(), Quantity.FREQUENCY), hoffset, voffset + line*vspacing, g); line++;
		}
		if (carrier_diffusion_warning_timer > 0) {
			drawError("Error: Metal cannot touch simulation boundary when carrier diffusion view is enabled.", hoffset, voffset + line*vspacing, g); line ++;
		}
		if (e.numerical_overflow) {
			drawError("Error: Numerical overflow detected. Please reset simulation.", hoffset, voffset + line*vspacing, g); line ++;
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
			System.out.println("Graphics thread " + n_thread + " initialized.");
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

						
						boolean conductors_only = vector_view.isConductorOnly();
						
						double[][][] vf = {null, null};
						e.computeVectorField(vf, vector_view);
						double[][] vf_x = vf[0];
						double[][] vf_y = vf[1];

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
									if (conductors_only && Utils.bilinearinterp(e.conducting, d.x, d.y, e.nx, e.ny) == 0) {
										d.lifespan = 0;
										d.time = 0;
									}
								}
							}

							setColorFloat(1, 1, 1);
							int dot_offset = (scalefactor-1)/2;
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
									drawRectangle((int)((d.x+0.5)*scalefactor)-dot_offset, (int)((d.y+0.5)*scalefactor)-dot_offset, scalefactor, scalefactor,
									1f, 1f, 1f, (float)alphaFG, 1f);
								}
							}
						}
					}
					

					if (scalar_view != ScalarView.NONE && (scalar_display_mode == ScalarMode.CONTOUR_COLORS || scalar_display_mode == ScalarMode.CONTOUR)) {
						graphics_mid_barrier.await();
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
					
					
					if (show_carriers) {

						if ((!e.opts.gui_paused.isSelected() || delta_t > 0)) {
							if (n_thread == 0 && carrier_diffusion_warning_timer > 0) {
								carrier_diffusion_warning_timer--;
							}
							
							double C = cc_default_dot_density*Math.pow(10.0, e.opts.gui_carrier_density.getValue()/20.0);

							if (n_thread == 0) {
								rho_n_dist.prepare(e.rho_n);
								rho_p_dist.prepare(e.rho_p);
								rho_G_dist.prepare(e.G);
							}
							
							int i_low = lower(ccdots.size());
							int i_high = upper(ccdots.size());

							graphics_mid_barrier.await();

							double A = e.ds*e.ds;
							double N_G = rho_G_dist.getTotalAmount()*A*e.e_charge*C*delta_t;
							double N_n = rho_n_dist.getTotalAmount()*A*C*delta_t/tau;
							double N_p = rho_p_dist.getTotalAmount()*A*C*delta_t/tau;

							double N_n_excess = Math.max(C-C_prev, 0)*rho_n_dist.getTotalAmount()*A;
							double N_p_excess = Math.max(C-C_prev, 0)*rho_p_dist.getTotalAmount()*A;

							double P_deficit = Math.max(-(C-C_prev)/C_prev, 0);

							boolean show_gen_recomb = e.opts.menu_gen_recomb.isSelected();
							
							for (int i = i_low; i < i_high; i++) {
								ChargeCarrierDot d = ccdots.get(i);
								if (d == null) continue;
								
								d.time -= delta_t;
								if (d.time < 0) {
									ccdots.remove(i);
									i--;
								}
								else {
									double R_tmp = Utils.bilinearinterp(e.R, d.x, d.y, e.nx, e.ny)*e.e_charge;
									if (d.type == DotType.HOLE) {
										if (R_tmp > 0 && frand.next() < R_tmp/Utils.bilinearinterp(e.rho_p, d.x, d.y, e.nx, e.ny)*delta_t) {
											if (show_gen_recomb) {
												ccdots.replace(i, new ChargeCarrierDot(d.x, d.y, tau_events, tau_events*0.5, DotType.RECOMBINATION, 1));
											} else {
												ccdots.remove(i);
												i--;
											}
											continue;
										}
									} else if (d.type == DotType.ELECTRON) {
										if (R_tmp > 0 && frand.next() < -R_tmp/Utils.bilinearinterp(e.rho_n, d.x, d.y, e.nx, e.ny)*delta_t) {
											if (show_gen_recomb) {
												ccdots.replace(i, new ChargeCarrierDot(d.x, d.y, tau_events, tau_events*0.5, DotType.RECOMBINATION, 1));
											} else {
												ccdots.remove(i);
												i--;
											}
											continue;
										}
									} else {
										if (Utils.bilinearinterp(e.semiconducting, d.x, d.y, e.nx, e.ny) == 0) {
											d.time -= 5*delta_t;
										}
									}

									if (P_deficit > 0 && frand.next() < P_deficit) {
										ccdots.remove(i);
										i--;
										continue;
									}
								}
							}

							graphics_mid_barrier.await();

							rho_n_dist.generateSamples(N_n/n_threads, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau, DotType.ELECTRON, 0)); });
							rho_p_dist.generateSamples(N_p/n_threads, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau, DotType.HOLE, 0)); });
							rho_G_dist.generateSamples(N_G/n_threads, (c) -> {
								ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), DotType.ELECTRON, 0));
								if (show_gen_recomb) ccdots.add(new ChargeCarrierDot(c.x, c.y, tau_events, tau_events*0.5, DotType.GENERATION, 1));
							});
							rho_G_dist.generateSamples(N_G/n_threads, (c) -> {
								ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), DotType.HOLE, 0));
								if (show_gen_recomb) ccdots.add(new ChargeCarrierDot(c.x, c.y, tau_events, tau_events*0.5, DotType.GENERATION, 1));
							});
							rho_n_dist.generateSamples(N_n_excess/n_threads, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), DotType.ELECTRON, 0)); });
							rho_p_dist.generateSamples(N_p_excess/n_threads, (c) -> { ccdots.add(new ChargeCarrierDot(c.x, c.y, tau, tau*frand.next(), DotType.HOLE, 0)); });

							C_prev = C;
						}

						graphics_mid_barrier.await();

						boolean fast = e.opts.menu_hide_carriers_metal.isSelected();
						boolean show_diffusion = e.opts.menu_carrier_diffusion.isSelected();
						
						int lower = lower(ccdots.size());
						int upper = upper(ccdots.size());

						int steps = fast? 2 : 10;
						double dt_dot = delta_t/steps;

						try {
							double factor = Math.sqrt(24*dt_dot);
							int dot_offset = (scalefactor-1)/2;
							
							for (int i = lower; i < upper; i++) {
								ChargeCarrierDot d = ccdots.get(i);

								if (d.time > 0) {

									boolean dorender = !fast || Utils.bilinearinterp(e.semiconducting, d.x, d.y, e.nx, e.ny) > 0;

									if (delta_t > 0) {
										if (!show_diffusion) {
											if (d.type == DotType.ELECTRON) {
												for (int k = 0; k < steps; k++) {
													double s = dt_dot/(e.ds*Utils.bilinearinterp(e.rho_n, d.x, d.y, e.nx, e.ny));
													d.x += s*Utils.bilinearinterp(e.Jx_n, d.x-0.5, d.y, e.nx, e.ny);
													d.y += s*Utils.bilinearinterp(e.Jy_n, d.x, d.y-0.5, e.nx, e.ny);
												}
											} else if (d.type == DotType.HOLE) {
												for (int k = 0; k < steps; k++) {
													double s = dt_dot/(e.ds*Utils.bilinearinterp(e.rho_p, d.x, d.y, e.nx, e.ny));
													d.x += s*Utils.bilinearinterp(e.Jx_p, d.x-0.5, d.y, e.nx, e.ny);
													d.y += s*Utils.bilinearinterp(e.Jy_p, d.x, d.y-0.5, e.nx, e.ny);
												}
											}
										} else {
											if (d.type == DotType.ELECTRON) {
												double sqrt_D_n = Utils.bilinearinterp(e.sqrt_D_eff_n, d.x, d.y, e.nx, e.ny);
												double s = dt_dot/e.ds;
												double t = factor*sqrt_D_n/e.ds; //Random walk PDF obeys diffusion equation. Variance of uniform dist is 12L^2 and variance of heat kernel is 2Dt
												for (int k = 0; k < steps; k++) {
													double dx_diff = t*(frand.next()-0.5) + s*Utils.bilinearinterp(e.vel_x_n, d.x-0.5, d.y, e.nx, e.ny);
													double dy_diff = t*(frand.next()-0.5) + s*Utils.bilinearinterp(e.vel_y_n, d.x, d.y-0.5, e.nx, e.ny);
													if (e.conducting[(int)(d.x+dx_diff+0.5)][(int)(d.y+dy_diff+0.5)] == 1) {
														d.x += dx_diff;
														d.y += dy_diff;
													}
												}
											} else if (d.type == DotType.HOLE) {
												double sqrt_D_p = Utils.bilinearinterp(e.sqrt_D_eff_p, d.x, d.y, e.nx, e.ny);
												double s = dt_dot/e.ds;
												double t = factor*sqrt_D_p/e.ds;
												
												for (int k = 0; k < steps; k++) {
													double dx_diff = t*(frand.next()-0.5) + s*Utils.bilinearinterp(e.vel_x_p, d.x-0.5, d.y, e.nx, e.ny);
													double dy_diff = t*(frand.next()-0.5) + s*Utils.bilinearinterp(e.vel_y_p, d.x, d.y-0.5, e.nx, e.ny);
													if (e.conducting[(int)(d.x+dx_diff+0.5)][(int)(d.y+dy_diff+0.5)] == 1) {
														d.x += dx_diff;
														d.y += dy_diff;
													}
												}
											}
										}

										double p = d.time/d.lifespan;
										d.brightness = 1.0*bump(p, 1/3.0);
									}

									double alphaFG = d.brightness;
									if (dorender) {
										if (d.type == DotType.ELECTRON)
											drawRectangle((int)((d.x+0.5)*scalefactor)-dot_offset, (int)((d.y+0.5)*scalefactor)-dot_offset, scalefactor, scalefactor,
											0.25f, 0.25f, 1f, (float)alphaFG, 1-(float)alphaFG, d.random_id);
										else if (d.type == DotType.HOLE)
											drawRectangle((int)((d.x+0.5)*scalefactor)-dot_offset, (int)((d.y+0.5)*scalefactor)-dot_offset, scalefactor, scalefactor,
											1f, 0.25f, 0.25f, (float)alphaFG, 1-(float)alphaFG, d.random_id);
										else if (d.type == DotType.GENERATION) {
											drawRectangle((int)((d.x+0.5)*scalefactor)-dot_offset, (int)((d.y+0.5)*scalefactor)-dot_offset, scalefactor, scalefactor,
											0f, 0f, 0f, (float)alphaFG, 1-(float)alphaFG, d.random_id);
										} else if (d.type == DotType.RECOMBINATION) {
											drawRectangle((int)((d.x+0.5)*scalefactor)-dot_offset, (int)((d.y+0.5)*scalefactor)-dot_offset, scalefactor, scalefactor,
											1f, 1f, 1f, (float)alphaFG, 1-(float)alphaFG, d.random_id);
										}
									}

								}
							}
						} catch (ArrayIndexOutOfBoundsException e) {
							carrier_diffusion_warning_timer = 20;
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
			if (text.isRightJustified) text.x -= (text.width-4);
			if (text.isBottomJustified) text.y -= (text.height-4);
			if (text.isHorizontalCentered) text.x -= (text.width/2-2);
			if (text.isVerticalCentered) text.y -= (text.height/2-2);
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

	public void drawMonospacedStringSimCoords(String str1, int x, int y, Graphics g) {

		double sf_x = (double)e.canvas.zoom_bound_x/(e.controls.zoom_i2-e.controls.zoom_i1+1);
		double sf_y = (double)e.canvas.zoom_bound_y/(e.controls.zoom_j2-e.controls.zoom_j1+1);

		drawMonospacedString(str1, (int) ((x-e.controls.zoom_i1+0.5)*sf_x+e.canvas.offset_x), (int) ((y-e.controls.zoom_j1+0.5)*sf_y+e.canvas.offset_y), g);
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
		boolean isRightJustified = false;
		boolean isBottomJustified = false;
		boolean isHorizontalCentered = false;
		boolean isVerticalCentered = false;

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
		NONE("No scalar overlay",															"None",		Quantity.DIMENSIONLESS,				ColorScheme.OTHER,			1),
		E_FIELD("View: E field magnitude",													"E",		Quantity.ELECTRIC_FIELD,			ColorScheme.GREEN,			1e5),
		B_FIELD("View: B field",															"B",		Quantity.MAGNETIC_FLUX_DENSITY,		ColorScheme.CYAN_YELLOW,	1e-5),
		H_FIELD("View H field",																"H",		Quantity.MAGNETIC_FIELD_STRENGTH,	ColorScheme.CYAN_YELLOW,	1e-5/1.257e-6),
		POTENTIAL("View \u03d5: Electric scalar potential", 								"\u03d5",	Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE,		1),
		ENERGY("View u: Electromagnetic energy density",									"u",		Quantity.ENERGY_DENSITY,			ColorScheme.GREEN,			1),
		CURRENT("View J: Total current magnitude",											"J",		Quantity.CURRENT_DENSITY,			ColorScheme.GREEN,			1e7),
		CHARGE("View \u03c1: Net charge density",											"\u03c1",	Quantity.CHARGE_DENSITY,			ColorScheme.OTHER,			1e1),
		ELECTRON_CHARGE("View \u03c1\u2099: Electron charge density",						"\u03c1\u2099",	Quantity.CHARGE_DENSITY,		ColorScheme.RED_BLUE,		1e1),
		HOLE_CHARGE("View \u03c1\u209A: Hole charge density",								"\u03c1\u209A",	Quantity.CHARGE_DENSITY,		ColorScheme.RED_BLUE,		1e1),
		BACKGROUND_CHARGE("View \u03c1\u2080: Static charge density (doping)",				"\u03c1\u2080",	Quantity.CHARGE_DENSITY,		ColorScheme.RED_BLUE,		1e1),
		COMBINED_CHARGE("View: Combined electron+hole charge density",						"\u03c1\u2099 and \u03c1\u209A",	Quantity.DIMENSIONLESS,	ColorScheme.OTHER,			1),
		HEAT("View Q: Heat dissipation",													"q",		Quantity.POWER_DENSITY,				ColorScheme.RED_BLUE,		1e12),
		ENTROPY("View s: Entropy generation rate",											"s",		Quantity.ENTROPY_DENSITY_RATE,		ColorScheme.RED_BLUE,		3.33e9),
		ELECTRON_POTENTIAL("View F\u2099: Electron quasi Fermi level",						"μ\u2099",	Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE, 		1),
		HOLE_POTENTIAL("View F\u209A: Hole quasi Fermi level",								"μ\u209A",	Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE, 		1),
		AVERAGE_POTENTIAL("View V: Voltage",												"V",		Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE,		1),
		GENERATION("View G: Carrier generation rate",										"G",		Quantity.RATE_DENSITY,				ColorScheme.GREEN,			1e31),
		RECOMBINATION("View R: Carrier recombination rate",									"R",		Quantity.RATE_DENSITY,				ColorScheme.GREEN,			1e31),
		LIGHT("View: Emitted light",														"Light",	Quantity.DIMENSIONLESS,				ColorScheme.OTHER,			1e30),
		ELECTRON_DENSITY("View: Electron density",											"ne",		Quantity.NUMBER_DENSITY,			ColorScheme.GREEN,			1e23),
		HOLE_DENSITY("View: Hole density",													"nh",		Quantity.NUMBER_DENSITY,			ColorScheme.GREEN,			1e23),
		ELECTRON_VOLTAGE("View: Electron voltage",											"Ve",		Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE,		1),
		HOLE_VOLTAGE("View: Hole voltage",													"Vh",		Quantity.ELECTRIC_POTENTIAL,		ColorScheme.RED_BLUE,		1),
		ELECTRON_VEL("View: Electron drift velocity",										"ve",		Quantity.VELOCITY,					ColorScheme.GREEN,			1e6),
		HOLE_VEL("View: Hole drift velocity",												"vh",		Quantity.VELOCITY,					ColorScheme.GREEN,			1e6),
		RECOMB_RAD("View: Radiative recombination rate",									"R (rad)",	Quantity.RATE_DENSITY,				ColorScheme.GREEN,			1e31),
		RECOMB_SRH("View: SRH recombination rate",											"R (SRG)",	Quantity.RATE_DENSITY,				ColorScheme.GREEN,			1e31),
		RECOMB_AUGER("View: Auger recombination rate",										"R (aug)",	Quantity.RATE_DENSITY,				ColorScheme.GREEN,			1e31),
		DEBUG("Debug",																		"Debug",	Quantity.DIMENSIONLESS,				ColorScheme.RED_BLUE,		1e-3);
	
		enum ColorScheme {
			RED_BLUE, CYAN_YELLOW, GREEN, WHITE, OTHER;
		}
	
		public String name;
		public Quantity unit;
		public String shorthand;
		public ColorScheme colorscheme;
		public double scale; //Typical order of magnitude of the quantity
	
		ScalarView(String name, String shorthand, Quantity unit, ColorScheme colorScheme, double scale)
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
		NONE("No vector overlay",								"None",				Quantity.DIMENSIONLESS,			1),
		E_FIELD("View E: Electric field",						"E",				Quantity.ELECTRIC_FIELD,		10),
		D_FIELD("View D: Displacement field",					"D",				Quantity.ELECTRIC_FLUX_DENSITY,	10*8.85e-12),
		ELECTRON_CURRENT("View J\u2099: Electron current",		"J\u2099",			Quantity.CURRENT_DENSITY,		1),
		HOLE_CURRENT("View J\u209A: Hole current",				"J\u209A",			Quantity.CURRENT_DENSITY,		1),
		TOTAL_CURRENT("View J: Total current",					"J",				Quantity.CURRENT_DENSITY,		1),
		EMF("View \u2130: External electromotive force",		"\u2130",			Quantity.ELECTRIC_FIELD,		1),
		POYNTING("View S: Poynting vector",						"S",				Quantity.INTENSITY,				1),
		ELECTRON_DRIFT("View: Electron drift current",			"Jn (drift)",		Quantity.CURRENT_DENSITY,		1),
		ELECTRON_DIFFUSION("View: Electron diffusion current",	"Jn (diffusion)",	Quantity.CURRENT_DENSITY,		1),
		ELECTRON_VELOCITY("View: Electron drift velocity",		"vn",				Quantity.VELOCITY,				1),
		HOLE_DRIFT("View: Hole drift current",					"Jp (drift)",		Quantity.CURRENT_DENSITY,		1),
		HOLE_DIFFUSION("View: Hole diffusion current",			"Jp (diffusion)",	Quantity.CURRENT_DENSITY,		1),
		HOLE_VELOCITY("View: Hole drift velocity",				"vp",				Quantity.VELOCITY,				1);
	
		public String name;
		public Quantity unit;
		public String shorthand;
		public double scale;
	
		VectorView(String name, String shorthand, Quantity unit, double scale)
		{
			this.name = name;
			this.shorthand = shorthand;
			this.unit = unit;
			this.scale = scale;
		}
	
		@Override
		public String toString() {
			return name;
		}
		
		public boolean isConductorOnly() {
			return (!(this == E_FIELD || this == D_FIELD || this == POYNTING));
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
		DotType type;
		int random_id;
		
		public ChargeCarrierDot(double x, double y, double lifespan, double time, DotType type, int random_id) {
			super(x, y, lifespan, time);
			this.type = type;
			this.random_id = random_id;
		}
	}
	
	public enum DotType {
		ELECTRON, HOLE, GENERATION, RECOMBINATION;
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

			int canvas_x = this.getWidth();
			int canvas_y = this.getHeight();
			
			int xw = e.controls.zoom_i2 - e.controls.zoom_i1 + 1;
			int yw = e.controls.zoom_j2 - e.controls.zoom_j1 + 1;
			
			if (xw/(double)yw >= canvas_x/(double) canvas_y) {
				int dim2 = (int) (canvas_x*yw/(double)xw);
				zoom_bound_x = canvas_x - 1;
				zoom_bound_y = dim2 - 1;
				offset_x = 0;
				offset_y = (canvas_y - dim2)/2;
			} else {
				int dim2 = (int) (canvas_y*xw/(double)yw);
				zoom_bound_x = dim2 - 1;
				zoom_bound_y = canvas_y - 1;
				offset_x = (canvas_x - dim2)/2;
				offset_y = 0;
			}
			
			g.drawImage(e.renderer.img_front, offset_x, offset_y, zoom_bound_x+1+offset_x, zoom_bound_y+1+offset_y, 
				e.controls.zoom_i1*e.renderer.scalefactor, e.controls.zoom_j1*e.renderer.scalefactor, (e.controls.zoom_i2+1)*e.renderer.scalefactor, (e.controls.zoom_j2+1)*e.renderer.scalefactor, e.opts);

			e.renderer.drawText((Graphics2D)g);
		}

		public RenderCanvas(Simulation w) {
			e = w;
			setPreferredSize(new Dimension(768, 768));
		}
	}
}
