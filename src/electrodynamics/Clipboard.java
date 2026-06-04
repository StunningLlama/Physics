// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.util.ArrayList;
import java.util.List;

import electrodynamics.probe.Probe;
import electrodynamics.util.Utils;

public class Clipboard {
	
	ClipboardMaterial[][] mat;
	List<Probe> probes;
	Simulation e;
	
	public Clipboard(Simulation e) {
		this.e = e;
		mat = new ClipboardMaterial[e.nx][e.ny];
		probes = new ArrayList<Probe>();
	}

	public void clear() {
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (mat[i][j] == null)
					mat[i][j] = new ClipboardMaterial();
				else
					mat[i][j].erase();
			}
		}
		
		probes.clear();
	}
	
	public void cut() {
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (mat[i][j] == null)
					mat[i][j] = new ClipboardMaterial();
				mat[i][j].erase();
				if (e.controls.selected[i][j]) {
					mat[i][j] = new ClipboardMaterial(e, i, j);
					e.controls.selected[i][j] = false;
					e.eraseMaterial(i, j);
				}
			}
		}

		probes = Utils.cloneList(e.probes, Probe::clone, (p) -> p.selected);
		e.probes.removeIf((p) -> p.selected);
	}

	public void paste(int delta_mx, int delta_my) {
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				int si = i-delta_mx;
				int sj = j-delta_my;
				if (si >= 0 && sj >= 0 && si < e.nx && sj < e.ny && mat[si][sj].m.type != MaterialType.VACUUM) {
					e.eraseMaterial(i, j);
					mat[si][sj].paste(e, i, j);
					e.controls.selected[i][j] = true;
				}
			}
		}


		for (Probe p : probes) {
			p.translate(delta_mx, delta_my);
			p.selected = true;
			if (p.checkInBounds(e)) {
				e.probes.add(p.clone());
			}
		}
	}
	
	public void transfer(Clipboard target) {
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				e.controls.selected[i][j] = false;
			}
		}
		
		for (Probe p : probes)
			p.selected = false;

		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				target.mat[i][j] = mat[i][j].clone();
			}
		}
		
		target.probes = Utils.cloneList(probes, Probe::clone);
	}
	
	public boolean copy(boolean deleteselection) {
		int i_min = e.nx-1;
		int j_min = e.ny-1;
		boolean clipboardempty = true;
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				mat[i][j].erase();
				if (e.controls.selected[i][j] && e.materials[i][j].type != MaterialType.VACUUM) {
					if (i < i_min) i_min = i;
					if (j < j_min) j_min = j;
				}
			}
		}
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (e.controls.selected[i][j] && e.materials[i][j].type != MaterialType.VACUUM) {
					
					mat[i-i_min][j-j_min] = new ClipboardMaterial(e, i, j);
					
					if (deleteselection) {
						e.eraseMaterial(i, j);
					}
					
					clipboardempty = false;
				}
				if (deleteselection) {
					e.controls.selected[i][j] = false;
				}
			}
		}
		
		probes = Utils.cloneList(e.probes, Probe::clone, (p) -> p.selected);
		
		if (clipboardempty && !probes.isEmpty()) {
			i_min = e.controls.mx_start;
			j_min = e.controls.my_start;
			
			clipboardempty = false;
		}
		
		for (Probe p : probes)
			p.translate(-i_min, -j_min);
		

		if (deleteselection) {
			e.probes.removeIf((p) -> p.selected);
		}
		
		return clipboardempty;
	}
	
	public void flip_h() {
		int i_max = 0;
		int j_max = 0;
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (mat[i][j].m.type != MaterialType.VACUUM) {
					if (i > i_max) i_max = i;
					if (j > j_max) j_max = j;
				}
			}
		}
		
		ClipboardMaterial[][] new_selection = new ClipboardMaterial[i_max+1][j_max+1];
		for (int i = 0; i <= i_max; i++)
		{
			for (int j = 0; j <= j_max; j++)
			{
				new_selection[i_max-i][j] = mat[i][j].clone();
				new_selection[i_max-i][j].flip_h();
				mat[i][j].erase();
			}
		}
		
		for (int i = 0; i <= i_max; i++)
		{
			for (int j = 0; j <= j_max; j++)
			{
				mat[i][j] = new_selection[i][j];
			}
		}
		
		for (Probe p : probes) p.flip_h(0, i_max);
	}
	
	public void flip_v() {
		int i_max = 0;
		int j_max = 0;
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (mat[i][j].m.type != MaterialType.VACUUM) {
					if (i > i_max) i_max = i;
					if (j > j_max) j_max = j;
				}
			}
		}
		
		ClipboardMaterial[][] new_selection = new ClipboardMaterial[i_max+1][j_max+1];
		for (int i = 0; i <= i_max; i++)
		{
			for (int j = 0; j <= j_max; j++)
			{
				new_selection[i][j_max - j] = mat[i][j].clone();
				new_selection[i][j_max - j].flip_v();
				mat[i][j].erase();
			}
		}
		
		for (int i = 0; i <= i_max; i++)
		{
			for (int j = 0; j <= j_max; j++)
			{
				mat[i][j] = new_selection[i][j];
			}
		}

		for (Probe p : probes) p.flip_v(0, j_max);
	}

	public void rotate90() {
		int i_max = 0;
		int j_max = 0;
		for (int i = 0; i < e.nx; i++)
		{
			for (int j = 0; j < e.ny; j++)
			{
				if (mat[i][j].m.type != MaterialType.VACUUM) {
					if (i > i_max) i_max = i;
					if (j > j_max) j_max = j;
				}
			}
		}

		ClipboardMaterial[][] new_selection = new ClipboardMaterial[j_max+1][i_max+1];
		for (int i = 0; i <= i_max; i++)
		{
			for (int j = 0; j <= j_max; j++)
			{
				new_selection[j_max-j][i] = mat[i][j].clone();
				new_selection[j_max-j][i].rotate90();
				mat[i][j].erase();
			}
		}

		for (int i = 0; i <= j_max; i++)
		{
			for (int j = 0; j <= i_max; j++)
			{
				mat[i][j] = new_selection[i][j];
			}
		}

		for (Probe p : probes) p.rotate90(i_max, j_max);
	}
}

class ClipboardMaterial implements Cloneable {
	double rho_n = 0;
	double rho_p = 0;
	Material m = new Material();
	
	public ClipboardMaterial() {};

	public ClipboardMaterial(Material m) {
		this.m.copyFrom(m);
	}
	
    public ClipboardMaterial(Simulation e, int i, int j) {
    	this(e.materials[i][j]);
    	rho_n = e.rho_n[i][j];
    	rho_p = e.rho_p[i][j];
    }
    
    public void paste(Simulation e, int i, int j) {
    	e.materials[i][j].copyFrom(m);
    	e.rho_n[i][j] = rho_n;
    	e.rho_p[i][j] = rho_p;
    }
	
	public void erase() {
		m.initialize();
		rho_n = 0;
		rho_p = 0;
	}
	
    @Override
    public ClipboardMaterial clone() {
        ClipboardMaterial mat = new ClipboardMaterial(m);
        mat.rho_n = rho_n;
        mat.rho_p = rho_p;
        return mat;
    }
    
    public void rotate90() {
		m.emf_direction += Math.PI/2.0;
    }
    
    public void flip_h() {
		m.emf_direction = Math.PI - m.emf_direction;
    }
    
    public void flip_v() {
		m.emf_direction = -m.emf_direction;
    }
}