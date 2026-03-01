package electrodynamics;

public class ClipboardMaterial implements Cloneable {
	double rho_n = 0;
	double rho_p = 0;
	Material m = new Material();
	
	public ClipboardMaterial() {};

	public ClipboardMaterial(Material m) {
		this.m = m.clone();
	}
	
    public ClipboardMaterial(Simulation e, int i, int j) {
    	this(e.materials[i][j]);
    	rho_n = e.rho_n[i][j];
    	rho_p = e.rho_p[i][j];
    }
    
    public void paste(Simulation e, int i, int j) {
    	e.materials[i][j] = m.clone();
    	e.rho_n[i][j] = rho_n;
    	e.rho_p[i][j] = rho_p;
    }
	
	public void erase() {
		m.erase();
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
}