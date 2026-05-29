package electrodynamics.units;

public class GaussianUnits extends UnitSystem {
	
	UnitSystem emu = new EMU();
	UnitSystem esu = new ESU();

	@Override
	public double toSI(double val, Quantity q) {
		return getSystem(q).toSI(val, q);
	}

	@Override
	public double fromSI(double val, Quantity q) {
		return getSystem(q).fromSI(val, q);
	}
	
	@Override
	public String getSymbol(Quantity q) {
		return getSystem(q).getSymbol(q);
	}
	

	@Override
	public double getDefaultMagnitude(Quantity q) {
		return getSystem(q).getDefaultMagnitude(q);
	}
	
	public UnitSystem getSystem(Quantity q) {
		if (q == Quantity.MAGNETIC_FIELD_STRENGTH || q == Quantity.MAGNETIC_FLUX || q == Quantity.MAGNETIC_FLUX_DENSITY)
			return emu;
		else
			return esu;
	}

	@Override
	public String getName() {
		return "CGS-Gaussian";
	}
}
