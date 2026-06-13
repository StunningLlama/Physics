// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

public class GeneralMaterialType {
	public MaterialType type = null;
	public int cust_id = -1;
	public String name = "";
	public static GeneralMaterialType EMPTY = new GeneralMaterialType(MaterialType.VACUUM);
	
	public GeneralMaterialType(MaterialType type) {
		this.type = type;
	}
	
	public GeneralMaterialType(int cust_id, String name) {
		this.cust_id = cust_id;
		this.name = name;
		this.type = MaterialType.CUSTOM;
	}
	
	public GeneralMaterialType(Material mat) {
		this.type = mat.type;
		this.name = mat.name;
		this.cust_id = mat.cust_id;
	}

	public String toString() {
		if (cust_id == -1)
			return type.getName();
		else
			return name + "*";
	}
	
	@Override
	public boolean equals(Object m) {
		if (!(m instanceof GeneralMaterialType)) return false;
		GeneralMaterialType mat = (GeneralMaterialType) m;
		return type == mat.type && cust_id == mat.cust_id;
	}
}