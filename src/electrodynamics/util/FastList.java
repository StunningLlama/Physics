// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

import java.util.ArrayList;

public class FastList<T> {
	
	ArrayList<T> list = new ArrayList<T>();

	public synchronized void add(T obj) {
		list.add(obj);
	}
	
	public synchronized void remove(int index) {
		if (index < list.size()) {
			list.set(index, list.get(list.size()-1));
			list.remove(list.size()-1);
		}
	}
	
	public synchronized void replace(int index, T obj) {
		if (index < list.size()) {
			list.set(index, obj);
		}
	}
	
	public T get(int index) {
		if (index < list.size())
			return list.get(index);
		return null;
	}
	
	public synchronized int size() {
		return list.size();
	}

	public synchronized void clear() {
		list.clear();
	}
}
