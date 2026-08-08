// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.event.MouseEvent;

import javax.swing.JMenuItem;

public class CustJMenuItem extends JMenuItem {

	private static final long serialVersionUID = -8521591420235407202L;

	public CustJMenuItem() {
		super();
	}
	
	public CustJMenuItem(String text) {
		super(text);
	}

	@Override
	protected void processMouseEvent(MouseEvent evt) {
		if (evt.getID() == MouseEvent.MOUSE_RELEASED  && contains(evt.getPoint())) {
			doClick();
			setArmed(true);
		} else
			super.processMouseEvent(evt);
	}

}