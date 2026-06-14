package electrodynamics.gui;

import java.awt.event.MouseEvent;

import javax.swing.JRadioButtonMenuItem;

public class CustJRadioButtonMenuItem extends JRadioButtonMenuItem {

	private static final long serialVersionUID = 8471715336310341301L;

	public CustJRadioButtonMenuItem() {
		super();
	}
	
	public CustJRadioButtonMenuItem(String text) {
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