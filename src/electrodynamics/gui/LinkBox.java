// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.net.URISyntaxException;

import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.HyperlinkEvent.EventType;

public class LinkBox extends JEditorPane {
	private static final long serialVersionUID = -3582016301930083217L;

	public LinkBox(String text) {
	    JLabel label = new JLabel();
	    Font font = label.getFont();
	    Color color = label.getBackground();
	    
	    StringBuffer style = new StringBuffer("font-family:" + font.getFamily() + ";");
	    style.append("font-weight:" + (font.isBold() ? "bold" : "normal") + ";");
	    style.append("font-size:" + font.getSize() + "pt;");
	    style.append("background-color: rgb("+color.getRed()+","+color.getGreen()+","+color.getBlue()+");");
	    
		setEditable(false);
		setContentType("text/html");
		setText("<html><body style=\"" + style + "\">" + text + "</body></html>");
		setBorder(null);
		
		addHyperlinkListener(new HyperlinkListener() {
			@Override
			public void hyperlinkUpdate(HyperlinkEvent ev) {
				if (ev.getEventType() == EventType.ACTIVATED) {
					try {
						java.awt.Desktop.getDesktop().browse(ev.getURL().toURI());
					} catch (IOException | URISyntaxException ex) {
						ex.printStackTrace();
					}
				}
			}
		});
	}
}
