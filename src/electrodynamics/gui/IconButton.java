// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;

import javax.swing.ImageIcon;
import javax.swing.JButton;

public class IconButton extends JButton {

	private static final long serialVersionUID = 7046425063110958439L;

	BufferedImage icon;
	ImageIcon normal;
	ImageIcon highlighted;
	Color fc;
	Color bc;

	public IconButton(BufferedImage icon, int size) {
		super();

		this.icon = icon;
		this.setIcon(new ImageIcon(icon));
		setPreferredSize(new Dimension(size+8, size+8));
		setMaximumSize(new Dimension(size+8, size+8));
		
		updateUI();
	}

	@Override
	public void updateUI() {
		super.updateUI();

		fc = this.getForeground();
		bc = this.getBackground();

		if (icon != null) {
			normal = makeIcon(fc, Color.BLACK, true);
			highlighted = makeIcon(bc, fc, false);
			setIcon(normal);
		}
	}

	public void setHighlighted(boolean value) {
		if (value) {
			setIcon(highlighted);
			this.setBackground(fc);
		}
		else {
			setIcon(normal);
			this.setBackground(bc);
		}
	}
	
	private ImageIcon makeIcon(Color fc, Color bc, boolean applyAlpha) {
		int fr = fc.getRed();
		int fg = fc.getGreen();
		int fb = fc.getBlue();
		
		int br = bc.getRed();
		int bg = bc.getGreen();
		int bb = bc.getBlue();

		BufferedImage image = new BufferedImage(icon.getWidth(), icon.getHeight(), icon.getType());
		Graphics gr = image.getGraphics();
		gr.drawImage(icon, 0, 0, null);
		for(int y = 0; y < icon.getHeight(); y++) {
			for(int x = 0; x < icon.getWidth(); x++)
			{
				int argb = icon.getRGB(x, y);

				int a = ((argb>>24)&255);
				double af = a/255.0;
				
				int r = (int)(fr*af+br*(1-af));
				int g = (int)(fg*af+bg*(1-af));
				int b = (int)(fb*af+bb*(1-af));
				if (r > 255)
					r = 255;
				if (g > 255)
					g = 255;
				if (b > 255)
					b = 255;
				
				if (!applyAlpha)
					a = 255;
				
				image.setRGB(x, y, a<<24 | r << 16 | g << 8 | b);
			}
		}
		
		return new ImageIcon(image);
	}
}
