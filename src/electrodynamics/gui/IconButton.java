package electrodynamics.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.image.BufferedImage;

import javax.swing.ImageIcon;
import javax.swing.JButton;

public class IconButton extends JButton {

	private static final long serialVersionUID = 7046425063110958439L;

	BufferedImage icon;

	public IconButton(BufferedImage icon, int size) {
		super();

		this.icon = icon;
		this.setIcon(new ImageIcon(icon));
		setPreferredSize(new Dimension(size+8, size+8));
		setMaximumSize(new Dimension(size+8, size+8));
	}
	
	@Override
	public void updateUI() {
		super.updateUI();

		Color fc = this.getForeground();
		
		int fr = fc.getRed();
		int fg = fc.getGreen();
		int fb = fc.getBlue();
		
		if (icon != null) {
			BufferedImage image = icon.getSubimage(0, 0, icon.getWidth(), icon.getHeight());
			for(int y = 0; y < icon.getHeight(); y++)
				for(int x = 0; x < icon.getWidth(); x++)
				{
					int argb = icon.getRGB(x, y);

					int a = ((argb>>24)&255);
					
					int r = (int)(fr*a/256.0);
					int g = (int)(fg*a/256.0);
					int b = (int)(fb*a/256.0);
					if (r > 255)
						r = 255;
					if (g > 255)
						g = 255;
					if (b > 255)
						b = 255;
					
					image.setRGB(x, y, a<<24 | r << 16 | g << 8 | b);
				}

			setIcon(new ImageIcon(image));
		}
	}
}
