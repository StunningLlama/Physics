// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.border.EmptyBorder;

public class HelpDialog extends JFrame {
	
	private static final long serialVersionUID = -58876812350742907L;
	
	private JPanel contentPane;

	public HelpDialog() {
		setTitle("Help dialog");
		setBounds(100, 100, 458, 388);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		
		JPanel panel = new JPanel();
		contentPane.add(panel, BorderLayout.CENTER);
		panel.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane = new JScrollPane();
		panel.add(scrollPane);
		
		JTextPane txtpnThisProgramDemonstrates = new JTextPane();
		txtpnThisProgramDemonstrates.setContentType("text/html");
		txtpnThisProgramDemonstrates.setText("<h1>Simulation details</h1>"
				+ "The main way to interact with circuits is to change the strength of voltage sources and turn switches on and off. The quickest way to get started is to load one of the examples,\r\n"
				+ "uncheck the pause button and click on one of the voltage sources. You can then adjust the voltage using a slider located on the right panel.\r\n"
				+ "<h2>Tools</h2>\r\n"
				+ "<ul>\r\n"
				+ "    <li><strong>Interact:</strong> Click on a voltage source to set its strength.</li>\r\n"
				+ "    <li><strong>Draw:</strong> Add material to the field.</li>\r\n"
				+ "    <li><strong>Voltage:</strong> Add a voltage probe.</li>\r\n"
				+ "    <li><strong>Current:</strong> Click and drag to add a current probe that measures current across a wire.\r\n"
				+ "    </li>\r\n"
				+ "    <li><strong>Ground:</strong> Specifies the point relative to which probes measure voltage (optional).</li>\r\n"
				+ "    <li><strong>Delete probe:</strong> Click to delete a probe.</li>\r\n"
				+ "    <li><strong>Replace:</strong> Draw over other materials.</li>\r\n"
				+ "    <li><strong>Line:</strong> Click and drag to make a line.</li>\r\n"
				+ "    <li><strong>Fill:</strong> Fill a region.</li>\r\n"
				+ "    <li><strong>Erase:</strong> Erase.</li>\r\n"
				+ "    <li><strong>Select:</strong> Click and drag to make a rectangular selection or move it.</li>\r\n"
				+ "    <li><strong>Select region:</strong> Click to select a contiguous region.</li>\r\n"
				+ "    <li><strong>Text:</strong> Click to place a text cursor and type text on the screen.</li>\r\n"
				+ "</ul>\r\n"
				+ "\r\n"
				+ "<h2>Controls</h2>\r\n"
				+ "<ul>\r\n"
				+ "    <li><strong>P/Space:</strong> Pause &amp; unpause</li>\r\n"
				+ "    <li><strong>F:</strong> Advance frame</li>\r\n"
				+ "    <li><strong>Q:</strong> Change brush shape</li>\r\n"
				+ "    <li><strong>C:</strong> Toggle material color</li>\r\n"
				+ "    <li><strong>V:</strong> Toggle vectors</li>\r\n"
				+ "    <li><strong>S:</strong> Toggle scalar colors</li>\r\n"
				+ "    <li><strong>T:</strong> Toggle tooltip</li>\r\n"
				+ "    <li><strong>G:</strong> Toggle text background</li>\r\n"
				+ "    <li><strong>Mouse wheel:</strong> Change brush size</li>\r\n"
				+ "    <li><strong>Shift:</strong> Draw straight lines</li>\r\n"
				+ "    <li><strong>Ctrl:</strong> Fill area</li>\r\n"
				+ "    <li><strong>Alt/Option:</strong> Pick material</li>\r\n"
				+ "    <li><strong>Ctrl-X:</strong> Cut</li>\r\n"
				+ "    <li><strong>Ctrl-C:</strong> Copy</li>\r\n"
				+ "    <li><strong>Ctrl-V:</strong> Paste</li>\r\n"
				+ "    <li><strong>Left mouse:</strong> Draw material</li>\r\n"
				+ "    <li><strong>Right mouse:</strong> Erase material</li>\r\n"
				+ "    <li><strong>Middle mouse:</strong> Pick material</li>\r\n"
				+ "</ul>\r\n"
				+ "\r\n"
				+ "<h2>Materials</h2>\r\n"
				+ "<ul>\r\n"
				+ "    <li><strong>Voltage source:</strong> Generates a voltage that can be used to power circuits.</li>\r\n"
				+ "    <li><strong>Switch:</strong> Conductivity can be switched on and off by the user.</li>\r\n"
				+ "    <li><strong>Metal:</strong> Material that conducts electricity very well.</li>\r\n"
				+ "    <li><strong>Conductive metal:</strong> More conductive than regular metal.</li>\r\n"
				+ "    <li><strong>Resistive metal:</strong> Less conductive than regular metal.</li>\r\n"
				+ "    <li><strong>High workfunction metal:</strong> Metal that forms an ohmic contact with p-type semiconductor.</li>\r\n"
				+ "    <li><strong>Low workfunction metal:</strong> Metal that forms an ohmic contact with n-type semiconductor.</li>\r\n"
				+ "    <li><strong>Intrinsic semiconductor:</strong> Undoped, with equal number of electrons and holes.</li>\r\n"
				+ "    <li><strong>P-type semiconductor:</strong> Represents a semiconductor doped with holes.</li>\r\n"
				+ "    <li><strong>N-type semiconductor:</strong> Represents a semiconductor doped with electrons.</li>\r\n"
				+ "    <li><strong>Heavily doped P-type semiconductor:</strong> Has a large concentration of holes.</li>\r\n"
				+ "    <li><strong>Heavily doped N-type semiconductor:</strong> Has a large concentration of electrons.</li>\r\n"
				+ "    <li><strong>Lightly doped P-type semiconductor:</strong> Has a small concentration of holes.</li>\r\n"
				+ "    <li><strong>Lightly doped N-type semiconductor:</strong> Has a small concentration of electrons.</li>\r\n"
				+ "    <li><strong>Dielectric:</strong> Material with a large permittivity/dielectric constant.</li>\r\n"
				+ "    <li><strong>Ferromagnet:</strong> Magnetic material with high relative permeability.</li>\r\n"
				+ "    <li><strong>Positive static charge:</strong> Positively charged insulating material.</li>\r\n"
				+ "    <li><strong>Negative static charge:</strong> Negatively charged insulating material.</li>\r\n"
				+ "    <li><strong>Decoration:</strong> Used for text or circuit symbols, has no effect otherwise.</li>\r\n"
				+ "    <li><strong>Vacuum:</strong> Empty space.</li>\r\n"
				+ "</ul>\r\n"
				+ "\r\n"
				+ "<h2>What do the colors mean?</h2>\r\n"
				+ "In general, the color red is associated with either holes or a positive charge. Blue represents electrons or negative charge. White means both electrons and holes exist a location.\r\n"
				+ "In the rest of the cases, yellow represents a positve quantity (eg. chemical potential or magnetic field), while cyan is negative. Finally, green is used for quantites that are always positive (eg. energy density).\r\n"
				+ "Note: Each material also has its own color which is unrelated to the aforementioned color scheme.\r\n"
				+ "\r\n"
				+ "<h2>What do voltmeters actually measure?</h2>\r\n"
				+ "You might notice that the reading from a voltage probe doesn't match the electric potential Φ. In reality, voltmeters do not measure Φ but rather differences in electrochemical potential of charge carriers.\r\n"
				+ "Things get a bit trickier when we ask what the voltage is in a piece of semiconductor, becuase now there are multiple charge carriers! In this case we can try to define voltage as the reading we get when we stick a small metallic\r\n"
				+ "probe at a certain point. This can actually be performed in the simulation, and the result is that the electrochemical potential of the metal lies between that of electrons and holes, closer to whichever one has a larger density.\r\n"
				+ "I approximate this with a simple weighted average, the result of which is displayed on the voltage probe.\r\n"
				+ "\r\n"
				+ "<h2>Why does the magnetic field vanish outside of circuits?</h2>\r\n"
				+ "Because the simulation is in 2D, circuits actually extend infinitely in the z-direction (out of the page), so current flowing through a closed circuit has the same effect as current flowing through a 3D solenoid.\r\n"
				+ "If you recall from E&amp;M class, the magnetic field within an infinitely long solenoid is entirely contained within it. This is certainly a point of departure from how we expect circuits to behave. It means that each current loop has its own inductance,\r\n"
				+ "and trying to create \"inductors\" that behave like their 3d counterparts is quite tricky.\r\n"
				+ "<p>\r\n"
				+ "    Copyright (c) 2025 Brandon Li<br>\r\n"
				+ "    <a href=\"mailto:brandonli.lex@gmail.com\">brandonli.lex@gmail.com</a>\r\n"
				+ "</p>");
		scrollPane.setViewportView(txtpnThisProgramDemonstrates);
	}

}