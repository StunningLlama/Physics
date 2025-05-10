// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;
import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import java.awt.Insets;
import java.awt.Font;

public class HelpDialog extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -58876812350742907L;
	
	private JPanel contentPane;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					HelpDialog frame = new HelpDialog();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
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
		txtpnThisProgramDemonstrates.setText("<h1>Semiconductor Simulation Program</h1>\n\n    <p>\n        This program demonstrates the behavior of semiconductor devices.\n        Press <strong>\"Open\"</strong> to load one of the demonstrations and use the mouse to interact.\n    </p>\n\n    <h2>Tools</h2>\n    <ul>\n        <li><strong>Interact:</strong> Click on a voltage source to set its strength.</li>\n        <li><strong>Draw:</strong> Add material to the field.</li>\n        <li><strong>Voltage:</strong> Add a voltage probe.</li>\n        <li><strong>Current:</strong> Click and drag to add a current probe that measures current across a wire.</li>\n        <li><strong>Ground:</strong> Specifies the point relative to which probes measure voltage (optional).</li>\n        <li><strong>Delete probe:</strong> Click to delete a probe.</li>\n        <li><strong>Replace:</strong> Draw over other materials.</li>\n        <li><strong>Line:</strong> Click and drag to make a line.</li>\n        <li><strong>Fill:</strong> Fill a region.</li>\n        <li><strong>Erase:</strong> Erase.</li>\n        <li><strong>Select:</strong> Click and drag to make a rectangular selection or move it.</li>\n        <li><strong>Select region:</strong> Click to select a contiguous region.</li>\n        <li><strong>Text:</strong> Click to place a text cursor and type text on the screen.</li>\n    </ul>\n\n    <h2>Controls</h2>\n    <ul>\n        <li><strong>P/Space:</strong> Pause & unpause</li>\n        <li><strong>F:</strong> Advance frame</li>\n        <li><strong>Q:</strong> Change brush shape</li>\n        <li><strong>C:</strong> Toggle material color</li>\n        <li><strong>V:</strong> Toggle vectors</li>\n        <li><strong>S:</strong> Toggle scalar colors</li>\n        <li><strong>T:</strong> Toggle tooltip</li>\n        <li><strong>G:</strong> Toggle text background</li>\n        <li><strong>Mouse wheel:</strong> Change brush size</li>\n        <li><strong>Shift:</strong> Draw straight lines</li>\n        <li><strong>Ctrl:</strong> Fill area</li>\n        <li><strong>Alt/Option:</strong> Pick material</li>\n        <li><strong>Ctrl-X:</strong> Cut</li>\n        <li><strong>Ctrl-C:</strong> Copy</li>\n        <li><strong>Ctrl-V:</strong> Paste</li>\n        <li><strong>Left mouse:</strong> Draw material</li>\n        <li><strong>Right mouse:</strong> Erase material</li>\n        <li><strong>Middle mouse:</strong> Pick material</li>\n    </ul>\n\n    <h2>Materials</h2>\n    <ul>\n        <li><strong>Voltage source:</strong> Generates a voltage that can be used to power circuits.</li>\n        <li><strong>Metal:</strong> Conducts electricity very well.</li>\n        <li><strong>Conductive metal:</strong> More conductive than normal metal.</li>\n        <li><strong>Resistive metal:</strong> Less conductive than normal metal.</li>\n        <li><strong>High workfunction metal:</strong> Forms an ohmic contact with p-type semiconductor.</li>\n        <li><strong>Low workfunction metal:</strong> Forms an ohmic contact with n-type semiconductor.</li>\n        <li><strong>Intrinsic semiconductor:</strong> Undoped, has equal number of electrons and holes.</li>\n        <li><strong>P-type semiconductor:</strong> Has more holes than electrons.</li>\n        <li><strong>N-type semiconductor:</strong> Has more electrons than holes.</li>\n        <li><strong>Heavily doped P-type semiconductor:</strong> Has many more holes than electrons.</li>\n        <li><strong>Heavily doped N-type semiconductor:</strong> Has many more electrons than holes.</li>\n        <li><strong>Lightly doped P-type semiconductor:</strong> Has slightly more holes than electrons.</li>\n        <li><strong>Lightly doped N-type semiconductor:</strong> Has slightly more electrons than holes.</li>\n        <li><strong>Dielectric:</strong> Has high relative permittivity, can be used for capacitors.</li>\n        <li><strong>Ferromagnet:</strong> Has high relative permeability, can be used for inductors.</li>\n        <li><strong>Decoration:</strong> Inert, can be used for text or circuit symbols.</li>\n        <li><strong>Vacuum:</strong> Empty space.</li>\n    </ul>\n\n    <p>\n        Copyright (c) 2025 Brandon Li<br>\n        <a href=\"mailto:brandonli.lex@gmail.com\">brandonli.lex@gmail.com</a>\n    </p>");
		scrollPane.setViewportView(txtpnThisProgramDemonstrates);
	}

}