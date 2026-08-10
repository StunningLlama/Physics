// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.border.EmptyBorder;
import electrodynamics.Simulation;

public class HtmlBox extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	private JPanel contentPane;
	
	private JButton btn_close;
	public JTextPane textPane;
	public JScrollPane scrollPane;
	private JPanel panel_1;
	private JPanel panel_2;
	private JPanel panel_3;

	public HtmlBox() {
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 705, 595);
		contentPane = new JPanel();
		contentPane.setPreferredSize(new Dimension(650, 500));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new BorderLayout(0, 0));
		
		panel_1 = new JPanel();
		panel_1.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel_1, BorderLayout.CENTER);
		panel_1.setLayout(new BorderLayout(0, 0));
		
		scrollPane = new JScrollPane();
		panel_1.add(scrollPane);
		
		textPane = new JTextPane();
		scrollPane.setViewportView(textPane);
		
		panel_3 = new JPanel();
		contentPane.add(panel_3, BorderLayout.SOUTH);
		panel_3.setLayout(new BorderLayout(0, 0));
		
		panel_2 = new JPanel();
		panel_3.add(panel_2, BorderLayout.EAST);
		panel_2.setLayout(new BoxLayout(panel_2, BoxLayout.X_AXIS));
		
		btn_close = new JButton("Close");
		panel_2.add(btn_close);
		
		pack();
		
		initialize();
	}
	
	public void initialize() {
		this.btn_close.addActionListener(this);
		setLocationRelativeTo(null);

		textPane.setContentType("text/html");
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == btn_close) {
			this.setVisible(false);
			this.dispose();
		}
	}
}
