// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import com.codedisaster.steamworks.SteamUGCDetails;

import electrodynamics.Simulation;
import electrodynamics.Steam;

public class SteamDownloadUI extends JFrame implements ActionListener, ListSelectionListener {

	private static final long serialVersionUID = 1L;
	Simulation e;
	private JPanel contentPane;
	
	public JList<String> list;
	public List<SteamUGCDetails> details;
	private JButton btn_cancel;
	private JButton btn_load;
	private JTextPane textPane;
	private JPanel panel;
	private JPanel panel_1;
	private JPanel panel_2;
	private JPanel panel_3;
	private JScrollPane scrollPane_1;
	private PreviewCanvas panel_4;

	public SteamDownloadUI() {
		setTitle("Open workshop item");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setBounds(100, 100, 705, 595);
		contentPane = new JPanel();
		contentPane.setPreferredSize(new Dimension(650, 500));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new BorderLayout(0, 0));
		
		panel = new JPanel();
		panel.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel, BorderLayout.WEST);
		panel.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane = new JScrollPane();
		panel.add(scrollPane);
		
		list = new JList<>();
		list.setPreferredSize(new Dimension(225, 0));
		list.setMinimumSize(new Dimension(225, 0));
		scrollPane.setViewportView(list);
		
		panel_1 = new JPanel();
		panel_1.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel_1, BorderLayout.CENTER);
		panel_1.setLayout(new BorderLayout(0, 0));
		
		scrollPane_1 = new JScrollPane();
		panel_1.add(scrollPane_1);
		
		textPane = new JTextPane();
		scrollPane_1.setViewportView(textPane);
		
		panel_4 = new PreviewCanvas();
		panel_4.setPreferredSize(new Dimension(200, 200));
		panel_4.setMinimumSize(new Dimension(200, 200));
		panel_1.add(panel_4, BorderLayout.NORTH);
		panel_4.setLayout(null);
		
		panel_3 = new JPanel();
		contentPane.add(panel_3, BorderLayout.SOUTH);
		panel_3.setLayout(new BorderLayout(0, 0));
		
		panel_2 = new JPanel();
		panel_3.add(panel_2, BorderLayout.EAST);
		panel_2.setLayout(new BoxLayout(panel_2, BoxLayout.X_AXIS));
		
		btn_load = new JButton("Load");
		panel_2.add(btn_load);
		
		btn_cancel = new JButton("Cancel");
		panel_2.add(btn_cancel);
		
		pack();
		
		initialize();
	}
	
	public void initialize() {
		this.btn_cancel.addActionListener(this);
		this.btn_load.addActionListener(this);
		this.list.addListSelectionListener(this);
		setLocationRelativeTo(null);

		textPane.setContentType("text/html");
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == btn_cancel) {
			this.setVisible(false);
		} else if (e.getSource() == btn_load) {
			int index = list.getSelectedIndex();
			if (index != -1) {
				Steam.openUGC(details.get(index).getPublishedFileID());
				this.setVisible(false);
			}
		}
	}

	@Override
	public void valueChanged(ListSelectionEvent ev) {
		int index = list.getSelectedIndex();
		if (index != -1) {
			textPane.setText(details.get(index).getDescription());
			textPane.setCaretPosition(0);
			panel_4.bi = Steam.loadThumbnail(details.get(index).getPublishedFileID());
			panel_4.repaint();
		}
	}
	

	public class PreviewCanvas extends JPanel {
		
		private static final long serialVersionUID = -87273824460710054L;
		public int zoom_bound_x = 0;
		public int zoom_bound_y = 0;
		public int offset_x = 0;
		public int offset_y = 0;
		public BufferedImage bi = null;

		@Override
		public void paintComponent(Graphics real) {
			draw((Graphics2D)real);
		}

		public void draw(Graphics2D g) {
			g.setBackground(Color.BLACK);
			
			if (g.getClipBounds() != null)
				g.clearRect(0, 0, g.getClipBounds().width, g.getClipBounds().height);
			
			if (bi == null)
				return;

			int canvas_x = this.getWidth();
			int canvas_y = this.getHeight();
			
			int zoom_i1 = 0;
			int zoom_j1 = 0;
			
			int zoom_i2 = bi.getWidth()-1;
			int zoom_j2 = bi.getHeight()-1;

			int xw = zoom_i2 - zoom_i1 + 1;
			int yw = zoom_j2 - zoom_j1 + 1;
			
			if (xw/(double)yw >= canvas_x/(double) canvas_y) {
				int dim2 = (int) (canvas_x*yw/(double)xw);
				zoom_bound_x = canvas_x - 1;
				zoom_bound_y = dim2 - 1;
				offset_x = 0;
				offset_y = (canvas_y - dim2)/2;
			} else {
				int dim2 = (int) (canvas_y*xw/(double)yw);
				zoom_bound_x = dim2 - 1;
				zoom_bound_y = canvas_y - 1;
				offset_x = (canvas_x - dim2)/2;
				offset_y = 0;
			}
			
			g.drawImage(bi, offset_x, offset_y, zoom_bound_x+1+offset_x, zoom_bound_y+1+offset_y, 
				zoom_i1, zoom_j1, (zoom_i2+1), (zoom_j2+1), SteamDownloadUI.this);
		}
	}
}
