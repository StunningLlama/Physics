package electrodynamics.gui;

import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JTextPane;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkEvent.EventType;
import javax.swing.event.HyperlinkListener;

import com.codedisaster.steamworks.SteamRemoteStorage.WorkshopFileType;

import electrodynamics.Steam;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.awt.Font;
import java.awt.BorderLayout;
import java.awt.Dimension;

public class SteamUploadUI extends JFrame implements ActionListener, HyperlinkListener {
	private static final long serialVersionUID = -5574535390957201856L;
	public JTextField text_title;
	public JTextArea desc;
	public JButton btn_upload;
	public JButton btn_cancel;
	public int result = 0;
	private JPanel contentPane;
	private JPanel panel;
	private JPanel panel_1;
	private JPanel panel_2;
	private JPanel panel_3;
	public SteamUploadUI() {
		setBounds(100, 100, 438, 467);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(new BorderLayout(0, 0));
		
		panel = new JPanel();
		panel.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.add(panel, BorderLayout.CENTER);
		panel.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane = new JScrollPane();
		panel.add(scrollPane);
		
		desc = new JTextArea();
		desc.setWrapStyleWord(true);
		desc.setLineWrap(true);
		desc.setFont(new Font("SansSerif", Font.PLAIN, 13));
		scrollPane.setViewportView(desc);
		
		JLabel lblNewLabel = new JLabel("Description");
		panel.add(lblNewLabel, BorderLayout.NORTH);
		
		panel_1 = new JPanel();
		panel_1.setBorder(new EmptyBorder(0, 5, 5, 5));
		contentPane.add(panel_1, BorderLayout.SOUTH);
		panel_1.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane_1 = new JScrollPane();
		scrollPane_1.setPreferredSize(new Dimension(10, 60));
		scrollPane_1.setMinimumSize(new Dimension(10, 60));
		panel_1.add(scrollPane_1, BorderLayout.CENTER);
		
		JTextPane txtpnbySubmittingThis = new JTextPane();
		txtpnbySubmittingThis.setEditable(false);
		txtpnbySubmittingThis.setContentType("text/html");
		txtpnbySubmittingThis.setText("<html>By submitting this item, you agree to the <a href=\"http://steamcommunity.com/sharedfiles/workshoplegalagreement\">workshop terms of service.</a></html>");
		scrollPane_1.setViewportView(txtpnbySubmittingThis);
		
		panel_2 = new JPanel();
		panel_1.add(panel_2, BorderLayout.SOUTH);
		
		btn_upload = new JButton("Upload");
		panel_2.add(btn_upload);
		SwingUtilities.getRootPane(this).setDefaultButton(btn_upload);
		
		btn_cancel = new JButton("Cancel");
		panel_2.add(btn_cancel);
		btn_cancel.addActionListener(this);
		
		panel_3 = new JPanel();
		panel_3.setBorder(new EmptyBorder(2, 2, 2, 2));
		contentPane.add(panel_3, BorderLayout.NORTH);
		panel_3.setLayout(new BorderLayout(0, 0));
		
		text_title = new JTextField();
		panel_3.add(text_title);
		text_title.setText("(Simulation title)");
		text_title.setColumns(10);
		
		btn_upload.addActionListener(this);
		txtpnbySubmittingThis.addHyperlinkListener(this);
	}
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == btn_upload) {

			Steam.ws_title = text_title.getText();
			Steam.ws_description = desc.getText();
			Steam.UGC.createItem(Steam.Utils.getAppID(), WorkshopFileType.Community);
			this.setVisible(false);
			
			result = 1;
		} else if (e.getSource() == btn_cancel) {
			this.setVisible(false);
			
			result = 0;
		}
		
		/*synchronized(this) {
		    this.notify();
		}*/
	}
	@Override
	public void hyperlinkUpdate(HyperlinkEvent ev) {
		if (ev.getEventType() == EventType.ACTIVATED) {
			try {
				java.awt.Desktop.getDesktop().browse(new URI("http://steamcommunity.com/sharedfiles/workshoplegalagreement"));
			} catch (IOException | URISyntaxException ex) {
				ex.printStackTrace();
			}
		}
	}
}
