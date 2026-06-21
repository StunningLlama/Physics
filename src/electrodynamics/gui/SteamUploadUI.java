package electrodynamics.gui;

import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JTextPane;
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

public class SteamUploadUI extends JFrame implements ActionListener, HyperlinkListener {
	private static final long serialVersionUID = -5574535390957201856L;
	public JTextField text_title;
	public JTextArea desc;
	public JButton btn_upload;
	public JButton btn_cancel;
	public int result = 0;
	private JPanel contentPane;
	public SteamUploadUI() {
		setBounds(100, 100, 438, 467);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));

		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		text_title = new JTextField();
		text_title.setText("(Simulation title)");
		text_title.setBounds(6, 6, 410, 26);
		getContentPane().add(text_title);
		text_title.setColumns(10);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 64, 410, 268);
		getContentPane().add(scrollPane);
		
		desc = new JTextArea();
		desc.setFont(new Font("SansSerif", Font.PLAIN, 11));
		scrollPane.setViewportView(desc);
		
		JLabel lblNewLabel = new JLabel("Description");
		lblNewLabel.setBounds(16, 44, 116, 16);
		getContentPane().add(lblNewLabel);
		
		btn_upload = new JButton("Upload");
		btn_upload.setBounds(182, 402, 117, 23);
		getContentPane().add(btn_upload);
		
		btn_cancel = new JButton("Cancel");
		btn_cancel.setBounds(299, 402, 117, 23);
		getContentPane().add(btn_cancel);
		
		JScrollPane scrollPane_1 = new JScrollPane();
		scrollPane_1.setBounds(6, 339, 410, 57);
		getContentPane().add(scrollPane_1);
		
		JTextPane txtpnbySubmittingThis = new JTextPane();
		txtpnbySubmittingThis.setEditable(false);
		txtpnbySubmittingThis.setContentType("text/html");
		txtpnbySubmittingThis.setText("<html>By submitting this item, you agree to the <a href=\"http://steamcommunity.com/sharedfiles/workshoplegalagreement\">workshop terms of service.</a></html>");
		scrollPane_1.setViewportView(txtpnbySubmittingThis);
		txtpnbySubmittingThis.addHyperlinkListener(this);
		
		btn_upload.addActionListener(this);
		btn_cancel.addActionListener(this);
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
