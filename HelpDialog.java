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
		setBounds(100, 100, 378, 315);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		
		JPanel panel = new JPanel();
		contentPane.add(panel, BorderLayout.CENTER);
		panel.setLayout(new BorderLayout(0, 0));
		
		JScrollPane scrollPane = new JScrollPane();
		panel.add(scrollPane);
		
		JTextArea textPane = new JTextArea();
		textPane.setWrapStyleWord(true);
		textPane.setText("This program demonstrates Maxwell's equations in a variety of different scenarios. Press \"Open\" to load one of the demonstrations and use mouse to interact.\r\n\r\nMaterials:\r\n- Conductor: Realistic conducting material that follows Ohm's law.\r\n- EMF: Acts as either an AC or DC voltage source. Can be turned on or off by clicking.\r\n- Dielectric: Material with adjustable permittivity.\r\n- Paramagnet/Diamagnet: Material with adjustable permeability.\r\n- Switch: Conductor that can be turned into insulator by clicking.\r\n\r\nOptions:\r\n- Clear fields: Resets EM fields while leaving objects intact\r\n- Reset all: Resets everything\r\n- Scalar/Vector view: Determines what colors/arrows represent, respectively.\r\n- EMF: Switches between steady/alternating field.\r\n- Direction: Sets direction of EMF\r\n\r\nControls:\r\n- P: Pause/Unpause\r\n- F: Advance frame (Single time step)\r\n- C: Clear fields\r\n- R: Reset all\r\n- TAB: Change brush shape\r\n- Mouse wheel: Change brush size\r\n- SHIFT: Draw straight lines\r\n- Left mouse: Add objects\r\n- Right mouse: Delete objects\r\n\r\nCopyright (c) 2022 Brandon Li\r\nbrandonli.lex@gmail.com");
		textPane.setFont(new Font("SansSerif", Font.PLAIN, 13));
		textPane.setMargin(new Insets(4, 4, 4, 4));
		textPane.setLineWrap(true);
		textPane.setEditable(false);
		scrollPane.setViewportView(textPane);
		textPane.setCaretPosition(0);
	}

}
