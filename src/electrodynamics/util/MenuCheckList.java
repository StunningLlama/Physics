package electrodynamics.util;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.ButtonGroup;
import javax.swing.JMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;

public class MenuCheckList<T extends Enum<?>> implements ActionListener {

	public HashMap<T, JRadioButtonMenuItem> buttonmap = new HashMap<T, JRadioButtonMenuItem>();
	public HashMap<Integer, JRadioButtonMenuItem> buttonlist = new HashMap<Integer, JRadioButtonMenuItem>();
	public ButtonGroup buttongroup = new ButtonGroup();
	JMenu menu;
	
	public void initialize(T[] values, JMenu menu, ActionListener a, T default_option, T[] sep) {
		this.menu = menu;
		
		int i = 0;
		for (T b : values) {
			if (sep != null && Arrays.asList(sep).contains(b))
				menu.add(new JSeparator());
			
			JRadioButtonMenuItem button = new JRadioButtonMenuItem(b.toString());
			buttonmap.put(b, button);
			buttonlist.put(i, button);
			button.addActionListener(a);
			button.addActionListener(this);
			buttongroup.add(button);
			menu.add(button);
			
			i++;
		}
		
		if (default_option != null)
			buttongroup.setSelected(buttonmap.get(default_option).getModel(), true);
	}
	
	public void removeOption(T t) {
		menu.remove(buttonmap.get(t));
	}

	public T getOption() {
		for (T t : buttonmap.keySet()) {
			if (buttonmap.get(t).getModel() == buttongroup.getSelection())
				return t;
		}
		
		return null;
	}


	public void setOption(T t) {
		buttongroup.setSelected(buttonmap.get(t).getModel(), true);
	}
	
	public void setOption(int i) {
		buttongroup.setSelected(buttonlist.get(i).getModel(), true);
	}

	@Override
	public void actionPerformed(ActionEvent ev) {
		//for (T b : buttonmap.keySet())
		//	if (ev.getSource() == buttonmap.get(b))
		//		selected = b;
	}
	
	public T containsButton(JRadioButtonMenuItem src) {
		for (T b : buttonmap.keySet())
			if (src == buttonmap.get(b))
				return b;
		
		return null;
	}
}
