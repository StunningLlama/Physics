package electrodynamics.util;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.Supplier;

import javax.swing.ButtonGroup;
import javax.swing.JMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;

public class MenuCheckList<T extends Enum<?>, U extends JRadioButtonMenuItem> implements ActionListener {

	public HashMap<T, U> buttonmap = new HashMap<T, U>();
	public HashMap<Integer, U> buttonlist = new HashMap<Integer, U>();
	public ButtonGroup buttongroup = new ButtonGroup();
	JMenu menu;
	
	public void initialize(T[] values, JMenu menu, ActionListener a, T default_option, T[] sep, Supplier<U> constructor) {
		this.menu = menu;
		
		int i = 0;
		for (T b : values) {
			if (sep != null && Arrays.asList(sep).contains(b))
				menu.add(new JSeparator());
			
			U button = constructor.get();
			button.setText(b.toString());
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

	public void addOption(T t) {
		menu.add(buttonmap.get(t));
	}

	public T getOption() {
		for (T t : buttonmap.keySet()) {
			if (buttonmap.get(t).getModel() == buttongroup.getSelection())
				return t;
		}
		
		return null;
	}


	public void setOption(T t) {
		if (buttonmap.get(t) != null)
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
