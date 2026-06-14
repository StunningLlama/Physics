package electrodynamics.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.function.Supplier;

import javax.swing.ButtonGroup;
import javax.swing.JMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;

public class MenuCheckList<T extends Enum<?>, Button extends JRadioButtonMenuItem> implements ActionListener {

	public HashMap<T, Button> buttonmap = new HashMap<T, Button>();
	public List<Button> buttonlist = new ArrayList<Button>();
	public List<T> optionlist = new ArrayList<T>();
	public ButtonGroup buttongroup = new ButtonGroup();
	public JMenu menu;
	public HashSet<T> separators;
	
	public void initialize(T[] values, JMenu menu, ActionListener a, T default_option, T[] sep, Supplier<Button> constructor) {
		this.menu = menu;
		this.separators = new HashSet<T>();
		if (sep != null)
			separators.addAll(Arrays.asList(sep));
		
		for (T b : values) {
			if (separators.contains(b))
				menu.add(new JSeparator());
			
			Button button = constructor.get();
			button.setText(b.toString());
			buttonmap.put(b, button);
			buttonlist.add(button);
			optionlist.add(b);
			button.setActionCommand(values[0].getClass().getName());
			button.addActionListener(a);
			button.addActionListener(this);
			buttongroup.add(button);
			menu.add(button);
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
	
	public Button getButton(T t) {
		return buttonmap.get(t);
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
