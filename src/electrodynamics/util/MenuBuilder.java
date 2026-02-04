package electrodynamics.util;
import javax.swing.*;

import electrodynamics.Simulation;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;

public class MenuBuilder {
	
    public static void addDirectoryToMenu(JMenu menu, File directory, Simulation e) {
        if (directory == null || !directory.isDirectory()) {
        	return;
        }

        File[] entries = directory.listFiles(new FilenameFilter() {
        	public boolean accept(File dir, String name) {
        		return new File(dir, name).isDirectory() || name.toLowerCase().endsWith(".semisim");
        	}
        });
        if (entries == null) {
            return;
        }

        // Optional: sort alphabetically
        Arrays.sort(entries, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));

        for (File entry : entries) {
            if (entry.isDirectory()) {
                JMenu subMenu = new JMenu(entry.getName());
                menu.add(subMenu);
                addDirectoryToMenu(subMenu, entry, e); // recursion
            } else {
                JMenuItem item = new JMenuItem(entry.getName().split("\\.semisim")[0]);

                item.addActionListener(ev -> {
            		SwingUtilities.invokeLater(() -> {
            			e.savemanager.readfile(entry);
            		});
                });

                menu.add(item);
            }
        }
    }
}