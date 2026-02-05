package electrodynamics.util;
import javax.swing.*;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.function.Consumer;

public class MenuBuilder {
	
    public static void addDirectoryToMenu(JMenu menu, File directory, String extension, Consumer<File> readfileFunc) {
        if (directory == null || !directory.isDirectory()) {
        	return;
        }

        File[] entries = directory.listFiles(new FilenameFilter() {
        	public boolean accept(File dir, String name) {
        		return new File(dir, name).isDirectory() || name.toLowerCase().endsWith(extension);
        	}
        });
        if (entries == null) {
            return;
        }

        Arrays.sort(entries, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));

        for (File entry : entries) {
            if (entry.isDirectory()) {
                JMenu subMenu = new JMenu(entry.getName());
                menu.add(subMenu);
                addDirectoryToMenu(subMenu, entry, extension, readfileFunc); // recursion
            } else {
                JMenuItem item = new JMenuItem(entry.getName().split("\\" + extension)[0]);

                item.addActionListener(ev -> {
            		SwingUtilities.invokeLater(() -> {
            			readfileFunc.accept(entry);
            		});
                });

                menu.add(item);
            }
        }
    }
}