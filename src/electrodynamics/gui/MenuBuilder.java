// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.gui;
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

                item.addActionListener((ev) -> {
            		new Thread(() -> readfileFunc.accept(entry)).start();
                });

                menu.add(item);
            }
        }
    }
}