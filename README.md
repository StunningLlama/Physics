Semiconductor Simulation Program
================================

Running the Program
-------------------

1.  Make sure you have the latest version of Java installed.
2.  If SemiSim crashes when loading a file:
    *   Make sure you have a 64-bit version of Java installed.

Mac OS
------

1.  If you see "_SemiSim.jar cannot be opened because it is from an unidentified developer_":
    *   Right click "SemiSim.jar" and click **Open**.
    *   Click **Open** again on the popup window.
2.  If you are unable to see and open files:
    *   Go into **System Preferences → Security → Privacy → Full Disk Access**.
    *   Add `/System/Library/CoreServices/Jar Launcher.app` to the list and give it disk access.

* * *

This program demonstrates the behavior of semiconductor devices. Press **"Open"** to load one of the demonstrations and use the mouse to interact.

Tools
-----

*   **Interact:** Click on a voltage source to set its strength.
*   **Draw:** Add material to the field.
*   **Voltage:** Add a voltage probe.
*   **Current:** Click and drag to add a current probe that measures current across a wire.
*   **Ground:** Specifies the point relative to which probes measure voltage (optional).
*   **Delete probe:** Click to delete a probe.
*   **Replace:** Draw over other materials.
*   **Line:** Click and drag to make a line.
*   **Fill:** Fill a region.
*   **Erase:** Erase.
*   **Select:** Click and drag to make a rectangular selection or move it.
*   **Select region:** Click to select a contiguous region.
*   **Text:** Click to place a text cursor and type text on the screen.

Controls
--------

*   **P/Space:** Pause & unpause
*   **F:** Advance frame
*   **Q:** Change brush shape
*   **C:** Toggle material color
*   **V:** Toggle vectors
*   **S:** Toggle scalar colors
*   **T:** Toggle tooltip
*   **G:** Toggle text background
*   **Mouse wheel:** Change brush size
*   **Shift:** Draw straight lines
*   **Ctrl:** Fill area
*   **Alt/Option:** Pick material
*   **Ctrl-X:** Cut
*   **Ctrl-C:** Copy
*   **Ctrl-V:** Paste
*   **Left mouse:** Draw material
*   **Right mouse:** Erase material
*   **Middle mouse:** Pick material

Materials
---------

*   **Voltage source:** Generates a voltage that can be used to power circuits.
*   **Metal:** Conducts electricity very well.
*   **Conductive metal:** More conductive than normal metal.
*   **Resistive metal:** Less conductive than normal metal.
*   **High workfunction metal:** Forms an ohmic contact with p-type semiconductor.
*   **Low workfunction metal:** Forms an ohmic contact with n-type semiconductor.
*   **Intrinsic semiconductor:** Undoped, has equal number of electrons and holes.
*   **P-type semiconductor:** Has more holes than electrons.
*   **N-type semiconductor:** Has more electrons than holes.
*   **Heavily doped P-type semiconductor:** Has many more holes than electrons.
*   **Heavily doped N-type semiconductor:** Has many more electrons than holes.
*   **Lightly doped P-type semiconductor:** Has slightly more holes than electrons.
*   **Lightly doped N-type semiconductor:** Has slightly more electrons than holes.
*   **Dielectric:** Has high relative permittivity, can be used for capacitors.
*   **Ferromagnet:** Has high relative permeability, can be used for inductors.
*   **Decoration:** Inert, can be used for text or circuit symbols.
*   **Vacuum:** Empty space.

Copyright (c) 2025 Brandon Li  
[brandonli.lex@gmail.com](mailto:brandonli.lex@gmail.com)