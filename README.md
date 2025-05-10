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

Background
==========

I created this simulator because I wanted to get a deeper understanding of how semiconductors work. It's been my experience that there's been a lack of good simulations that demonstrate advanced topics in physics. There certaintly exists many educational physics simulations, but they're all either aimed at lower educational levels, or they are very restricted in how the user can interact with the system. The only examples I've seen of simulations that combine advanced topics with a generous amount of interactivity are the physics applets written by Paul Falstad, to whom I'm also grateful for looking over my project and helping to convert it to Javascript. I've tried to give users many different ways of interacting with the simulation. Circuits can be drawn with just a few clicks of the mouse, allowing users to easily experiment with their own circuits. There are also many different ways of visualizing the underlying physics that I've incorporated into the settings. Each one gives a different perspective on the physics that is happening. At the end of the day, I think the best way to use my program is to just start playing around with it.

Physics
=======

My program simulates Maxwell's equations on a two-dimensional grid. The electric field is tangent to the screen and the magnetic field points out of it. To evolve the E and B fields forward in time, I use Yee's method. On top of this, I've added in semiconductors, which have two kinds of charge carriers, electrons and holes. Both types experience electric and chemical forces that determine their motion. The charge carrier density is determined by the continuity equation, with an extra term that describes recombination. When all this is put together, the result is a simulation that demonstrates many important properties of semiconductors. These include:

*   PN junctions: The depletion region and built-in potential are clearly visible.
*   Metal-semiconductor junctions: Metals and semiconductors will form either an ohmic contact or Schottky junction depending on their workfunctions.
*   Field effect: Charge carriers will respond to an external electric field, either moving towards or away from it.
*   Recombination: When electrons and holes recombine, light is produced (the simulation assumes all semiconductors are direct bandgap)
*   Galvani potential: Two metals with different chemical potentials will acquire an electrical potential difference on contact.
*   Thermoelectricity: Current flowing between different materials can result in cooling of one material and heating of another.

Limitations
-----------

A simulation is always an incomplete representation of physical reality. My software is intended to be educational, and used as a tool for gaining intuition, but no more than that. It's good for demonstrating how basic semiconductor devices operate, but it leaves out many physical effects that become important when the parameters of the system are pushed to the extremes. Some of the inaccuracies include:

*   Metals are modelled as semiconductors with a huge equilibrium carrier concentration. Real metals contain electrons at all energy levels in the conduction band, which contribute to conductivity in different amounts.
*   Electrical breakdown, which is important for the functioning of certain devices, is not modeled.
*   Quantum tunelling, which is important at small length scales and is used in tunnel diodes, isn't modelled.
*   Velocity saturation, which limits the amount of current that can flow in a semiconductor, isn't properly modelled.
*   Fermi level pinning isn't taken into account. The properties of the metal-semiconductor junction are influenced more by surface effects than the metal's workfunction.
*   The Hall effect isn't simulated.
*   The model for carrier generation and recombination is a simplified version of those found in textbooks: different recombination mechanisms are combined into a single rate constant.
*   Finally, certain properties of materials differ from their real life counterparts. For example, the charge carrier mobility of the semiconductor material is about 1500x greater than that of Silicon. This was done to demonstrate the properties of semiconductors as clearly as possible.

Simulation details
==================

The main way to interact with circuits is to change the strength of voltage sources and turn switches on and off. The quickest way to get started is to load one of the examples, uncheck the pause button and click on one of the voltage sources. You can then adjust the voltage using a slider located on the right panel.

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
*   **Switch:** Conductivity can be switched on and off by the user.
*   **Metal:** Material that conducts electricity very well.
*   **Conductive metal:** More conductive than regular metal.
*   **Resistive metal:** Less conductive than regular metal.
*   **High workfunction metal:** Metal that forms an ohmic contact with p-type semiconductor.
*   **Low workfunction metal:** Metal that forms an ohmic contact with n-type semiconductor.
*   **Intrinsic semiconductor:** Undoped, with equal number of electrons and holes.
*   **P-type semiconductor:** Represents a semiconductor doped with holes.
*   **N-type semiconductor:** Represents a semiconductor doped with electrons.
*   **Heavily doped P-type semiconductor:** Has a large concentration of holes.
*   **Heavily doped N-type semiconductor:** Has a large concentration of electrons.
*   **Lightly doped P-type semiconductor:** Has a small concentration of holes.
*   **Lightly doped N-type semiconductor:** Has a small concentration of electrons.
*   **Dielectric:** Material with a large permittivity/dielectric constant.
*   **Ferromagnet:** Magnetic material with high relative permeability.
*   **Positive static charge:** Positively charged insulating material.
*   **Negative static charge:** Negatively charged insulating material.
*   **Decoration:** Used for text or circuit symbols, has no effect otherwise.
*   **Vacuum:** Empty space.

What do the colors mean?
------------------------

In general, the color red is associated with either holes or a positive charge. Blue represents electrons or negative charge. White means both electrons and holes exist a location. In the rest of the cases, yellow represents a positve quantity (eg. chemical potential or magnetic field), while cyan is negative. Finally, green is used for quantites that are always positive (eg. energy density). Note: Each material also has its own color which is unrelated to the aforementioned color scheme.

What do voltmeters actually measure?
------------------------------------

You might notice that the reading from a voltage probe doesn't match the electric potential Φ. In reality, voltmeters do not measure Φ but rather differences in electrochemical potential of charge carriers. Things get a bit trickier when we ask what the voltage is in a piece of semiconductor, becuase now there are multiple charge carriers! In this case we can try to define voltage as the reading we get when we stick a small metallic probe at a certain point. This can actually be performed in the simulation, and the result is that the electrochemical potential of the metal lies between that of electrons and holes, closer to whichever one has a larger density. I approximate this with a simple weighted average, the result of which is displayed on the voltage probe.

Why does the magnetic field vanish outside of circuits?
-------------------------------------------------------

Because the simulation is in 2D, circuits actually extend infinitely in the z-direction (out of the page), so current flowing through a closed circuit has the same effect as current flowing through a 3D solenoid. If you recall from E&M class, the magnetic field within an infinitely long solenoid is entirely contained within it. This is certainly a point of departure from how we expect circuits to behave. It means that each current loop has its own inductance, and trying to create "inductors" that behave like their 3d counterparts is quite tricky.

Copyright (c) 2025 Brandon Li  
[brandonli.lex@gmail.com](mailto:brandonli.lex@gmail.com)