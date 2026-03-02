# Installation Instructions and Troubleshooting

## Installation

1.  Make sure you have the latest version of Java installed. Java can be found [here](https://www.java.com/en/download/manual.jsp).
2.  Extract the contents of SemiSim.zip into a new folder.
3.  Double click SemiSim.jar to run it.
4.  If SemiSim crashes when loading a file:
    *   Make sure you have a 64-bit version of Java installed.

## Mac OS

1.  If you see "_SemiSim.jar cannot be opened because it is from an unidentified developer_":
    *   Right click "SemiSim.jar" and click "Open".
    *   Click "Open" again on the popup window.
2.  If you are unable to see and open files:
    *   Go into System Preferences → Security → Privacy → Full Disk Access.
    *   Add "/System/Library/CoreServices/Jar Launcher.app" to the list and give it disk access.

# Build Instructions

Clone the repository using `bash git clone https://github.com/StunningLlama/SemiSim.git`  
Note: Building requires JRE 1.8 and Maven for dependency management.

## Eclipse

1.  Go to File → Import → Existing Maven Projects
2.  Select the cloned repository folder
3.  Click Finish

# Introduction

Brandon's semiconductor simulator (SemiSim) is an educational tool made with the original purpose of helping its creator understand semiconductor devices. It is fully interactive, letting users draw circuits and create their own devices in a manner similar to painting software. There is a wide variety of different materials to choose from and many ways to visualize the electromagnetic phenomena associated with semiconductors. Users can either load one of the many premade simulations or create their own.

## Physics

The simulation is set on a two-dimensional grid over which it solves the 2D Maxwell equations. The electric field lies in the same plane as the surface of the screen whereas the magnetic field is perpendicular to it. On top of this, there are two types of charge carriers, electrons and holes, which feel electric and chemical forces that determine their dynamics. The simulation uses a FDTD (finite-difference time domain) scheme obtained by discretizing the Maxwell equations and the drift-diffusion equations and coupling the two together. This results in a simulation that demonstrates many of the important properties of semiconductors and semiconductor devices, such as:

*   PN junctions
*   Metal-semiconductor junctions
*   Field effect transistors
*   Galvani potential
*   Thermoelectricity

## Limitations

Because the simulation uses a very simplified model of how semiconductors work, it is not able capture some phenomena that show up in real systems. Things that the simulation cannot handle include:

*   Metal band structures
*   Electrical breakdown
*   Kinetic effects
*   Quantum tunelling
*   Velocity saturation
*   Fermi level pinning
*   Different recombination mechanisms
*   Surface effects

Finally, certain properties of materials differ from their real life counterparts (for example, the charge carrier mobility of the semiconductor material is about 1500x greater than that of Silicon).

# Simulation features

The interface consists of the simulation area which the user can interact with, simulation controls that allow various parameters to be adjusted, and a menu bar containing different settings and tools.

![Application](images/app.png)

Simulation area (purple), controls (blue), and menu bar (green).

The main way to interact with circuits is through voltage sources and switches. A good way to get started with this software is to load one of the many examples.

## View options

**Colors:** The value of a field is mapped to a color.

![Colors](images/colors.png)

**Contour:** Contour lines indicate a region on which the scalar field has a constant value.

![Contour lines](images/contour.png)

**Arrows:** The direction and brightness of arrows corresponds to the direction and magnitude of the vector field.

![Arrows](images/arrows.png)

**Lines:** The brightness and density of lines indicates the magnitude of the field.

![Field lines](images/lines.png)

**Dots:** The velocity of the dots is proportional to the value of the vector field.

![Dots](images/dots.png)

**Charge carriers:** Electrons and holes are represented as blue and red dots, respectively. Both dots move at the actual drift velocities of their corresponding charge carriers and their density accurately reflects the density of charge carriers.

![Charge carriers](images/carriers.png)

## Options

|     |     |
| --- | --- |
| Boundary condition | Choose between a boundary that absorbs outgoing radiation or a perfectly conductive boundary that reflects it. |
| Timestep | Sets the simulation timestep. The maximum timestep is determined by the CFL condition for the wave equation and diffusion equations for each charge carrier. |
| Sim steps/frame | Sets the number of iterations performed during each frame. Most of the examples require at least 10 steps/frame to run responsively. The maximum number depends on how good the user's computer is. |
| Set fields to zero | Sets all fields to their default values, leaving the materials unchanged. |
| Set image size | Sets the size of the simulation display area. |
| Advanced settings | Opens advanced simulation settings window, contains settings to modify physics and material constants. |

## Tools

|     |     |
| --- | --- |
| Interact | Allows user to control voltage sources and turn switches on and off by clicking. |
| Flashlight | Shines a light on the region under the cursor, creating pairs of electrons and holes. |
| Zoom | Click and drag to zoom into a region. Click to zoom out. |
| Draw | Adds material to the field. |
| Replace | Similar to the draw tool, but overwrites occupied areas. |
| Line | Draws a line of material. |
| Fill | Fills a region with a certain material, similar to the bucket tool. |
| Eraser | Erases material. |
| Select and move | Makes a rectangular selection which can be dragged around and moved. |
| Flood select | Selects a contiguous region, similar to the bucket tool. |
| Text | Place a text cursor allowing text to be typed on the screen. |
| Voltage probe | Adds a voltage probe that measures electrochemical potential at a certain point (See "What do voltmeters actually measure"). |
| Current probe | Click and drag to add a current probe that measures current across a wire. |
| Charge probe | Click and drag to add a charge probe that measures electrical charge within a given region. |
| Ground | Specifies the point relative to which probes measure voltage (optional). |
| Delete probe | Click to delete a probe. |
| Move label | Moves labels attached to probes. |
| Plot bands | Click and drag to plot energy bands along a line. The plot contains the conduction and valence energy band edges as well as the quasi-Fermi levels of electrons and holes.<br><br>![Band plot](images/bandplot.png)<br><br>Valence and conduction band energies shown as solid lines. Quasi-Fermi levels are dashed lines. |
| Plot scalar field | Click and drag to make a plot of the currently selected scalar field. |
| Plot carriers | Click an drag to make a logarithmic plot of the number density of electrons and holes. |
| Plot probe data | Selecting this will create a plot of probe measurements over time, similar to an oscilloscope. |

## Keyboard/Mouse Controls

|     |     |
| --- | --- |
| P or Space | Pause & unpause |
| F   | Advance frame |
| Q   | Change brush shape |
| Mouse wheel | Change brush size |
| Alt or Option | Pick material |
| Left mouse | Draw material |
| Right mouse | Erase material |
| Middle mouse | Pick material |
| R   | Record probe data (saves to probedata.txt) |
| \[  | Previous tool |
| \]  | Next tool |
| Ctrl-X | Cut a region made with the select tool |
| Ctrl-C | Copy a region made with the select tool |
| Ctrl-V | Paste the clipboard |
| Ctrl-R | Rotate the clipboard after using the paste command |
| Ctrl-F | Flip the clipboard horizontally after using the paste command |
| Ctrl-G | Flip the clipboard vertically |

## Materials

|     |     |
| --- | --- |
| Voltage source | Generates a voltage that can be used to power circuits. |
| AC voltage source | Voltage source that oscillates sinusoidally at a fixed frequency. |
| Current source | Approximation of an ideal current source. Used to generate a constant current density. |
| Switch | Conductivity can be switched on and off by the user. |
| Metal | Material that conducts electricity very well. |
| Conductive metal | More conductive than regular metal. |
| Resistive metal | Less conductive than regular metal. |
| High workfunction metal | Metal that forms an ohmic contact with p-type semiconductor. |
| Low workfunction metal | Metal that forms an ohmic contact with n-type semiconductor. |
| Intrinsic semiconductor | Undoped, with equal number of electrons and holes. |
| P-type semiconductor | Represents a semiconductor doped with holes. |
| N-type semiconductor | Represents a semiconductor doped with electrons. |
| Heavily doped P-type semiconductor | Has a large concentration of holes. |
| Heavily doped N-type semiconductor | Has a large concentration of electrons. |
| Lightly doped P-type semiconductor | Has a small concentration of holes. |
| Lightly doped N-type semiconductor | Has a small concentration of electrons. |
| Dielectric | Material with a large permittivity/dielectric constant. |
| Ferromagnet | Magnetic material with high relative permeability. |
| Positive static charge | Positively charged insulating material. |
| Negative static charge | Negatively charged insulating material. |
| Absorber | Absorbs incoming electromagnetic radiation very effectively. |
| Decoration | Used for text or circuit symbols, has no effect otherwise. |
| Vacuum | Empty space. |

![Palette](images/palette.png)

Colors of all the materials.

# Miscellaneous questions and answers

## Why does this simulation exist?

I created this simulator because I wanted to get a deeper understanding of how semiconductors work. It's been my experience that there's been a lack of good simulations that demonstrate advanced topics in physics. There certainly exists many educational physics simulations, but they're all either aimed at lower educational levels, or they are very restricted in how the user can interact with the system. The only examples I've seen of simulations that combine advanced topics with a generous amount of interactivity are the physics applets written by Paul Falstad, to whom I'm also grateful for looking over my project and helping to convert it to Javascript. I've tried to give users many different ways of interacting with the simulation. Circuits can be drawn with just a few clicks of the mouse, allowing users to easily experiment with their own circuits. There are also many different ways of visualizing the underlying physics that I've incorporated into the settings. Each one gives a different perspective on the physics that is happening.

## What do the colors mean?

In general, the color red is associated with either holes or a positive field. Blue represents electrons or negative fields. White means both electrons and holes exist a location. For magnetic fields, cyan/yellow represent fields going in to and out of the page, respectively. Finally, green is used for quantities that are always positive (eg. energy density). Note: Each material also has its own color which is unrelated to the aforementioned color scheme.

## What do voltmeters actually measure?

You might notice that the reading from a voltage probe doesn't match the electric potential Φ. In reality, voltmeters do not measure Φ but rather differences in electrochemical potential of charge carriers. Things get a bit trickier when we ask what the voltage is in a piece of semiconductor, because now there are multiple charge carriers! In this case we can try to define voltage as the reading we get when we stick a small metallic probe at a certain point. This can actually be performed in the simulation, and the result is that the electrochemical potential of the metal lies between that of electrons and holes, closer to whichever one has a larger density. I approximate this with a simple weighted average, the result of which is displayed on the voltage probe.

## Why does the magnetic field vanish outside of circuits?

Because the simulation is in 2D, circuits actually extend infinitely in the z-direction (out of the page), so current flowing through a closed circuit has the same effect as current flowing through a 3D solenoid. If you recall from E&M class, the magnetic field within an infinitely long solenoid is entirely contained within it. This is certainly a point of departure from how we expect circuits to behave. It means that each current loop has its own inductance, and trying to create "inductors" that behave like their 3D counterparts is quite tricky.

Copyright (c) 2026 Brandon Li  
[brandonli.lex@gmail.com](mailto:brandonli.lex@gmail.com)