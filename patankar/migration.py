#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Migration Solver
===============================================================================
Implements a **1D finite-volume mass transfer solver (`senspatankar`)** for multilayer structures.
Uses a modified Patankar scheme with exact solutions to handle partitioning at interfaces.

**Main Components:**
- **`senspatankar`** (Main solver)
    - Computes time evolution of a migrating substance in a multilayer structure
    - Supports **Robin, impervious, and periodic** boundary conditions
    - Stores simulation results in `SensPatankarResult`
- **`SensPatankarResult`** (Stores simulation outputs)
    - Concentration profiles in packaging (`Cx`) and food (`CF`)
    - Time-dependent fluxes
    - Includes interpolation and visualization methods

**Integration with SFPPy Modules:**
- Requires `layer.py` to define multilayer structures.
- Uses `food.py` to set food contact conditions.
- Relies on `property.py` for migration parameters (D, K).
- Calls `geometry.py` when volume/surface area calculations are needed.

Example:
```python
from patankar.migration import senspatankar
solution = senspatankar(multilayer, medium)
solution.plotCF()
```


===============================================================================
Details
===============================================================================

This module provides a solver (``senspatankar``) to simulate in 1D the mass transfer of a substance
initially distributed into a multilayer packaging structure (``layer``) into a contacting medium (``foodlayer``).
It uses a finite-volume method adapted from the Patankar scheme to handle partition coefficients between all layers,
as well as between the food and the contact layer (food is on the left). The right boundary condition is assumed
impervious (no mass transfer at the right edge).

The numerical method has been published here:
    Nguyen, P.-M., Goujon, A., Sauvegrain, P. and Vitrac, O. (2013),
    A computer-aided methodology to design safe food packaging and related systems.
    AIChE J., 59: 1183-1212. https://doi.org/10.1002/aic.14056

The module offers :
    - methods to simulate mass transfer under various boudary conditions (Robin, impervious, periodic),
    - simulation chaining
    - result management (merging, edition...)
    - plotting and printing to disk capabilities


Classes
-------
- SensPatankarResult

Functions
---------
- senspatankar(multilayer, medium, t=None, autotime=True, timescale="sqrt", ntimes=1e4, RelTol=1e-4, AbsTol=1e-4)

Example
-------
```python

    from patankar.food import ethanol
    from patankar.layer import layer

    # Create medium and layers
    medium = ethanol()
    A = layer(layername="layer A")
    B = layer(layername="layer B")
    multilayer = A + B

    # Run solver
    sol = senspatankar(multilayer, medium)

    # Plot results
    sol.plotCF()
    sol.plotC()
```

@version: 1.40
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2022-01-17
@rev: 2025-03-26

"""
# Dependencies
import os
import shutil
import random
import re
from datetime import datetime
from copy import deepcopy as duplicate
# math libraries
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, coo_matrix
from scipy.interpolate import interp1d
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.optimize import minimize
# plot libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
# data libraries
import pandas as pd

# Local dependencies
from patankar.layer import layer, check_units, layerLink
from patankar.food import foodphysics,foodlayer
from patankar.useroverride import useroverride # useroverride is already an instance (not a class)

__all__ = ['CFSimulationContainer', 'Cprofile', 'PrintableFigure', 'SensPatankarResult', 'autoname', 'check_units', 'cleantex', 'colormap', 'compute_fc_profile_PBC', 'compute_fv_profile', 'create_plotmigration_widget', 'create_simulation_widget', 'custom_plt_figure', 'custom_plt_subplots', 'foodlayer', 'foodphysics', 'is_latex_available', 'is_valid_figure', 'layer', 'layerLink', 'print_figure', 'print_pdf', 'print_png', 'restartfile', 'restartfile_senspantakar', 'rgb', 'senspatankar', 'tooclear', 'useroverride']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.40"

# Plot configuration (preferred units)
plotconfig_default = {
    "tscale": 24.0 * 3600, # days used as time scale
    "tunit": "days",
    "lscale": 1e-6, # µm
    "lunit": "µm",
    "Cscale": 1.0,
    "Cunit": "a.u."
    }
_fig_metadata_atrr_ = "__filename__"
# %% Private functions and classes
def autoname(nchars=6, charset="a-zA-Z0-9"):
    """
    Generates a random simulation name.

    Parameters:
    - nchars (int): Number of characters in the name (default: 6).
    - charset (str): Character set pattern (e.g., "a-zA-Z0-9").

    Returns:
    - str: A randomly generated name.
    """

    # Expand regex-like charset pattern
    char_pool = []
    # Find all ranges (e.g., "a-z", "A-Z", "0-9")
    pattern = re.findall(r'([a-zA-Z0-9])\-([a-zA-Z0-9])', charset)
    for start, end in pattern:
        char_pool.extend(chr(c) for c in range(ord(start), ord(end) + 1))
    # Include any explicit characters (e.g., "ABC" in "ABC0-9")
    explicit_chars = re.sub(r'([a-zA-Z0-9])\-([a-zA-Z0-9])', '', charset)  # Remove ranges
    char_pool.extend(explicit_chars)
    # Remove duplicates and sort (just for readability)
    char_pool = sorted(set(char_pool))
    # Generate random name
    return ''.join(random.choices(char_pool, k=nchars))

def is_valid_figure(fig):
    """
    Checks if `fig` is a valid and open Matplotlib figure.

    Parameters:
    - fig: object to check

    Returns:
    - bool: True if `fig` is a valid, open Matplotlib figure.
    """
    return isinstance(fig, Figure) and hasattr(fig, 'canvas') and fig.canvas is not None


def _generate_figname(fig, extension):
    """
    Generate a clean filename based on metadata or current date/time.

    Parameters:
    - fig: Matplotlib figure object.
    - extension: File extension ('.pdf' or '.png').

    Returns:
    - str: Cleaned filename with correct extension.
    """
    # Try to retrieve the hidden filename metadata
    if hasattr(fig, _fig_metadata_atrr_):
        filename = getattr(fig, _fig_metadata_atrr_)
    else:
        # Default: Use date-time format if metadata is missing
        filename = "fig" + datetime.now().strftime("%Y%m%d_%H%M%S")
    # Clean filename (replace spaces, trim, remove special characters)
    filename = filename.strip().replace(" ", "_")
    # Ensure correct file extension
    if not filename.lower().endswith(extension):
        filename += extension
    return filename

def tooclear(color, threshold=0.6, correction=0.15):
    """
    Darkens a too-bright RGB(A) color tuple.

    Parameters:
    -----------
    color : tuple (3 or 4 elements)
        RGB or RGBA color in [0,1] range.
    threshold : float, optional (default=0.6)
        Grayscale threshold above which colors are considered too bright.
    correction : float, optional (default=0.15)
        Amount by which to darken too bright colors.

    Returns:
    --------
    tuple
        Adjusted RGB(A) color tuple with too bright colors darkened.

    Example:
    --------
    corrected_color = tooclear((0.9, 0.9, 0.7, 1.0))
    """
    if not isinstance(color, tuple) or len(color) not in [3, 4]:
        raise ValueError("Input must be an RGB or RGBA tuple.")
    rgb = color[:3]  # Extract RGB values
    # Compute grayscale brightness (mean of RGB channels)
    brightness = sum(rgb) / 3
    # Darken if brightness exceeds the threshold
    if brightness > threshold:
        rgb = tuple(max(0, c - correction) for c in rgb)
    return rgb + (color[3],) if len(color) == 4 else rgb  # Preserve alpha if present


def print_pdf(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a given figure as a PDF.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, PDF filename (auto-generated if empty).
    - destinationfolder: str, folder to save the file.
    - overwrite: bool, overwrite existing file.
    - dpi: int, resolution (default=300).
    """
    if not is_valid_figure(fig):
        print("no valid figure")
        return
    # Generate filename if not provided
    if not filename:
        filename = _generate_figname(fig, ".pdf")
    # add extension if missing
    if not filename.endswith(".pdf"):
        filename = filename+".pdf"
    # Ensure full path
    filename = os.path.join(destinationfolder, filename)
    # Prevent overwriting unless specified
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use overwrite=True to replace it.")
        return
    # Save figure as PDF
    fig.savefig(filename, format="pdf", dpi=dpi, bbox_inches="tight")
    print(f"Saved PDF: {filename}")


def print_png(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a given figure as a PNG.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, PNG filename (auto-generated if empty).
    - destinationfolder: str, folder to save the file.
    - overwrite: bool, overwrite existing file.
    - dpi: int, resolution (default=300).
    """
    if not is_valid_figure(fig):
        print("no valid figure")
        return
    # Generate filename if not provided
    if not filename:
        filename = _generate_figname(fig, ".png")
    # add extension if missing
    if not filename.endswith(".png"):
        filename = filename+".png"
    # Ensure full path
    filename = os.path.join(destinationfolder, filename)
    # Prevent overwriting unless specified
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use overwrite=True to replace it.")
        return
    # Save figure as PNG
    fig.savefig(filename, format="png", dpi=dpi, bbox_inches="tight")
    print(f"Saved PNG: {filename}")


def print_figure(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save the figure in both PDF and PNG formats.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, base filename (auto-generated if empty).
    - destinationfolder: str, folder to save the files.
    - overwrite: bool, overwrite existing files.
    - dpi: int, resolution (default=300).
    """
    if is_valid_figure(fig):
        print_pdf(fig, filename, destinationfolder, overwrite, dpi)
        print_png(fig, filename, destinationfolder, overwrite, dpi)
    else:
        print("no valid figure")


# Categorized colors with headers and spacing
COLOR_CATEGORIES = [
    ("White & Gray", ["White", "Snow", "Honeydew", "MintCream", "Azure", "AliceBlue", "GhostWhite", "WhiteSmoke",
                      "Seashell", "Beige", "OldLace", "FloralWhite", "Ivory", "AntiqueWhite", "Linen",
                      "LavenderBlush", "MistyRose", "Gray", "Gainsboro", "LightGray", "Silver", "DarkGray",
                      "DimGray", "LightSlateGray", "SlateGray", "DarkSlateGray", "Black"], 2),

    ("Red, Pink & Orange", ["Red", "LightSalmon", "Salmon", "DarkSalmon", "LightCoral", "IndianRed", "Crimson",
                            "FireBrick", "DarkRed", "", "Pink", "LightPink", "HotPink", "DeepPink", "PaleVioletRed",
                            "MediumVioletRed", "", "Orange", "DarkOrange", "Coral", "Tomato", "OrangeRed"], 1),

    ("Yellow & Brown", ["Yellow", "LightYellow", "LemonChiffon", "LightGoldenrodYellow", "PapayaWhip", "Moccasin",
                        "PeachPuff", "PaleGoldenrod", "Khaki", "DarkKhaki", "Gold", "", "Brown", "Cornsilk",
                        "BlanchedAlmond", "Bisque", "NavajoWhite", "Wheat", "BurlyWood", "Tan", "RosyBrown",
                        "SandyBrown", "Goldenrod", "DarkGoldenrod", "Peru", "Chocolate", "SaddleBrown",
                        "Sienna", "Maroon"], 2),

    ("Green", ["Green", "PaleGreen", "LightGreen", "YellowGreen", "GreenYellow", "Chartreuse", "LawnGreen", "Lime",
               "LimeGreen", "MediumSpringGreen", "SpringGreen", "MediumAquamarine", "Aquamarine", "LightSeaGreen",
               "MediumSeaGreen", "SeaGreen", "DarkSeaGreen", "ForestGreen", "DarkGreen", "OliveDrab", "Olive",
               "DarkOliveGreen", "Teal"], 0),

    ("Blue", ["Blue", "LightBlue", "PowderBlue", "PaleTurquoise", "Turquoise", "MediumTurquoise", "DarkTurquoise",
              "LightCyan", "Cyan", "Aqua", "DarkCyan", "CadetBlue", "LightSteelBlue", "SteelBlue", "LightSkyBlue",
              "SkyBlue", "DeepSkyBlue", "DodgerBlue", "CornflowerBlue", "RoyalBlue", "MediumBlue", "DarkBlue",
              "Navy", "MidnightBlue"], 0),

    ("Purple", ["Purple", "Lavender", "Thistle", "Plum", "Violet", "Orchid", "Fuchsia", "Magenta", "MediumOrchid",
                "MediumPurple", "Amethyst", "BlueViolet", "DarkViolet", "DarkOrchid", "DarkMagenta", "SlateBlue",
                "DarkSlateBlue", "MediumSlateBlue", "Indigo"], 0)
]
# Extract colors from Matplotlib
CSS_COLORS = {k.lower(): v for k, v in mcolors.CSS4_COLORS.items()}

def rgb():
    """Displays a categorized color chart with properly aligned headers."""
    ncols = len(COLOR_CATEGORIES)
    max_rows = max(len(colors) + spacing for _, colors, spacing in COLOR_CATEGORIES)
    fig, ax = plt.subplots(figsize=(ncols * 2.5, max_rows * 0.6))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    x_spacing = 1.8  # Horizontal spacing between columns
    y_spacing = 1.0  # Vertical spacing between color patches
    text_size = 13   # Increased text size by 50%
    for col_idx, (category, colors, extra_space) in enumerate(COLOR_CATEGORIES):
        y_pos = max_rows  # Start at the top
        ax.text(col_idx * x_spacing + (x_spacing - 0.2) / 2, y_pos + 1.2, category,
                fontsize=text_size + 2, fontweight='bold', ha='center')
        y_pos -= y_spacing  # Move down after title
        for color in colors:
            if color == "":  # Empty string is a spacer
                y_pos -= y_spacing * 0.5
                continue
            hexval = CSS_COLORS.get(color.lower(), "white")
            y_pos -= y_spacing  # Move down before drawing
            ax.add_patch(plt.Rectangle((col_idx * x_spacing, y_pos), x_spacing - 0.2, y_spacing - 0.2, facecolor=hexval))
            r, g, b = mcolors.to_rgb(hexval)
            brightness = (r + g + b) / 3
            text_color = 'white' if brightness < 0.5 else 'black'
            ax.text(col_idx * x_spacing + (x_spacing - 0.2) / 2, y_pos + y_spacing / 2, color, ha='center',
                    va='center', fontsize=text_size, color=text_color)
        y_pos -= extra_space * y_spacing
    ax.set_xlim(-0.5, ncols * x_spacing)
    ax.set_ylim(-0.5, max_rows * y_spacing + 2)
    plt.tight_layout()
    plt.show()

# return colormaps
def colormap(name="viridis", ncolors=16, tooclearflag=True, reverse=False):
    """
    Generates a list of `ncolors` colors from the specified colormap.

    Parameters:
    -----------
    name : str, optional (default="viridis")
        Name of the Matplotlib colormap to use.
    ncolors : int, optional (default=16)
        Number of colors to generate.
    tooclearflag : bool, optional (default=True)
        If True, applies `tooclear` function to adjust brightness.
    reverse : bool, optional (default=False)
        If True, reverses the colormap.

    Supported colormaps:
    --------------------
    - "viridis"
    - "jet"
    - "plasma"
    - "inferno"
    - "magma"
    - "cividis"
    - "turbo"
    - "coolwarm"
    - "spring"
    - "summer"
    - "autumn"
    - "winter"
    - "twilight"
    - "rainbow"
    - "hsv"

    Returns:
    --------
    list of tuples
        List of RGB(A) colors in [0,1] range.

    Raises:
    -------
    ValueError
        If the colormap name is not recognized.
    """
    cmap_name = name + "_r" if reverse else name  # Append "_r" to reverse colormap
    # Check if the colormap exists
    if cmap_name not in plt.colormaps():
        raise ValueError(f"Invalid colormap name '{cmap_name}'. Use one from: {list(plt.colormaps())}")

    cmap = plt.colormaps.get_cmap(cmap_name)  # Fetch the colormap
    colors = [cmap(i / (ncolors - 1)) for i in range(ncolors)]  # Normalize colors
    return [tooclear(c) if tooclearflag else c[:3] for c in colors]  # Apply tooclear if enabled

# Latex fixer on systems without LaTeX
def is_latex_available():
    """
    Check whether LaTeX is available in the system PATH.

    Returns:
        bool: True if LaTeX is found, False otherwise.
    """
    # 'latex' should be available in the PATH if a LaTeX distribution is installed.
    return shutil.which("latex") is not None

_LaTeXavailable = is_latex_available()

# Clean tex
def cleantex(text,islatexavailable=_LaTeXavailable):
    r"""
    Process a LaTeX string to guess the plain text by performing substitutions
    and removing formatting characters, while preserving inner content as much as possible.

    Replacements performed:
      - Replace '\,' with a space.
      - Replace '\frac{num}{denom}' with 'num/denom'.
      - Replace '\sum' with 'sum'.
      - Replace '\int' with 'int'.

    Afterwards, it removes:
      - Dollar signs ($$ and $) used as math delimiters.
      - Inline math delimiters: \( and \).
      - Display math delimiters: \[ and \].
      - Any remaining LaTeX commands (a backslash followed by letters).
      - Underscores (_).
      - Curly braces ({ and }).

    Parameters:
        text (str): The input text containing LaTeX code.

    Returns:
        str: The processed text with LaTeX formatting stripped out.
    """
    if islatexavailable:
        return text  # LaTeX is available, return the original text
    # Replace \, with space
    text = text.replace(r'\,', ' ')
    # Replace \frac{num}{denom} with num/denom
    text = re.sub(r'\\frac\{([^{}]+)\}\{([^{}]+)\}', r'\1/\2', text)
    # Replace \sum with "sum", \int with "int"
    text = re.sub(r'\\sum', 'sum', text)
    text = re.sub(r'\\int', 'int', text)
    # Remove dollar delimiters
    text = text.replace('$$', '')
    text = text.replace('$', '')
    # Remove inline math delimiters \( and \)
    text = re.sub(r'\\\(|\\\)', '', text)
    # Remove display math delimiters \[ and \]
    text = re.sub(r'\\\[|\\\]', '', text)
    # Remove any remaining LaTeX commands (backslash followed by letters)
    text = re.sub(r'\\[A-Za-z]+', '', text)
    # Remove underscores and curly braces
    text = text.replace('_', '')
    text = text.replace('{', '')
    text = text.replace('}', '')
    return text


# Define PrintableFigure class
class PrintableFigure(Figure):
    """Custom Figure class with show and print methods."""

    def show(self, display_mode=None):
        """
        Show figure based on the environment:
        - In Jupyter: Uses display(fig)
        - In scripts or GUI mode: Uses plt.show()
        """
        try:
            get_ipython  # Check if running inside Jupyter
            if display_mode is None or display_mode == "auto":
                display(self)  # Preferred for Jupyter
            else:
                super().show()  # Use Matplotlib’s default show
        except NameError:
            super().show()  # Use Matplotlib’s default show in scripts

    def print(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_figure(self, filename, destinationfolder, overwrite, dpi)

    def print_png(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_png(self, filename, destinationfolder, overwrite, dpi)

    def print_pdf(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_pdf(self, filename, destinationfolder, overwrite, dpi)

# ✅ Store original references
original_plt_figure = plt.figure
original_plt_subplots = plt.subplots

# ✅ Override `plt.figure()` to use PrintableFigure
def custom_plt_figure(*args, **kwargs):
    """Ensure all figures use PrintableFigure."""
    kwargs.setdefault("FigureClass", PrintableFigure)
    return original_plt_figure(*args, **kwargs)

# ✅ Override `plt.subplots()` to return PrintableFigure
def custom_plt_subplots(*args, **kwargs):
    """Ensure plt.subplots() returns a PrintableFigure."""
    fig, ax = original_plt_subplots(*args, **kwargs)  # Create normal figure
    fig.__class__ = PrintableFigure  # Explicitly convert to PrintableFigure
    return fig, ax

# ✅ Apply overrides
plt.figure = custom_plt_figure
plt.subplots = custom_plt_subplots
plt.FigureClass = PrintableFigure  # Ensure all figures default to PrintableFigure globally
plt.rcParams['figure.figsize'] = (8, 6)  # Optional: Default size

# %% Widget
def create_simulation_widget():
    """
    Creates a widget interface for launching a migration simulation.

    The UI consists of four horizontally arranged dropdowns representing:
      - Substance: selected from mysubstances.keys()
      - Contact (medium): selected from mycontacts.keys()
      - Packaging geometry: selected from mypackaging.keys()
      - Material (layer): selected from mymaterials.keys()

    These dropdowns are displayed in a row with operator labels:
        substance % medium << geom >> layers

    Below, a text field allows the user to give a simulation name (default "mig1"),
    and a "Launch Simulation" button triggers the simulation, which computes:

        result = substance % medium << geom >> layers >> medium

    The resulting simulation object is stored in the global dictionary mymigration.

    Returns:
      An ipywidgets.VBox instance containing the complete UI.
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError("ipywidgets and IPython are required for this widget interface.") from e

    from patankar.loadpubchem import migrant
    from patankar.geometry import Packaging3D
    import builtins
    # Ensure global dictionaries exist:
    if not hasattr(builtins, "mysubstances"):
        builtins.mysubstances = {}
    if not hasattr(builtins, "mycontacts"):
        builtins.mycontacts = {}
    if not hasattr(builtins, "mypackaging"):
        builtins.mypackaging = {}
    if not hasattr(builtins, "mymaterials"):
        builtins.mymaterials = {}
    if not hasattr(builtins, "mymigration"):
        builtins.mymigration = {}
    global mysubstances, mycontacts, mypackaging, mymaterials, mymigration
    mysubstances = builtins.mysubstances
    mycontacts = builtins.mycontacts
    mypackaging = builtins.mypackaging
    mymaterials = builtins.mymaterials
    mymigration = builtins.mymigration
    # flag for preheated gui interface (widgets should be initialized manually, instead of being empty)
    _preheatedGUI_ = hasattr(builtins, "_PREHEATED_") and getattr(builtins, "_PREHEATED_") is True

    # Create dropdowns for each parameter.
    substance_dropdown = widgets.Dropdown(
         options=list(mysubstances.keys()),
         description="Substance:",
         layout=widgets.Layout(width="200px")
    )
    medium_dropdown = widgets.Dropdown(
         options=list(mycontacts.keys()),
         description="Contact:",
         layout=widgets.Layout(width="200px")
    )
    geom_dropdown = widgets.Dropdown(
         options=list(mypackaging.keys()),
         description="Packaging:",
         layout=widgets.Layout(width="200px")
    )
    layers_dropdown = widgets.Dropdown(
         options=list(mymaterials.keys()),
         description="Material:",
         layout=widgets.Layout(width="200px")
    )

    # Create labels for operator symbols.
    op1 = widgets.Label(value=" % ")
    op2 = widgets.Label(value=" << ")
    op3 = widgets.Label(value=" >> ")

    # Arrange the four panels and operators horizontally.
    panels = widgets.HBox([substance_dropdown, op1, medium_dropdown, op2, geom_dropdown, op3, layers_dropdown])

    # Simulation name input.
    sim_name_text = widgets.Text(
         value="mig1",
         description="Sim Name:",
         layout=widgets.Layout(width="300px")
    )

    # Launch button.
    launch_button = widgets.Button(
         description="Launch Simulation",
         button_style="success",
         tooltip="Click to launch the simulation"
    )

    # Output area.
    sim_output = widgets.Output()

    # Callback to launch simulation.
    def launch_simulation(b):
         with sim_output:
             sim_output.clear_output()
             try:
                 substance = mysubstances[substance_dropdown.value]
                 medium = mycontacts[medium_dropdown.value]
                 geom = mypackaging[geom_dropdown.value]
                 layers = mymaterials[layers_dropdown.value]
             except KeyError as e:
                 print("Error: one of the selections is missing:", e)
                 return
             if not isinstance(substance,migrant):
                 raise TypeError(f"mysubstances['{substance_dropdown.value}'] must be a migrant not a {type(substance).__name__}")
             if not isinstance(medium,foodphysics):
                 raise TypeError(f"mycontacts['{medium_dropdown.value}'] must be a foodphysics not a {type(medium).__name__}")
             if not isinstance(geom,Packaging3D):
                 raise TypeError(f"mypackaging['{geom_dropdown.value}'] must be a Packaging3D not a {type(geom).__name__}")
             if not isinstance(layers,layer):
                 raise TypeError(f"mymaterials['{layers_dropdown.value}'] must be a layer not a {type(layers).__name__}")
             try:
                 # Compute the simulation result using overloaded operators.
                 result = substance % medium << geom >> layers >> medium
             except Exception as e:
                 print("Error during simulation computation:", e)
                 return
             sim_name = sim_name_text.value.strip()
             if not sim_name:
                 print("Please provide a simulation name.")
                 return
             mymigration[sim_name] = result
             print(f"Simulation '{sim_name}' launched successfully.")
             print("Result:", result)
             print("Current simulations:", list(mymigration.keys()))
    launch_button.on_click(launch_simulation)

    if _preheatedGUI_:
        launch_simulation(None) # We lauch the simulation manually

    # Arrange the complete UI.
    ui = widgets.VBox([panels, sim_name_text, launch_button, sim_output])
    return ui


def create_plotmigration_widget():
    """
    Creates a widget interface for managing simulation plots, point evaluations,
    tabulated values, and for saving simulation data. It uses the global dictionary
    mymigration (keys such as "mig1") to select a simulation instance.

    Panels:
     1. Plot Panel:
         - Dropdown to select simulation.
         - Checkbox (plotSML, default True).
         - Button "Plot CF" calls: simulation.plotCF(plotSML=plotSML)
         - Slider for nmax (5–30, default 15).
         - Button "Plot Cx" calls: simulation.plotCx(nmax=nmax)

     2. Point Evaluation Panel:
         - FloatText for ttarget.
         - Dropdown for tunit.
         - Button "Evaluate Point" computes and prints CF and CF/SML at ttarget.

     3. Tabulated Values Panel:
         - FloatTexts for trange_min and trange_max.
         - Slider for number of values (3–500, default 5).
         - Button "Generate Table" builds an interactive table (DataFrame) of t and CF.

     4. Save Panel:
         - Text field for filename (default: simulation key).
         - Text field for destination folder (default: os.getcwd()).
         - Two buttons: "Save CF" and "Save Cx" that trigger saving to CSV (first) and then to Excel.

    The simulation instance is retrieved as: simulation = mymigration[selected_key]
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display, clear_output
        import numpy as np
        import pandas as pd
        import os
    except ImportError as e:
        raise ImportError("ipywidgets, IPython, numpy and pandas are required for this widget.") from e

    # Global dictionary mymigration must exist.
    import builtins
    if not hasattr(builtins, "mymigration"):
        builtins.mymigration = {}
    global mymigration
    mymigration = builtins.mymigration
    # flag for preheated gui interface (widgets should be initialized manually, instead of being empty)
    _preheatedGUI_ = hasattr(builtins, "_PREHEATED_") and getattr(builtins, "_PREHEATED_") is True

    # ---- Interactive Table helper function ----
    def interactive_table(df):
        range_slider = widgets.IntRangeSlider(
             value=[0, len(df)],
             min=0,
             max=len(df),
             step=1,
             description='Rows:',
             continuous_update=False,
             layout=widgets.Layout(width='80%')
        )
        out = widgets.Output()
        def update_table(change):
             with out:
                  clear_output()
                  start, end = range_slider.value
                  display(df.iloc[start:end])
        range_slider.observe(update_table, names='value')
        update_table(None)
        return widgets.VBox([range_slider, out])

    # ---- Panel 1: Plot Panel ----
    sim_keys = list(mymigration.keys())
    if not sim_keys:
        sim_keys = ["<none>"]
    sim_dropdown = widgets.Dropdown(
         options=sim_keys,
         value=sim_keys[0],
         description="Sim:",
         layout=widgets.Layout(width="200px")
    )
    plotSML_checkbox = widgets.Checkbox(
         value=True,
         description="Plot SML"
    )
    plotCF_button = widgets.Button(
         description="Plot CF",
         button_style="info"
    )
    nmax_slider = widgets.IntSlider(
         value=15,
         min=5,
         max=30,
         step=1,
         description="nmax:",
         continuous_update=False,
         layout=widgets.Layout(width="200px")
    )
    plotCx_button = widgets.Button(
         description="Plot Cx",
         button_style="info"
    )
    plot_out = widgets.Output()
    plot_panel = widgets.VBox([
         widgets.HBox([sim_dropdown, plotSML_checkbox]),
         widgets.HBox([plotCF_button, nmax_slider, plotCx_button]),
         plot_out
    ])

    def on_plotCF(b):
        with plot_out:
            plot_out.clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            try:
                sim.plotCF(plotSML=plotSML_checkbox.value, subtitle=key)
            except Exception as e:
                print("Error in plotCF:", e)
    plotCF_button.on_click(on_plotCF)
    if _preheatedGUI_:
        on_plotCF(None) # we plot manually

    def on_plotCx(b):
        with plot_out:
            plot_out.clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            try:
                sim.plotCx(nmax=nmax_slider.value, subtitle=key)
            except Exception as e:
                print("Error in plotCx:", e)
    plotCx_button.on_click(on_plotCx)

    # ---- Panel 2: Point Evaluation Panel ----
    tunit_dropdown = widgets.Dropdown(
         options=["s", "min", "h", "days", "weeks", "months", "years"],
         value="days",
         description="tunit:"
    )
    ttarget_text = widgets.FloatText(
         value=1.0,
         description="ttarget:"
    )
    point_button = widgets.Button(
         description="Evaluate Point",
         button_style="success"
    )
    point_out = widgets.Output()
    point_panel = widgets.VBox([
         widgets.HBox([ttarget_text, tunit_dropdown]),
         point_button,
         point_out
    ])

    def on_eval_point(b):
        with point_out:
            point_out.clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            try:
                from patankar.layer import check_units
                tscale = check_units((1, tunit_dropdown.value))[0].item()
            except Exception as e:
                print("Error computing tscale:", e)
                return
            try:
                tmin = sim.t[0] / tscale
                tmax = sim.t[-1] / tscale
            except Exception:
                tmin, tmax = 0, 1
            ttarget = ttarget_text.value
            if ttarget < tmin or ttarget > tmax:
                print(f"Warning: ttarget ({ttarget}) is outside simulation time range [{tmin}, {tmax}].")
            try:
                CFtarget = sim.interp_CF(ttarget * tscale)
            except Exception as e:
                print("Error interpolating CF:", e)
                return
            # Concentration unit is assumed to be "a.u."
            Cunit = "a.u."
            try:
                SML = sim._SML
            except Exception:
                SML = None
            print(f"CF at t={ttarget} [{tunit_dropdown.value}]: {CFtarget:.6g} [{Cunit}]")
            if SML is not None:
                ratio = CFtarget / SML
                warn = ""
                if ratio > 10:
                    warn = " !!!"
                elif ratio > 1:
                    warn = " !"
                print(f"CF/SML at t={ttarget} [{tunit_dropdown.value}]: {ratio:.6g}{warn}")
    point_button.on_click(on_eval_point)

    # ---- Panel 3: Tabulated Values Panel ----
    trange_min_text = widgets.FloatText(
         value=0.0,
         description="t_min:"
    )
    trange_max_text = widgets.FloatText(
         value=10.0,
         description="t_max:"
    )
    num_values_slider = widgets.IntSlider(
         value=300,
         min=10,
         max=1000,
         step=10,
         description="n_vals:",
         continuous_update=False,
         layout=widgets.Layout(width="250px")
    )
    table_button = widgets.Button(
         description="Generate Table",
         button_style="success"
    )
    table_out = widgets.Output()
    table_panel = widgets.VBox([
         widgets.HBox([trange_min_text, trange_max_text, num_values_slider]),
         table_button,
         table_out
    ])

    def on_generate_table(b):
        with table_out:
            table_out.clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            try:
                from patankar.layer import check_units
                tscale = check_units((1, tunit_dropdown.value))[0].item()
            except Exception as e:
                print("Error computing tscale:", e)
                return
            try:
                tmin = sim.t[0] / tscale
                tmax = sim.t[-1] / tscale
            except Exception:
                tmin, tmax = 0, 1
            trange_min = trange_min_text.value if trange_min_text.value is not None else tmin
            trange_max = trange_max_text.value if trange_max_text.value is not None else tmax
            n_vals = num_values_slider.value
            trange = np.linspace(trange_min, trange_max, n_vals)
            try:
                CFrange = sim.interp_CF(trange * tscale)
            except Exception as e:
                print("Error interpolating CF for table:", e)
                return
            df = pd.DataFrame({
                f"Time ({tunit_dropdown.value})": trange,
                f"CF (a.u.)": CFrange
            })
            display(interactive_table(df))
    table_button.on_click(on_generate_table)

    # ---- Panel 4: Save Panel ----
    # Filename for saving (default: simulation key)
    filename_field = widgets.Text(
         value=sim_dropdown.value,
         description="Filename:",
         layout=widgets.Layout(width="300px")
    )
    # Destination folder (default: current working directory)
    destination_field = widgets.Text(
         value=os.getcwd(),
         description="Destination:",
         layout=widgets.Layout(width="400px")
    )
    save_cf_button = widgets.Button(
         description="Save CF",
         button_style="warning"
    )
    save_cx_button = widgets.Button(
         description="Save Cx",
         button_style="warning"
    )
    save_output = widgets.Output()

    def on_save_cf(b):
        with save_output:
            clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            fname = filename_field.value.strip() or key
            dest = destination_field.value.strip() or os.getcwd()
            # save CSV and XLS
            try:
                sim.save_as_csv_CF(fname + ".csv", destinationfolder=dest, overwrite=True)
                print("CF data saved as CSV.")
            except Exception as e:
                print("Error saving CF as CSV:", e)
            try:
                sim.save_as_excel_CF(fname + ".xlsx", destinationfolder=dest, overwrite=True)
                print("CF data saved as Excel.")
            except Exception as e:
                print("Error saving CF as Excel:", e)
            # print in PDF and PNG
            plt.ioff() # no interactive plots
            try:
                htmpCF = sim.plotCF(plotSML=plotSML_checkbox.value, subtitle=key, noshow=True)
                print_figure(htmpCF,fname, destinationfolder=dest, overwrite=True)
                plt.close(htmpCF)
                print("CF data printed as PNG and PDF.")
            except Exception as e:
                print("Error while printing CF in PNG or PDF", e)
            plt.ion() # reactivate interactive plots


    def on_save_cx(b):
        with save_output:
            clear_output()
            key = sim_dropdown.value
            if key not in mymigration:
                print("No simulation selected.")
                return
            sim = mymigration[key]
            fname = filename_field.value.strip() or key
            dest = destination_field.value.strip() or os.getcwd()
            try:
                sim.save_as_csv_Cx(fname + "_Cx.csv", destinationfolder=dest, overwrite=True,
                                   t=None, nmax=nmax_slider.value, long_format=False)
                print("Cx data saved as CSV.")
            except Exception as e:
                print("Error saving Cx as CSV:", e)
            try:
                sim.save_as_excel_Cx(fname + "_Cx.xlsx", destinationfolder=dest, overwrite=True,
                                     t=None, nmax=nmax_slider.value, long_format=False)
                print("Cx data saved as Excel.")
            except Exception as e:
                print("Error saving Cx as Excel:", e)
            # print in PDF and PNG
            plt.ioff() # no interactive plots
            try:
                htmpCx = sim.plotCx(nmax=nmax_slider.value, subtitle=key, noshow=True)
                print_figure(htmpCx,fname+"_Cx", destinationfolder=dest, overwrite=True)
                plt.close(htmpCx)
                print("Cx data printed as PNG and PDF.")
            except Exception as e:
                print("Error while printing Cx in PNG or PDF",e)
            plt.ion() # reactivate interactive plots


    save_cf_button.on_click(on_save_cf)
    save_cx_button.on_click(on_save_cx)

    save_panel = widgets.VBox([
         widgets.HTML("<h4>Save Simulation Data</h4>"),
         widgets.HBox([filename_field, destination_field]),
         widgets.HBox([save_cf_button, save_cx_button]),
         save_output
    ])

    # ---- Assemble Overall UI ----
    ui = widgets.VBox([
         widgets.HTML("<h3>Migration Simulation Evaluation</h3>"),
         plot_panel,
         widgets.HTML("<hr>"),
         point_panel,
         widgets.HTML("<hr>"),
         table_panel,
         widgets.HTML("<hr>"),
         save_panel
    ])
    return ui



# %% Generic Classes to manipulate results
class Cprofile:
    """
    A class representing a concentration profile C(x) for migration simulations.

    This class allows storing, interpolating, and analyzing the concentration of a
    migrating substance across a spatial domain.

    Attributes:
    -----------
    x : np.ndarray
        1D array of spatial positions.
    Cx : np.ndarray
        1D array of corresponding concentration values at `x`.

    Methods:
    --------
    interp(x_new)
        Interpolates the concentration at new spatial positions.
    integrate()
        Computes the integral of the concentration profile.
    mean_concentration()
        Computes the mean concentration over the spatial domain.
    find_indices_xrange(x_range)
        Returns indices where `x` falls within a specified range.
    find_indices_Cxrange(Cx_range)
        Returns indices where `Cx` falls within a specified concentration range.
    assign_values(indices, values)
        Assigns new concentration values at specified indices.

    Example:
    --------
    ```python
    x = np.linspace(0, 1, 10)
    Cx = np.exp(-x)
    profile = Cprofile(x, Cx)

    # Interpolating at new points
    new_x = np.linspace(0, 1, 50)
    interpolated_Cx = profile.interp(new_x)
    ```
    """

    def __init__(self, x=None, Cx=None):
        """Initialize the concentration profile Cx(x)."""
        if x is None or Cx is None:
            raise ValueError("Syntax: myprofile = Cprofile(x, Cx). Both x and Cx are mandatory.")
        self.x = np.array(x, dtype=float).reshape(-1)  # Ensure 1D NumPy array
        self.Cx = np.array(Cx, dtype=float).reshape(-1)  # Ensure 1D NumPy array
        # Check if x is strictly increasing
        if np.any(np.diff(self.x) <= 0):
            raise ValueError("x values must be strictly increasing.")
        # Create the interpolation function
        self._interp_func = interp1d(
            self.x, self.Cx, kind="linear", fill_value=0, bounds_error=False
        )

    def interp(self, x_new):
        """
        Interpolate concentration values at new x positions.

        Parameters:
            x_new (array-like): New positions where concentrations are needed.

        Returns:
            np.ndarray: Interpolated concentration values.
        """
        x_new = np.array(x_new, dtype=float)  # Ensure NumPy array
        return self._interp_func(x_new)

    def integrate(self):
        """
        Compute the integral of Cx over x using Simpson's rule.

        Returns:
            float: The integral ∫ Cx dx.
        """
        return simpson(self.Cx, self.x)

    def mean_concentration(self):
        """
        Compute the mean concentration using the integral.

        Returns:
            float: The mean value of Cx.
        """
        return self.integrate() / (self.x[-1] - self.x[0])

    def find_indices_xrange(self, x_range):
        """
        Find indices where x is within a specified range.

        Parameters:
            x_range (tuple): The (min, max) range of x.

        Returns:
            np.ndarray: Indices where x falls within the range.
        """
        xmin, xmax = x_range
        return np.where((self.x >= xmin) & (self.x <= xmax))[0]

    def find_indices_Cxrange(self, Cx_range=(0, np.inf)):
        """
        Find indices where Cx is within a specified range.

        Parameters:
            Cx_range (tuple): The (min, max) range of Cx.

        Returns:
            np.ndarray: Indices where Cx falls within the range.
        """
        Cmin, Cmax = Cx_range
        return np.where((self.Cx >= Cmin) & (self.Cx <= Cmax))[0]

    def assign_values(self, indices, values):
        """
        Assign new values to Cx at specified indices.

        Parameters:
            indices (array-like): Indices where values should be assigned.
            values (float or array-like): New values to assign.

        Raises:
            ValueError: If the number of values does not match the number of indices.
        """
        indices = np.array(indices, dtype=int)
        if np.isscalar(values):
            self.Cx[indices] = values  # Assign single value to all indices
        else:
            values = np.array(values, dtype=float)
            if values.shape[0] != indices.shape[0]:
                raise ValueError("Number of values must match the number of indices.")
            self.Cx[indices] = values

    def __repr__(self):
        """Representation of the profile."""
        stats_x = {
            "min": np.min(self.x),
            "max": np.max(self.x),
            "mean": np.mean(self.x),
            "median": np.median(self.x),
            "std": np.std(self.x),
        }
        stats_Cx = {
            "min": np.min(self.Cx),
            "max": np.max(self.Cx),
            "mean": np.mean(self.Cx),
            "median": np.median(self.Cx),
            "std": np.std(self.Cx),
        }

        print(
            f"Cprofile: {len(self.x)} points\n",
            f"x range: [{stats_x['min']:.4g}, {stats_x['max']:.4g}]\n",
            f"Cx range: [{stats_Cx['min']:.4g}, {stats_Cx['max']:.4g}]\n",
            f"x stats: mean={stats_x['mean']:.4g}, median={stats_x['median']:.4g}, std={stats_x['std']:.4g}\n",
            f"Cx stats: mean={stats_Cx['mean']:.4g}, median={stats_Cx['median']:.4g}, std={stats_Cx['std']:.4g}"
        )
        return str(self)

    def __str__(self):
        """Returns a formatted string representation of the profile."""
        return f"<{self.__class__.__name__}: including {len(self.x)} points>"


# -----------------------------------------------------
# Container for single simulations
# -----------------------------------------------------

class SensPatankarResult:
    """
    Container for the results of the 1D mass transfer simulation performed by ``senspatankar``.

    Attributes
    ----------
    ttarget : ndarray with shape (1,)
        target simulation time
        It is a duration not an absolute time.
    CFtarget : ndarray with shape (1,)
        CF value at ttarget
    Cxtarget : ndarray with shape (npoints,)
         Cx concentration profile at t=ttarget
    t : ndarray with shape (ntimes,)
        1D array of time points (in seconds) covering from 0 to 2*ttarget
        It is a duration not an absolute time.
    C : ndarray with shape (ntimes,)
        1D array of mean concentration in the packaging (averaged over all packaging nodes)
        at each time step. Shape: (ntimes,).
    CF : ndarray with shape (ntimes,)
        1D array of concentration in the food (left boundary) at each time step. Shape: (ntimes,).
    fc : ndarray with shape (ntimes,)
        1D array of the cumulative flux into the food. Shape: (ntimes,).
    f : ndarray with shape (ntimes,)
        1D array of the instantaneous flux into the food. Shape: (ntimes,).
    x : ndarray with shape (npoints,)
        1D array of the position coordinates of all packaging nodes (including sub-nodes).
        npoints = 3 * number of original FV elements (interfaces e and w are included).
    Cx : ndarray with shape (ntimes,npoints)
        2D array of the concentration profile across the packaging thickness for each time step.
        Shape: (ntimes, 3 * number_of_nodes). Each row corresponds to one time step.
    tC : ndarray with shape (ntimes,)
        1D array of the dimensionless time points
    C0eq : ndarray with shape (1,)
        Reference (equilibrium) concentration scaling factor.
    timebase : float
        Characteristic time scale (l_ref^2 / D_ref) used to normalize the solution.
    interp_CF : scipy.interpolate._interpolate.interp1d
        1D interpolant of CF vs time
    interp_Cx : scipy.interpolate._interpolate.interp1d
        1F interpolant of Cx vs time
    restart : restartfile_senspatankar object
        Restart object (see restartfile_senspatankar doc)

    ----- plotting parameters ------
    SML : None or float, optional
        SML value in final units (usually mg/kg, adapt the value accordingly)
    plotSML = True (default), bool, optional
        if True and SML is not None, plot SML limit as an horizontal line
    plotconfig : dict, optional
        Dictionary with plotting configuration, containing:
        - "tunit": Time unit label (e.g., 's').
        - "Cunit": Concentration unit label (e.g., 'mg/L').
        - "tscale": Time scaling factor.
        - "Cscale": Concentration scaling factor.

    """

    def __init__(self, name, description, ttarget, t, C, CF, fc, f, x, Cx, tC, C0eq, timebase,
                 restart,restart_unsecure,xi,Cxi,
                 SML = None, SMLunit=None, plotSML=True,
                 plotconfig=None, createcontainer=True, container=None, discrete=False):
        """Constructor for simulation results."""
        self.name = name
        self.description = description
        self.ttarget = ttarget
        self.t = t
        self.C = C
        self.CF = CF
        self.fc = fc
        self.f = f
        self.x = x
        self.Cx = Cx
        self.tC = tC
        self.C0eq = C0eq
        self.timebase = timebase
        self.discrete = discrete  # New flag for discrete data

        # Interpolation for CF and Cx
        self.interp_CF = interp1d(t, CF, kind="linear", fill_value="extrapolate")
        self.CFtarget = self.interp_CF(ttarget)
        self.interp_Cx = interp1d(t, Cx.T, kind="linear", axis=1, fill_value="extrapolate")
        self.Cxtarget = self.interp_Cx(ttarget)

        # Restart handling
        if xi is not None and Cxi is not None:
            Cxi_interp = interp1d(t, Cxi.T, kind="linear", axis=1, fill_value="extrapolate")
            Cxi_at_t = Cxi_interp(ttarget)
            restart.freezeCF(ttarget, self.CFtarget)
            restart.freezeCx(xi, Cxi_at_t)
        self.restart = restart # secure restart file (cannot be modified from outside)
        self.restart_unsecure = restart_unsecure # unsecure one (can be modified from outside)

        # Store state for simulation chaining
        self.savestate(self.restart.inputs["multilayer"], self.restart.inputs["medium"])

        # Plot configuration
        if plotconfig is not None:
            if not isinstance(plotconfig, dict):
                raise TypeError(f"plotconfig must be a dict not a {type(plotconfig).__name__}")
            self._plotconfig = {**plotconfig_default, **plotconfig}  # Merge without modifying plotconfig_default
        else:
            self._plotconfig = plotconfig_default.copy()  # Work with a copy to prevent accidental changes

        # SML_configuration
        if SML is None:
            if self.restart.inputs["medium"].hasSML:
                SML = self.restart.inputs["medium"].SML
                SMLunit = self.restart.inputs["medium"].SMLunit
        self._SML = SML
        self._SMLunit = self._plotconfig["Cunit"] if SMLunit is None else SMLunit
        self._plotSML = plotSML

        # Default container for results comparison
        if createcontainer:
            if container is None:
                self.comparison = CFSimulationContainer(name=name,plotconfig=self._plotconfig)
                currentname = "reference"
                color = "Crimson"
            elif isinstance(container, CFSimulationContainer):
                self.comparison = container
                currentname = name
                color = "Teal"
            else:
                raise TypeError(f"container must be a CFSimulationContainer, not {type(container).__name__}")
            self.comparison.add(self, label=currentname, color=color, linestyle="-", linewidth=2)

        # Distance pair
        self._distancepair = None


    def pseudoexperiment(self, npoints=25, std_relative=0.05, randomtime=False, autorecord=False, seed=None, t=None, CF=None, scale='linear'):
        """
        Generates discrete pseudo-experimental data from high-resolution simulated results.

        Parameters
        ----------
        npoints : int, optional
            Number of discrete time points to select (default: 25).
        std_relative : float, optional
            Relative standard deviation for added noise (default: 0.05).
        randomtime : bool, optional
            If True, picks random time points; otherwise, uses uniform spacing or a sqrt scale (default: False).
        autorecord : bool, optional
            If True, automatically adds the generated result to the container (default: False).
        seed : int, optional
            Random seed for reproducibility.
        t : list or np.ndarray, optional
            Specific time points to use instead of generated ones. If provided, `CF` must also be supplied.
        CF : list or np.ndarray, optional
            Specific CF values to use at the provided `t` time points. Must have the same length as `t`.
        scale : str, optional
            Determines how time points are distributed when `randomtime=False`:
            - "linear" (default): Uniformly spaced time points.
            - "sqrt": Time points are distributed more densely at the beginning using a square root scale.

        Returns
        -------
        SensPatankarResult
            A new SensPatankarResult object flagged as discrete.

        Raises
        ------
        ValueError
            If `t` and `CF` are provided but have mismatched lengths.
        """

        if seed is not None:
            np.random.seed(seed)

        if t is not None:
            t_discrete = np.array(t, dtype=float)
            if CF is None or len(CF) != len(t_discrete):
                raise ValueError("When providing t, CF values must be provided and have the same length.")
            CF_discrete_noisy = np.array(CF, dtype=float)
        else:
            if randomtime:
                t_discrete = np.sort(np.random.uniform(self.t.min(), self.t.max(), npoints))
            else:
                if scale == 'sqrt':
                    t_discrete = np.linspace(np.sqrt(self.t.min()), np.sqrt(self.t.max()), npoints) ** 2
                else:
                    t_discrete = np.linspace(self.t.min(), self.t.max(), npoints)

            CF_discrete = self.interp_CF(t_discrete)
            noise = np.random.normal(loc=0, scale=std_relative * CF_discrete)
            CF_discrete_noisy = CF_discrete + noise
            CF_discrete_noisy = np.clip(CF_discrete_noisy, a_min=0, a_max=None)

        discrete_result = SensPatankarResult(
            name=f"{self.name}_discrete",
            description=f"Discrete pseudo-experimental data from {self.name}",
            ttarget=self.ttarget,
            t=t_discrete,
            C=np.zeros_like(t_discrete),
            CF=CF_discrete_noisy,
            fc=np.zeros_like(t_discrete),
            f=np.zeros_like(t_discrete),
            x=self.x,
            Cx=np.zeros((len(t_discrete), len(self.x))),
            tC=self.tC,
            C0eq=self.C0eq,
            timebase=self.timebase,
            restart=self.restart,
            restart_unsecure=self.restart_unsecure,
            xi=None,
            Cxi=None,
            SML=self._SML,
            SMLunit=self._SMLunit,
            plotSML = self._plotSML,
            plotconfig=self._plotconfig,
            discrete=True
        )
        if autorecord:
            self.comparison.add(discrete_result, label="pseudo-experiment", color="black", marker='o', discrete=True)
        return discrete_result

    @property
    def currrentdistance(self):
        """returns the square distance to the last distance pair"""
        return self.distanceSq(self._distancepair) if self._distancepair is not None else None

    def __sub__(self, other):
        """Overloads the operator - for returning a square distance function"""
        return lambda: self.distanceSq(other)

    def distanceSq(self, other, std_relative=0.05, npoints=100, cum=True):
        """
        Compute the squared distance between two SensPatankarResult instances.

        Parameters
        ----------
        other : SensPatankarResult
            The other instance to compare against.
        std_relative : float, optional
            Relative standard deviation for normalization (default: 0.05).
        npoints : int, optional
            Number of points for interpolation if both are continuous (default: 100).
        cum : bool, optional
            If True, return the cumulative sum; otherwise, return pointwise values.

        Returns
        -------
        float or np.ndarray
            The squared normalized error.

        Raises
        ------
        TypeError
            If `other` is not an instance of SensPatankarResult.
        ValueError
            If the time ranges do not overlap or if discrete instances have different time points.
        """
        if not isinstance(other, SensPatankarResult):
            raise TypeError(f"other must be a SensPatankarResult not a {type(other).__name__}")

        # refresh
        self._distancepair = other # used for distance evaluation as self.currentdistance
        # Find common time range
        tmin, tmax = max(self.t.min(), other.t.min()), min(self.t.max(), other.t.max())
        if tmin >= tmax:
            raise ValueError("No overlapping time range between instances.")
        if not self.discrete and not other.discrete:
            # Case 1: Both are continuous
            t_common = np.linspace(tmin, tmax, npoints)
            CF_self = self.interp_CF(t_common)
            CF_other = other.interp_CF(t_common)
        elif self.discrete and not other.discrete:
            # Case 2: self is discrete, other is continuous
            t_common = self.t
            CF_self = self.CF
            CF_other = other.interp_CF(self.t)
        elif not self.discrete and other.discrete:
            # Case 3: self is continuous, other is discrete
            t_common = other.t
            CF_self = self.interp_CF(other.t)
            CF_other = other.CF
        else:
            # Case 4: Both are discrete
            if not np.array_equal(self.t, other.t):
                raise ValueError("Discrete instances must have the same time points.")
            t_common = self.t
            CF_self = self.CF
            CF_other = other.CF
        # Compute squared normalized error
        m = (CF_self + CF_other) / 2
        m[m == 0] = 1  # Avoid division by zero, results in zero error where both are zero
        e2 = ((CF_self - CF_other) / (m * std_relative)) ** 2
        return np.sum(e2) if cum else e2

    def fit(self,other,disp=True,std_relative=0.05,maxiter=100,xatol=1e-3,fatol=1e-3):
        """Fits simulation parameters D and k to fit a discrete CF data"""
        if not isinstance(other,SensPatankarResult):
            raise TypeError(f"other must be a SensPatankarResult not a {type(other).__name__}")
        if self.discrete:
            raise ValueError("the current instance contains discrete data, use it as other")
        if not other.discrete:
            raise ValueError("only discrete CF results can be fitted")
        # retrieve current Dlink and klink
        Dlink = self.restart_unsecure.inputs["multilayer"].Dlink
        klink = self.restart_unsecure.inputs["multilayer"].klink
        if Dlink is None and klink is None:
            raise ValueError("provide at least a Dlink or klink object")
        if Dlink is not None and not isinstance(Dlink,layerLink):
            raise TypeError(f"Dlink must be a layerLink not a {type(Dlink).__name__}")
        if klink is not None and not isinstance(klink,layerLink):
            raise TypeError(f"klink must be a layerLink not a {type(klink).__name__}")
        # options for the optimizer
        optimOptions = {"disp": disp, "maxiter": maxiter, "xatol": xatol, "fatol": fatol}
        # params is assembled by concatenating -log(Dlink.values) and log(klink.values)
        params_initial = np.concatenate((-np.log(Dlink.values),np.log(klink.values)))
        maskD = np.concatenate((np.ones(Dlink.nzlength, dtype=bool), np.zeros(klink.nzlength, dtype=bool)))
        maskk = np.concatenate((np.zeros(Dlink.nzlength, dtype=bool), np.ones(klink.nzlength, dtype=bool)))
        # distance criterion
        d2 = lambda: self.distanceSq(other, std_relative=0.05) # d2 = lambda: self - other works also
        def objective(params):
            """objective function, all parameters are passed via layerLink"""
            logD = params[maskD]
            logk = params[maskk]
            Dlink.values = np.exp(-logD)
            klink.values = np.exp(logk)
            self.rerun(name="optimizer",color="OrangeRed",linewidth=4)
            return d2()
        def callback(params):
            """Called at each iteration to display current values."""
            Dtmp, ktmp = np.exp(-params[maskD]), np.exp(params[maskk])
            print("Fitting Iteration:\n",f"D={Dtmp} [m²/s]\n",f"k={ktmp} [a.u.]\n")
        # do the optimization
        result = minimize(objective,
                          params_initial,
                          method='Nelder-Mead',
                          callback=callback,
                          options=optimOptions)
        # extract the solution, be sure it is updated
        Dlink.values, klink.values = np.exp(-result.x[maskD]), np.exp(result.x[maskk])
        return result


    def savestate(self,multilayer,medium):
        """Saves senspantankar inputs for simulation chaining"""
        self._lastmedium = medium
        self._lastmultilayer = multilayer
        self._isstatesaved = True

    def update(self, **kwargs):
        """
        Update modifiable parameters of the SensPatankarResult object.
        Parameters:
            - name (str): New name for the object.
            - description (str): New description.
            - tscale (float or tuple): Time scale (can be tuple like (1, "day")).
            - tunit (str): Time unit.
            - lscale (float or tuple): Length scale (can be tuple like (1e-6, "µm")).
            - lunit (str): Length unit.
            - Cscale (float or tuple): Concentration scale (can be tuple like (1, "a.u.")).
            - Cunit (str): Concentration unit.
        """
        def checkunits(value):
            """Helper function to handle unit conversion for scale/unit tuples."""
            if isinstance(value, tuple) and len(value) == 2:
                scale, unit = check_units(value)
                scale, unit = np.array(scale, dtype=float), str(unit)  # Ensure correct types
                return scale.item(), unit  # Convert numpy array to float
            elif isinstance(value, (int, float, np.ndarray)):
                value = np.array(value, dtype=float)  # Ensure float
                return value.item(), None  # Return as float with no unit change
            else:
                raise ValueError(f"Invalid value for scale/unit: {value}")

        # Update `name` and `description` if provided
        if "name" in kwargs:
            self.name = str(kwargs["name"])
        if "description" in kwargs:
            self.description = str(kwargs["description"])
        # Update `_plotconfig` parameters
        for key in ["tscale", "tunit", "lscale", "lunit", "Cscale", "Cunit"]:
            if key in kwargs:
                value = kwargs[key]
                if key in ["tscale", "lscale", "Cscale"]:
                    value, unit = checkunits(value)  # Process unit conversion
                    self._plotconfig[key] = value
                    if unit is not None:
                        self._plotconfig[key.replace("scale", "unit")] = unit  # Ensure unit consistency
                else:
                    self._plotconfig[key] = str(value)  # Convert unit strings directly
        return self  # Return self for method chaining if needed


    def rerun(self,name=None,color=None,linestyle=None,linewidth=None, container=None, **kwargs):
        """
        Rerun the simulation (while keeping everything unchanged)
            This function is intended to be used with layerLinks for updating internally the parameters.
            R.rerun() stores the updated simulation results in R
            Rupdate = R.rerun() returns a copy of R while updating R

        note: Use R.resume() to resume/continue a simulation not rerun, to be used for sensitivity analysis/fitting.
        """
        F = self._lastmedium
        P = self._lastmultilayer
        if not isinstance(F, foodphysics):
            raise TypeError(f"the current object is corrupted, _lastmedium is {type(self._lastmedium).__name__}")
        if not isinstance(P, layer):
            raise TypeError(f"the current object is corrupted, _lastmultilayer is {type(self._lastmultilayer).__name__}")
        container = self.comparison if container is None else container
        if not isinstance(container,CFSimulationContainer):
            raise TypeError(f"the container should be a CFSimulationContainer not a {type(CFSimulationContainer).__name__}")
        # rerun the simulation using unsecure restart data
        inputs = self.restart_unsecure.inputs # all previous inputs
        R = senspatankar(
                multilayer=inputs["multilayer"],
                medium=inputs["medium"],
                name=name if name is not None else inputs["name"],
                description=kwargs.get("description",inputs["description"]),
                t=kwargs.get("t",inputs["t"]),
                autotime=kwargs.get("autotime",inputs["autotime"]),
                timescale=useroverride("timescale",kwargs.get("timescale",inputs["timescale"]),valuelist=("linear","sqrt")),
                Cxprevious=inputs["Cxprevious"],
                ntimes=useroverride("ntimes",kwargs.get("ntimes",inputs["ntimes"]),valuemin=10,valuemax=20000),
                RelTol=useroverride("RelTol",kwargs.get("RelTol",inputs["RelTol"]),valuemin=1e-9,valuemax=1e-3),
                AbsTol=useroverride("AbsTol",kwargs.get("AbsTol",inputs["AbsTol"]),valuemin=1e-9,valuemax=1e-3),
                container=container)
        # Update numeric data in self whith those in R
        self.t = R.t
        self.C = R.C
        self.CF = R.CF
        self.fc = R.fc
        self.f = R.f
        self.x = R.x
        self.Cx = R.Cx
        self.tC = R.tC
        self.C0eq = R.C0eq
        self.timebase = R.timebase
        self.discrete = R.discrete
        self.interp_CF = R.interp_CF
        self.CFtarget = R.CFtarget
        self.interp_Cx = R.interp_Cx
        self.Cxtarget = R.Cxtarget
        # Update label, color, linestyle, linewidth for the new curve (-1: last in the container)
        # note if name already exists, the previous content is replaced
        self.comparison.update(-1, label=name, color=color, linestyle=linestyle, linewidth=linewidth)
        return self # for chaining


    def resume(self,t=None,**kwargs):
        """
        Resume simulation for a new duration (with all parameters are unchanged)

        For convenience user overrides are provided as:
            parameter = value
            with parameter = "name","description"..."RelTol","AbsTol" (see senspantankar)
        Use specifically:
            CF0 to assign a different concentration for the food
            Cx0 (Cprofile object) to assign a different concentration profile (not recommended)
            medium to set a different medium (food) in contact
        """

        # retrieve previous results
        previousCF = self.restart.CF # CF at at target
        previousCx = self.restart.Cprofile # corresponding profile
        previousmedium = self.restart.inputs["medium"].copy()
        previousmedium.CF0 = previousCF # we apply the concentration
        # CF override with CF=new value
        isCF0forced = "CF0" in kwargs
        newmedium = kwargs.get("medium",previousmedium)
        if isCF0forced:
            newCF0 = kwargs.get("CF0",previousCF)
            newmedium.CF0 = newCF0
        if t is None:
            ttarget = newmedium.get_param("contacttime",(10,"days"),acceptNone=False)
            t = 2*ttarget
        # Concentration profile override with Cx0=new profile
        newCx0 = kwargs.get("Cx0",previousCx)
        if not isinstance(newCx0,Cprofile):
            raise TypeError(f"Cx0 should be a Cprofile object not a {type(newCx0).__name__}")

        # extend the existing solution
        inputs = self.restart.inputs # all previous inputs
        newsol = senspatankar(
                multilayer=inputs["multilayer"],
                medium=newmedium,
                name=kwargs.get("name",inputs["name"]),
                description=kwargs.get("description",inputs["description"]),
                t=t,
                autotime=kwargs.get("autotime",inputs["autotime"]),
                timescale=useroverride("timescale",kwargs.get("timescale",inputs["timescale"]),valuelist=("linear","sqrt")),
                Cxprevious=newCx0,
                ntimes=useroverride("ntimes",kwargs.get("ntimes",inputs["ntimes"]),valuemin=10,valuemax=20000),
                RelTol=useroverride("RelTol",kwargs.get("RelTol",inputs["RelTol"]),valuemin=1e-9,valuemax=1e-3),
                AbsTol=useroverride("AbsTol",kwargs.get("AbsTol",inputs["AbsTol"]),valuemin=1e-9,valuemax=1e-3)
                )
        return newsol


    def copy(self):
        """
        Creates a deep copy of the current SensPatankarResult instance.

        Returns
        -------
        SensPatankarResult
            A new instance with identical attributes as the original.
        """
        return SensPatankarResult(
            name=self.name,
            description=self.description,
            ttarget=self.ttarget,
            t=self.t.copy(),
            C=self.C.copy(),
            CF=self.CF.copy(),
            fc=self.fc.copy(),
            f=self.f.copy(),
            x=self.x.copy(),
            Cx=self.Cx.copy(),
            tC=self.tC.copy(),
            C0eq=self.C0eq,
            timebase=self.timebase,
            restart=self.restart,
            restart_unsecure=self.restart_unsecure,
            xi=None,
            Cxi=None,
            SML=self._SML,
            SMLunit=self._SMLunit,
            plotSML=self._plotSML,
            plotconfig=self._plotconfig,
            discrete=self.discrete
        )

    def chaining(self,multilayer,medium,**kwargs):
        sim = self.resume(multilayer=multilayer,medium=medium,**kwargs)
        medium.lastsimulation = sim # store the last simulation result in medium
        medium.lastinput = multilayer # store the last input (in medium)
        sim.savestate(multilayer,medium) # store store the inputs in sim for chaining
        return sim

    # overloading operation
    def __rshift__(self, medium):
        """Overloads >> to propagate migration to food."""
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics object not a {type(medium).__name__}")
        if not self._isstatesaved:
            raise RuntimeError("The previous inputs were not saved within the instance.")
        # we update the contact temperature (see example3)
        return self.chaining(medium>>self._lastmultilayer,medium,CF0=self.restart.CF)

    def __add__(self, other):
        """Concatenate two solutions"""
        if not isinstance(other, SensPatankarResult):
            raise TypeError("Can only add two SensPatankarResult objects")

        # Ensure compatibility of x-axis
        if not np.isclose(self.x[0], other.x[0]) or not np.isclose(self.x[-1], other.x[-1]):
            raise ValueError("Mismatch in x-axis boundaries between solutions")

        # Interpolate other.Cx onto self.x
        interp_Cx_other = interp1d(other.x, other.Cx.T, kind="linear", fill_value=0, axis=0)
        Cx_other_interp = interp_Cx_other(self.x).T  # Ensuring shape (ntimes, npoints)

        # Restrict times for valid merging
        valid_indices_self = self.t <= self.ttarget
        valid_indices_other = (other.t > 0) #& (other.t <= other.ttarget)
        t_self = self.t[valid_indices_self]
        t_other = other.t[valid_indices_other] + self.ttarget  # Shift time

        # Merge time arrays without duplicates
        t_merged = np.unique(np.concatenate((t_self, t_other)))
        tC_merged = np.unique(np.concatenate((self.tC[valid_indices_self], other.tC[valid_indices_other])))

        # Merge concentration-related attributes
        C_merged = np.concatenate((self.C[valid_indices_self], other.C[valid_indices_other]))
        CF_merged = np.concatenate((self.CF[valid_indices_self], other.CF[valid_indices_other]))
        fc_merged = np.concatenate((self.fc[valid_indices_self], other.fc[valid_indices_other]))
        f_merged = np.concatenate((self.f[valid_indices_self], other.f[valid_indices_other]))

        # Merge concentration profiles
        Cx_merged = np.vstack((self.Cx[valid_indices_self], Cx_other_interp[valid_indices_other]))

        # Merged description
        if self.description and other.description:
            merged_description = f"Merged: {self.description} & {other.description}"
        elif self.description:
            merged_description = self.description
        elif other.description:
            merged_description = other.description
        else:
            merged_description = ""

        # Create new instance with merged data
        merged_result = SensPatankarResult(
            name=f"{self.name} + {other.name}" if self.name!=other.name else self.name,
            description=merged_description,
            ttarget=self.ttarget + other.ttarget,
            t=t_merged,
            C=C_merged,
            CF=CF_merged,
            fc=fc_merged,
            f=f_merged,
            x=self.x,  # Keep self.x as reference
            Cx=Cx_merged,
            tC=tC_merged,
            C0eq=self.C0eq,  # Keep self.C0eq
            timebase=other.timebase,  # Take timebase from other
            restart=other.restart,  # Take restart from other (the last valid one)
            restart_unsecure=other.restart_unsecure,  # Take restart from other (the last valid one)
            xi=None,  # xi and Cxi values are available
            Cxi=None,  # only from a fresh simulation
            SML=other._SML if self._SML is None else self._SML,
            SMLunit=other._SMLunit if self._SMLunit is None else self._SMLunit,
            plotSML=other._plotSML if self._plotSML is None else self._plotSML,
            plotconfig=other._plotconfig if self._plotconfig is None else self._plotconfig,
            discrete=self.discrete or other.discrete
        )

        return merged_result

    def interpolate_CF(self, t, kind="linear", fill_value="extrapolate"):
        """
        Interpolates the concentration in the food (CF) at given time(s).

        Parameters
        ----------
        t : float, list, tuple, or ndarray
            Time(s) at which to interpolate CF values.
            - If a tuple, it should be (value or list, unit) and will be converted to SI.
            - If a scalar or list, it is assumed to be in SI units already.
        kind : str, optional
            Interpolation method. Default is "linear".
            Possible values:
            - "linear": Piecewise linear interpolation (default).
            - "nearest": Nearest-neighbor interpolation.
            - "zero": Zero-order spline interpolation.
            - "slinear", "quadratic", "cubic": Spline interpolations of various orders.
        fill_value : str or float, optional
            Specifies how to handle values outside the given range.
            - "extrapolate" (default): Extrapolates values beyond available data.
            - Any float: Uses a constant value for out-of-bounds interpolation.

        Returns
        -------
        ndarray
            Interpolated CF values at the requested time(s).
        """
        # Convert time input to SI units if provided as a tuple
        if isinstance(t, tuple):
            t, _ = check_units(t)  # Convert to numeric array

        # Ensure t is a NumPy array for vectorized operations
        t = np.atleast_1d(t)

        # Create the interpolant on demand with user-defined settings
        interp_function = interp1d(self.t, self.CF, kind=kind, fill_value=fill_value, bounds_error=False)

        # Return interpolated values
        return interp_function(t)


    def __repr__(self):
        ntimes = len(self.t)
        nx = self.Cx.shape[1] if self.Cx.ndim > 1 else len(self.x)
        tmin, tmax = self.t.min(), self.t.max()
        xmin, xmax = self.x.min(), self.x.max()

        print(f"SensPatankarResult: {self.name}\n"
              f"\t {self.description if self.description != '' else '<no description>'}\n"
              f"\t - with {ntimes} time steps\n",
              f"\t - with {nx} spatial points\n"
              f"\t - Time range: [{tmin:.2e}, {tmax:.2e}] s\n"
              f"\t - Position range: [{xmin:.2e}, {xmax:.2e}] m")

        return str(self)


    def __str__(self):
        return (f'<{self.__class__.__name__}:{self.name}: '
            f'CF({(self.ttarget / self._plotconfig["tscale"]).item():.4g} [{self._plotconfig["tunit"]}]) = '
            f'{(self.CFtarget / self._plotconfig["Cscale"]).item():.4g} [{self._plotconfig["Cunit"]}]>')



    def plotCF(self, t=None, trange=None, SML=None, SMLunit=None, plotSML=None, plotconfig=None, title=None, subtitle=None, noshow=False):
        """
        Plot the concentration in the food (CF) as a function of time.

        - If `self.discrete` is True, plots discrete points.
        - If `self.discrete` is False, plots a continuous curve.
        - Highlights the target time(s).

        Parameters
        ----------
        t : float, list, or None, optional
            Specific time(s) for which the concentration should be highlighted.
            If None, defaults to `ttarget`.
        trange : None, float, or list [t_min, t_max], optional
            If None, the full profile is shown.
            If a float, it is treated as an upper bound (lower bound assumed 0).
            If a list `[t_min, t_max]`, the profile is limited to that range.
        SML : None or float, optional
            SML value in final units (usually mg/kg, adapt the value accordingly)
        plotSML = True (default), bool, optional
            if True and SML is not None, plot SML limit as an horizontal line
        plotconfig : dict, optional
            Dictionary with plotting configuration, containing:
            - "tunit": Time unit label (e.g., 's').
            - "Cunit": Concentration unit label (e.g., 'mg/L').
            - "tscale": Time scaling factor.
            - "Cscale": Concentration scaling factor.
        title, subtitle: str, optional
        noshow : bool, optional
            if True, the figure is not shown (used for printing along with plt.ioff() and plt.ion() )
        """

        # plot configuration
        plt.rc('text', usetex=False) # Enable LaTeX formatting for Matplotlib
        if plotconfig is not None:
            if not isinstance(plotconfig, dict):
                raise TypeError(f"plotconfig must be a dict not a {type(plotconfig).__name__}")
            plotconfig = {**self._plotconfig, **plotconfig}  # Merge without modifying self._plotconfig
        else:
            plotconfig = self._plotconfig.copy()  # Work with a copy to prevent accidental changes

        # SML configuration
        SML = self._SML if SML is None else SML
        SMLunit = self._SMLunit if SMLunit is None else SMLunit
        plotSML = self._plotSML if plotSML is None else plotSML

        # we apply eventual user overrides
        plotconfig = useroverride("plotconfig",plotconfig,expected_type=dict)
        plotSML = useroverride("plotSML",plotSML,expected_type=bool)

        # Ensure t is a list (even if a single value is given)
        if t is None:
            t_values = [self.ttarget]
        elif isinstance(t, (int, float)):
            t_values = [t]
        elif isinstance(t, np.ndarray):
            t_values = t.flatten()
        elif isinstance(t, tuple):
            t_values = check_units(t)[0]
        else:
            t_values = np.array(t)  # Convert to array
        # Interpolate CF values at given times
        CF_t_values = self.interp_CF(t_values)
        # Handle trange selection
        if trange is None:
            t_plot = self.t
            CF_plot = self.CF
        else:
            # Convert trange to a valid range
            if isinstance(trange, (int, float)):
                trange = [0, trange]  # Assume lower bound is 0
            elif len(trange) != 2:
                raise ValueError("trange must be None, a single float (upper bound), or a list of two values [t_min, t_max]")
            # Validate range
            t_min, t_max = trange
            if t_min < self.t.min() or t_max > self.t.max():
                print("Warning: trange values are outside the available time range and may cause extrapolation.")
            # Generate time values within range
            mask = (self.t >= t_min) & (self.t <= t_max)
            t_plot = self.t[mask]
            CF_plot = self.CF[mask]
        # Set up colormap for multiple target values
        cmap = plt.get_cmap('viridis', len(t_values))
        norm = mcolors.Normalize(vmin=min(t_values), vmax=max(t_values))
        # Create the figure
        fig, ax = plt.subplots(figsize=(8, 6))
        # Plot behavior depends on whether data is discrete
        if self.discrete:
            ax.scatter(t_plot / plotconfig["tscale"], CF_plot / plotconfig["Cscale"],
                       color='b', label='Concentration in Food (Discrete)', marker='o', alpha=0.7)
        else:
            ax.plot(t_plot / plotconfig["tscale"], CF_plot / plotconfig["Cscale"],
                    label='Concentration in Food', color='b')
        # Highlight each target time
        for i, tC in enumerate(t_values):
            color = tooclear(cmap(norm(tC))) if len(t_values) > 1 else 'r'  # Use colormap only if multiple t values

            # Vertical and horizontal lines
            ax.axvline(tC / plotconfig["tscale"], color=color, linestyle='--', linewidth=1)
            ax.axhline(CF_t_values[i] / plotconfig["Cscale"], color=color, linestyle='--', linewidth=1)
            # Highlight points
            ax.scatter(tC / plotconfig["tscale"], CF_t_values[i] / plotconfig["Cscale"],
                       color=color, edgecolor='black', zorder=3, marker='D')
            # Annotate time
            ax.text(tC / plotconfig["tscale"], min(CF_plot) / plotconfig["Cscale"],
                    f'{(tC / plotconfig["tscale"]).item():.2f} {plotconfig["tunit"]}',
                    verticalalignment='bottom', horizontalalignment='right', rotation=90, fontsize=10, color=color)
            # Annotate concentration
            ax.text(min(t_plot) / plotconfig["tscale"], CF_t_values[i] / plotconfig["Cscale"],
                    f'{(CF_t_values[i] / plotconfig["Cscale"]).item():.2f} {plotconfig["Cunit"]}',
                    verticalalignment='bottom', horizontalalignment='left', fontsize=10, color=color)
        # add SML values (we assume that the units of SML are already final, no need to use plotconfig["Cscale"])
        if plotSML:
            SML = self._SML if SML is None else SML
            SMLunit = self._SMLunit if SMLunit is None else SMLunit
            if SML is not None:
                SMLunit = plotconfig["Cunit"] if SMLunit is None else SMLunit
                ax.axhline(SML, color="ForestGreen", linestyle=(0, (3, 5, 1, 5)),
                           linewidth=2, label=f"SML={SML:.4g} {SMLunit}")
                ax.text(self.ttarget / plotconfig["tscale"], SML, "specific migration limit", fontsize=8, color='ForestGreen', ha='center', va='bottom')
        # Labels and title
        ax.set_xlabel(f'Time [{plotconfig["tunit"]}]')
        ax.set_ylabel(f'Concentration in Food [{plotconfig["Cunit"]}]')
        if title:
            title_main = title
        else:
            title_main = "Concentration in Food vs. Time"
            if self.discrete:
                title_main += " (Discrete Data)"
        if subtitle:
            title_sub = subtitle
        else:
            title_sub = rf"$\bf{{{self.name}}}$" + (f": {self.description}" if self.description else "")
        ax.set_title(f"{title_main}\n{title_sub}", fontsize=10)
        #ax.text(0.5, 1.05, title_sub, fontsize=8, ha="center", va="bottom", transform=ax.transAxes)
        ax.legend()
        ax.grid(True)
        if not noshow:
            plt.show()
        # Store metadata
        setattr(fig, _fig_metadata_atrr_, f"pltCF_{self.name}")
        return fig


    def plotCx(self, t=None, nmax=15, plotconfig=None, title=None, subtitle=None, noshow=False):
        """
        Plot the concentration profiles (Cx) in the packaging vs. position (x) for different times,
        using a color gradient similar to Parula, based on time values (not index order).
        Additionally, highlight the concentration profile at `ttarget` with a thick black line.

        Parameters
        ----------
        t : list, array-like, or None, optional
            List of specific times to plot. Only valid values (inside self.t) are used.
            If None, time values are selected using sqrt-spaced distribution.
        nmax : int, optional
            Maximum number of profiles to plot. The default is 15.
        plotconfig : dict, optional
            Dictionary with plotting configuration, containing:
            - "tunit": Time unit label (e.g., 's').
            - "Cunit": Concentration unit label (e.g., 'mg/L').
            - "tscale": Time scaling factor.
            - "Cscale": Concentration scaling factor.
            title, subtitle: str, optional
            noshow : bool, optional
                if True, the figure is not shown (used for printing along with plt.ioff() and plt.ion() )
        """
        # short circuit
        if self.discrete:
            print("discrete SensPatankarResult instance does not contain profile data, nothing to plot.")
            return None

        # plot configuration
        plt.rc('text', usetex=False) # Disable LaTeX for Matplotlib
        if plotconfig is not None:
            if not isinstance(plotconfig, dict):
                raise TypeError(f"plotconfig must be a dict not a {type(plotconfig).__name__}")
            plotconfig = {**self._plotconfig, **plotconfig}  # Merge without modifying self._plotconfig
        else:
            plotconfig = self._plotconfig.copy()  # Work with a copy to prevent accidental changes

        # we apply eventual user overrides
        plotconfig = useroverride("plotconfig",plotconfig,expected_type=dict)
        nmax = useroverride("nmax",nmax,valuemin=3,valuemax=50)

        # Ensure time values are within the available time range
        if t is None:
            # Default: Select `nmax` time values using sqrt-spacing
            nt = len(self.t)
            if nt <= nmax:
                t_values = self.t
            else:
                sqrt_t = np.sqrt(self.t)
                sqrt_t_values = np.linspace(sqrt_t[0], sqrt_t[-1], nmax)
                t_values = sqrt_t_values**2
        else:
            # Use user-specified time values
            if isinstance(t,tuple):
                t_values = check_units(t)[0]
            else:
                t_values = np.array(t)
            # Keep only valid times inside `self.t`
            t_values = t_values[(t_values >= self.t.min()) & (t_values <= self.t.max())]
            if len(t_values) == 0:
                print("Warning: No valid time values found in the specified range.")
                return
            # If more than `nmax`, keep the first `nmax` values
            t_values = t_values[:nmax]
        # Normalize time for colormap (Ensure at least one valid value)
        norm = mcolors.Normalize(vmin=t_values.min()/plotconfig["tscale"],
                                 vmax=t_values.max()/plotconfig["tscale"]) if len(t_values) > 1 \
            else mcolors.Normalize(vmin=self.t.min()/plotconfig["tscale"],
                                   vmax=self.t.max()/plotconfig["tscale"])
        cmap = plt.get_cmap('viridis', nmax)  # 'viridis' is similar to Parula
        # new figure
        fig, ax = plt.subplots(figsize=(8, 6))  # Explicitly create a figure and axis
        # Plot all valid concentration profiles with time-based colormap
        for tC in t_values:
            C = self.interp_Cx(tC)
            color = tooclear(cmap(norm(tC/plotconfig["tscale"])))  # Get color from colormap
            ax.plot(self.x / plotconfig["lscale"], C / plotconfig["Cscale"],
                    color=color, alpha=0.9, label=f't={tC / plotconfig["tscale"]:.3g} {plotconfig["tunit"]}')
        # Highlight concentration profile at `ttarget`
        ax.plot(self.x / plotconfig["lscale"], self.Cxtarget / plotconfig["Cscale"], 'k-', linewidth=3,
                label=f't={self.ttarget[0] / plotconfig["tscale"]:.2g} {plotconfig["tunit"]} (target)')
        # Create ScalarMappable and add colorbar
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Needed for colorbar
        cbar = fig.colorbar(sm, ax=ax)  # Explicitly associate colorbar with axis
        cbar.set_label(f'Time [{plotconfig["tunit"]}]')
        ax.set_xlabel(f'Position [{plotconfig["lunit"]}]')
        ax.set_ylabel(f'Concentration in Packaging [{plotconfig["Cunit"]}]')
        if title is not None:
            title_main = title
        else:
            title_main = "Concentration Profiles in Packaging vs. Position"
        if subtitle is not None:
            title_sub = subtitle
        else:
            title_sub = rf"$\bf{{{self.name}}}$" + (f": {self.description}" if self.description else "")
        ax.set_title(f"{title_main}\n{title_sub}", fontsize=10)
        ax.text(0.5, 1.05, title_sub, fontsize=8, ha="center", va="bottom", transform=ax.transAxes)
        ax.set_title(title_main)
        ax.grid(True)
        ax.legend()
        if not noshow:
            plt.show()
        # store metadata
        setattr(fig,_fig_metadata_atrr_,f"pltCx_{self.name}")
        return fig

    # export methods (dataframe, csv, xls) for t,CF and tC,x,Cx
    def to_CF_dataframe(self):
        """
        Returns a two-column DataFrame with columns ["time_s", "CF_a.u."],
        containing the raw time vector (self.t) and food concentration (self.CF).
        """
        return pd.DataFrame({
            "time_s": self.t,
            "CF_a.u.": self.CF
        })

    def to_Cx_dataframe(self, t=None, nmax=15, long_format=False):
        """
        Return a DataFrame of concentration profiles (Cx) vs. position (x),
        at user-specified or automatically sampled times.

        Parameters
        ----------
        t : list, array-like, tuple, or None, optional
            - If None, pick up to nmax time points using sqrt-spacing in self.t.
            - If array-like, filter to [t.min(), t.max()] within self.t range.
            - If tuple, e.g. (10, "hours"), convert to SI then use as array.
        nmax : int, optional
            Max # of time points to return if t is None or too large.
        long_format : bool, optional
            - False => wide: rows=x, columns=times
            - True  => tidy: columns=["time_s", "x_m", "C_a.u."]

        Returns
        -------
        pd.DataFrame
            See doc for details.
        """
        # Example from your previously posted logic...
        import numpy as np
        from patankar.layer import check_units

        if self.discrete or (self.Cx is None or self.x is None):
            print("No continuous profile data available.")
            return None

        # 1) pick times
        if t is None:
            nt = len(self.t)
            if nt <= nmax:
                t_values = self.t
            else:
                sqrt_t = np.sqrt(self.t)
                sqrt_t_values = np.linspace(sqrt_t[0], sqrt_t[-1], nmax)
                t_values = sqrt_t_values ** 2
        else:
            if isinstance(t, tuple):
                t_values = check_units(t)[0]
                t_values = np.atleast_1d(t_values)
            else:
                t_values = np.array(t, dtype=float).ravel()
            mask = (t_values >= self.t.min()) & (t_values <= self.t.max())
            t_values = t_values[mask][:nmax]

        if t_values.size == 0:
            print("Warning: no valid times in the specified range.")
            return None

        # 2) gather data
        c_matrix = []
        for tC in t_values:
            c_profile = self.interp_Cx(tC)  # shape => (npositions,)
            c_matrix.append(c_profile)
        c_matrix = np.vstack(c_matrix)  # shape => (#times, #positions)

        # 3) build the DataFrame
        if not long_format:
            # wide => rows=x, columns=time
            df = pd.DataFrame(
                data=c_matrix.T,
                index=self.x,
                columns=t_values
            )
            df.index.name = "x_m"
            df.columns.name = "time_s"
        else:
            # long => columns = time_s, x_m, C_a.u.
            t2d, x2d = np.meshgrid(t_values, self.x, indexing="ij")
            df = pd.DataFrame({
                "time_s": t2d.ravel(),
                "x_m": x2d.ravel(),
                "C_a.u.": c_matrix.ravel()
            })

        return df

    def save_as_excel_CF(self, filename, destinationfolder=None, overwrite=False):
        """
        Save the (t, CF) data to an Excel file.

        Parameters
        ----------
        filename : str
            Output Excel filename (e.g., 'myCF.xlsx').
        destinationfolder : str, optional
            Folder in which to save the file (default: os.getcwd()).
        overwrite : bool, optional
            If False (default), raise an error if file already exists.
        """
        df = self.to_CF_dataframe()
        if df is None:
            print("No CF data to save.")
            return

        if destinationfolder is None:
            destinationfolder = os.getcwd()
        os.makedirs(destinationfolder, exist_ok=True)
        fullpath = os.path.join(destinationfolder, filename)

        if not overwrite and os.path.isfile(fullpath):
            raise FileExistsError(f"File '{fullpath}' already exists. Use overwrite=True to replace it.")

        df.to_excel(fullpath, index=False)
        print(f"CF data saved to Excel: {fullpath}")

    def save_as_csv_CF(self, filename, destinationfolder=None, overwrite=False):
        """
        Save the (t, CF) data to a CSV file.

        Parameters
        ----------
        filename : str
            Output CSV filename (e.g., 'myCF.csv').
        destinationfolder : str, optional
            Folder in which to save the file (default: os.getcwd()).
        overwrite : bool, optional
            If False (default), raise an error if file already exists.
        """
        df = self.to_CF_dataframe()
        if df is None:
            print("No CF data to save.")
            return

        if destinationfolder is None:
            destinationfolder = os.getcwd()
        os.makedirs(destinationfolder, exist_ok=True)
        fullpath = os.path.join(destinationfolder, filename)

        if not overwrite and os.path.isfile(fullpath):
            raise FileExistsError(f"File '{fullpath}' already exists. Use overwrite=True to replace it.")

        df.to_csv(fullpath, index=False)
        print(f"CF data saved to CSV: {fullpath}")

    def save_as_excel_Cx(self, filename, destinationfolder=None, overwrite=False,
                         t=None, nmax=15, long_format=False):
        """
        Save the concentration profiles (Cx) to an Excel file.

        Parameters
        ----------
        filename : str
            Output Excel filename (e.g., 'profiles.xlsx').
        destinationfolder : str, optional
            Folder in which to save the file (default: os.getcwd()).
        overwrite : bool, optional
            If False (default), raise an error if file already exists.
        t : list, array-like, tuple, or None
            Time selection (see to_Cx_dataframe doc).
        nmax : int
            Max # times (see to_Cx_dataframe doc).
        long_format : bool
            True => tidy table; False => wide table.
        """
        df = self.to_Cx_dataframe(t=t, nmax=nmax, long_format=long_format)
        if df is None:
            print("No profile data to save.")
            return

        if destinationfolder is None:
            destinationfolder = os.getcwd()
        os.makedirs(destinationfolder, exist_ok=True)
        fullpath = os.path.join(destinationfolder, filename)

        if not overwrite and os.path.isfile(fullpath):
            raise FileExistsError(f"File '{fullpath}' already exists. Use overwrite=True to replace it.")

        df.to_excel(fullpath, index=not long_format)
        print(f"Profile data (Cx) saved to Excel: {fullpath}")

    def save_as_csv_Cx(self, filename, destinationfolder=None, overwrite=False,
                       t=None, nmax=15, long_format=False):
        """
        Save the concentration profiles (Cx) to a CSV file.

        Parameters
        ----------
        filename : str
            Output CSV filename (e.g., 'profiles.csv').
        destinationfolder : str, optional
            Folder in which to save the file (default: os.getcwd()).
        overwrite : bool, optional
            If False (default), raise an error if file already exists.
        t : list, array-like, tuple, or None
            Time selection (see to_Cx_dataframe doc).
        nmax : int
            Max # times (see to_Cx_dataframe doc).
        long_format : bool
            True => tidy table; False => wide table.
        """
        df = self.to_Cx_dataframe(t=t, nmax=nmax, long_format=long_format)
        if df is None:
            print("No profile data to save.")
            return

        if destinationfolder is None:
            destinationfolder = os.getcwd()
        os.makedirs(destinationfolder, exist_ok=True)
        fullpath = os.path.join(destinationfolder, filename)

        if not overwrite and os.path.isfile(fullpath):
            raise FileExistsError(f"File '{fullpath}' already exists. Use overwrite=True to replace it.")

        df.to_csv(fullpath, index=not long_format)
        print(f"Profile data (Cx) saved to CSV: {fullpath}")



# -----------------------------------------------------
# Container for multiple simulations
# -----------------------------------------------------

class CFSimulationContainer:
    """
    Container to store and compare multiple CF results from different simulations.

    Attributes
    ----------
    curves : dict
        Stores CF results with unique keys. Each entry contains:
        - 'label': Label used for legend.
        - 'tmin', 'tmax': Time range of the simulation.
        - 'interpolant': Interpolated CF function (if continuous).
        - 'times': Discrete time points (if discrete).
        - 'values': Discrete CF values (if discrete).
        - 'color': Assigned color for plotting.
        - 'linestyle': Line style (default is '-').
        - 'linewidth': Line width (default is 2).
        - 'marker': Marker style for discrete data.
        - 'markerfacecolor': Marker face color.
        - 'markersize': Marker size.
        - 'discrete': Boolean indicating discrete data.
    """

    def __init__(self,name="",description="", SML = None, SMLunit=None, plotSML=True, plotconfig=None):
        """
        Initialize an empty container for CF results.
        name: str, default =""
            Name of the container used to save results and shown in title
        description: str, default =""
            description of the container, shown in subtitle
        SML : None or float, optional
            SML value in final units (usually mg/kg, adapt the value accordingly)
        SMLunit : str, optional
            SML units
        plotSML = True (default), bool, optional
            if True and SML is not None, plot SML limit as an horizontal line
        plotconfig : dict, optional
            Dictionary with plotting configuration, containing:
            - "tunit": Time unit label (e.g., 's').
            - "Cunit": Concentration unit label (e.g., 'mg/L').
            - "tscale": Time scaling factor.
            - "Cscale": Concentration scaling factor.
        """
        self.curves = {}
        self._name = name
        self._description = description

        # SML configuration (if any)
        self._SML = SML
        self._SMLunit = SMLunit
        self._plotSML = plotSML

        # Plot configuration
        if plotconfig is not None:
            if not isinstance(plotconfig, dict):
                raise TypeError(f"plotconfig must be a dict not a {type(plotconfig).__name__}")
            self._plotconfig = {**plotconfig_default, **plotconfig}  # Merge without modifying plotconfig_default
        else:
            self._plotconfig = plotconfig_default.copy()  # Work with a copy to prevent accidental changes

    @property
    def name(self):
        return self._name or autoname(6)

    @property
    def description(self):
        return self._description or f"comparison of {len(self.curves)} curves"


    def add(self, simulation_result, label=None, color=None, linestyle="-", linewidth=2,
            marker='o', markerfacecolor='auto', markeredgecolor='black', markersize=6, discrete=False):
        """
        Add a CF result to the container.

        Parameters
        ----------
        simulation_result : SensPatankarResult
            The simulation result.
        discrete : bool, optional
            Whether the data is discrete.
        """
        if not isinstance(simulation_result, SensPatankarResult):
            raise TypeError(f"Expected SensPatankarResult, got {type(simulation_result).__name__}")
        label = label or f"plot{len(self.curves) + 1}"
        key = label[:80]
        if color is None:
            cmap = cm.get_cmap("tab10", len(self.curves) + 1)
            color = cmap(len(self.curves) % 10)
        if markerfacecolor == 'auto':
            markerfacecolor = color
        self.curves[key] = {
            "label": label,
            "color": color,
            "linestyle": linestyle,
            "linewidth": linewidth,
            "marker": marker,
            "markerfacecolor": markerfacecolor,
            "markeredgecolor": markeredgecolor,
            "markersize": markersize,
            "discrete": discrete
        }
        if discrete:
            self.curves[key].update({
                "times": simulation_result.t,
                "values": simulation_result.CF
            })
        else:
            self.curves[key].update({
                "tmin": simulation_result.t.min(),
                "tmax": simulation_result.t.max(),
                "interpolant": simulation_result.interp_CF
            })
        # update SML infos with data added to the container
        if self._SML is None:
            self._SML = simulation_result._SML
        if self._SMLunit is None:
            self._SMLunit = simulation_result._SMLunit

    def delete(self, identifier):
        """
        Remove a stored curve by its index (int) or label (str).

        Parameters
        ----------
        identifier : int or str
            - If `int`, removes the curve at the specified index.
            - If `str`, removes the curve with the matching label.
        """
        if isinstance(identifier, int):
            key = self._get_key_by_index(identifier)
        elif isinstance(identifier, str):
            key = identifier[:40]  # Match the label-based key
            if key not in self.curves:
                print(f"No curve found with label '{identifier}'")
                return
        else:
            raise TypeError("Identifier must be an integer (index) or a string (label).")

    def __repr__(self):
        """Return a summary of stored CF curves including index numbers."""
        if not self.curves:
            return "<CFSimulationContainer: No stored curves>"
        repr_str = "<CFSimulationContainer: Stored CF Curves>\n"
        repr_str += "--------------------------------------------------\n"
        for index, (key, data) in enumerate(self.curves.items()):
            repr_str += (f"[{index}] Label: {data['label']} | "
                         f"Time: [{data['tmin']:.2e}, {data['tmax']:.2e}] s | "
                         f"Color: {data['color']} | "
                         f"Style: {data['linestyle']} | "
                         f"Width: {data['linewidth']}\n")
        return repr_str

    def _validate_indices(self, indices):
        """Helper function to check if indices are valid."""
        if isinstance(indices, int):
            indices = [indices]
        if not all(isinstance(i, int) and 0 <= i < len(self.curves) for i in indices):
            raise IndexError(f"Invalid index in {indices}. Must be between 0 and {len(self.curves) - 1}.")
        return indices

    def _get_keys_by_indices(self, indices):
        """Helper function to retrieve keys based on indices."""
        if isinstance(indices, (int, str)):
            indices = [indices]
        keys = []
        all_keys = list(self.curves.keys())
        for idx in indices:
            if isinstance(idx, int):
                if idx < 0:
                    idx += len(all_keys)
                if idx < 0 or idx >= len(all_keys):
                    raise IndexError(f"Index {idx} is out of range for curves.")
                keys.append(all_keys[idx])
            elif isinstance(idx, str):
                if idx not in self.curves:
                    raise KeyError(f"Key '{idx}' does not exist in curves.")
                keys.append(idx)
            else:
                raise TypeError("Index must be an int, str, or a list of both.")
        return keys

    def update(self, index, label=None, linestyle=None, linewidth=None, color=None,
               marker=None, markersize=None, markerfacecolor=None, markeredgecolor=None):
        """
        Update properties of one or multiple curves.

        Parameters
        ----------
        index : int or list of int
            Index or indices of the curve(s) to update.
        label : str, optional
            New label for the curve(s).
        linestyle : str, optional
            New linestyle for the curve(s).
        linewidth : float, optional
            New linewidth for the curve(s).
        color : str or tuple, optional
            New color for the curve(s).
        marker : str, optional
            New marker style for discrete data.
        markersize : float, optional
            New marker size for discrete data.
        markerfacecolor : str or tuple, optional
            New marker face color.
        markeredgecolor : str or tuple, optional
            New marker edge color.
        """
        keys = self._get_keys_by_indices(index)

        for key in keys:
            if label is not None:
                self.curves[key]["label"] = label
            if linestyle is not None:
                self.curves[key]["linestyle"] = linestyle
            if linewidth is not None:
                self.curves[key]["linewidth"] = linewidth
            if color is not None:
                self.curves[key]["color"] = color
            if marker is not None:
                self.curves[key]["marker"] = marker
            if markersize is not None:
                self.curves[key]["markersize"] = markersize
            if markerfacecolor is not None:
                self.curves[key]["markerfacecolor"] = markerfacecolor
            if markeredgecolor is not None:
                self.curves[key]["markeredgecolor"] = markeredgecolor


    def label(self, index, new_label):
        """Change the label of one or multiple curves."""
        self.update(index, label=new_label)

    def linewidth(self, index, new_value):
        """Change the linewidth of one or multiple curves."""
        self.update(index, linewidth=new_value)

    def linestyle(self, index, new_style):
        """Change the linestyle of one or multiple curves."""
        self.update(index, linestyle=new_style)

    def color(self, index, new_color):
        """Change the color of one or multiple curves."""
        self.update(index, color=new_color)

    def marker(self, index, new_marker):
        """Change the marker style of one or multiple curves."""
        self.update(index, marker=new_marker)

    def markersize(self, index, new_size):
        """Change the marker size of one or multiple curves."""
        self.update(index, markersize=new_size)

    def markerfacecolor(self, index, new_facecolor):
        """Change the marker face color of one or multiple curves."""
        self.update(index, markerfacecolor=new_facecolor)

    def markeredgecolor(self, index, new_edgecolor):
        """Change the marker edge color of one or multiple curves."""
        self.update(index, markeredgecolor=new_edgecolor)

    def colormap(self, name="viridis", ncolors=16, tooclearflag=True, reverse=False):
        """
        Generates a list of `ncolors` colors from the specified colormap.

        Parameters:
        -----------
        name : str, optional (default="viridis")
            Name of the Matplotlib colormap to use.
        ncolors : int, optional (default=16)
            Number of colors to generate.
        tooclearflag : bool, optional (default=True)
            If True, applies `tooclear` function to adjust brightness.
        reverse : bool, optional (default=False)
            If True, reverses the colormap.

        Returns:
        --------
        list of tuples
            List of RGB(A) colors in [0,1] range.

        Raises:
        -------
        ValueError
            If the colormap name is not recognized.
        """
        return colormap(name, ncolors, tooclear, reverse)

    def viridis(self, ncolors=16, tooclear=True, reverse=False):
        """Generates colors from the Viridis colormap."""
        return colormap("viridis", ncolors, tooclear, reverse)

    def jet(self, ncolors=16, tooclear=True, reverse=False):
        """Generates colors from the Jet colormap."""
        return colormap("jet", ncolors, tooclear, reverse)


    def plotCF(self, t_range=None, SML = None, SMLunit=None, plotSML=None, plotconfig=None):
        """
        Plot all stored CF curves in a single figure.

        Parameters
        ----------
        t_range : tuple (t_min, t_max), optional
            Time range for plotting. If None, uses each curve's own range.
        SML : float, optional
            user SML
        SMLunit: str, optional
            user SML units
        plotSML: bool, default = True
            True plots SML line
        plotconfig : dict, optional
            Dictionary with plotting configuration, containing:
            - "tunit": Time unit label (e.g., 's').
            - "Cunit": Concentration unit label (e.g., 'mg/L').
            - "tscale": Time scaling factor.
            - "Cscale": Concentration scaling factor.
        """
        # Plot config
        # force LaTeX only on systems with latex installed
        plt.rc('text', usetex=_LaTeXavailable) # Enable LaTeX formatting for Matplotlib
        # Create a temporary plotconfig without modifying self._plotconfig
        if plotconfig is not None:
            if not isinstance(plotconfig, dict):
                raise TypeError(f"plotconfig must be a dict not a {type(plotconfig).__name__}")
            plotconfig = {**self._plotconfig, **plotconfig}  # Merge without modifying self._plotconfig
        else:
            plotconfig = self._plotconfig.copy()  # Work with a copy to prevent accidental changes

        # retrieve the SML plot choice (set at construction with possible user overrides)
        SML = self._SML if SML is None else SML
        SMLunit = self._SMLunit if SMLunit is None else SMLunit
        plotSML = self._plotSML if plotSML is None else plotSML

        # we apply eventual user overrides
        plotconfig = useroverride("plotconfig",plotconfig,expected_type=dict)
        plotSML = useroverride("plotSML",plotSML,expected_type=bool)

        if not self.curves:
            print("No curves to plot.")
            return

        fig, ax = plt.subplots(figsize=(8, 6))

        for data in self.curves.values():
            if data["discrete"]:
                # Discrete data plotting
                ax.scatter(data["times"], data["values"], label=data["label"],
                           color=data["color"], marker=data["marker"],
                           facecolor=data["markerfacecolor"], edgecolor=data["markeredgecolor"],
                           s=data["markersize"]**2)
            else:
                # Continuous data plotting
                t_min, t_max = data["tmin"], data["tmax"]
                if t_range:
                    t_min, t_max = max(t_min, t_range[0]), min(t_max, t_range[1])

                t_plot = np.linspace(t_min, t_max, 500)
                CF_plot = data["interpolant"](t_plot)
                ax.plot(t_plot, CF_plot, label=cleantex(data["label"]),
                        color=data["color"], linestyle=data["linestyle"], linewidth=data["linewidth"])

        # add SML values (we assume that the units of SML are already final, no need to use plotconfig["Cscale"])
        if plotSML:
            SML = self._SML if SML is None else SML
            SMLunit = self._SMLunit if SMLunit is None else SMLunit
            if SML is not None:
                tmiddle = (t_min+t_max)/2
                SMLunit = plotconfig["Cunit"] if SMLunit is None else SMLunit
                ax.axhline(SML, color="ForestGreen", linestyle=(0, (3, 5, 1, 5)),
                           linewidth=2, label=f"SML={SML:.4g} {SMLunit}")
                ax.text(tmiddle / plotconfig["tscale"], SML, "specific migration limit", fontsize=8, color='ForestGreen', ha='center', va='bottom')

        # Configure the plot
        ax.set_xlabel(f'Time [{plotconfig["tunit"]}]' if plotconfig else "Time")
        ax.set_ylabel(f'Concentration in Food [{plotconfig["Cunit"]}]' if plotconfig else "CF")
        title_main = "Concentration in Food vs. Time"
        title_sub = cleantex(rf"$\bf{{{self.name}}}$") + (f": {self.description}" if self.description else "")
        ax.set_title(f"{title_main}\n{title_sub}", fontsize=10)
        ax.text(0.5, 1.05, title_sub, fontsize=8, ha="center", va="bottom", transform=ax.transAxes)
        ax.set_title(title_main)
        ax.legend()
        ax.grid(True)
        plt.show()
        # store metadata
        setattr(fig,_fig_metadata_atrr_,f"cmp_pltCF_{self.name}")
        return fig


    def to_dataframe(self, t_range=None, num_points=1000, time_list=None):
        """
        Export interpolated CF data as a pandas DataFrame.
        Parameters:
        - t_range: tuple (t_min, t_max), optional
            The time range for interpolation (default: min & max of all stored results).
        - num_points: int, optional
            Number of points in the interpolated time grid (default: 100).
        - time_list: list or array, optional
            Explicit list of time points for interpolation (overrides t_range & num_points).
        Returns:
        - pd.DataFrame
            A DataFrame with time as index and CF values as columns (one per simulation).
        """
        if not self.curves:
            print("No data to export.")
            return pd.DataFrame()

        # Determine the time grid
        if time_list is not None:
            t_grid = np.array(time_list)
        else:
            all_t_min = min(data["tmin"] for data in self.curves.values())
            all_t_max = max(data["tmax"] for data in self.curves.values())
            # Default time range
            t_min, t_max = t_range if t_range else (all_t_min, all_t_max)
            # Create evenly spaced time grid
            t_grid = np.linspace(t_min, t_max, num_points)
        # Create DataFrame with time as index
        df = pd.DataFrame({"Time (s)": t_grid})

        # Interpolate each stored CF curve at the common time grid
        for key, data in self.curves.items():
            df[data["label"]] = data["interpolant"](t_grid)
        return df


    def save_as_excel(self, filename="CF_data.xlsx", destinationfolder=os.getcwd(), overwrite=False,
                      t_range=None, num_points=1000, time_list=None):
        """
        Save stored CF data to an Excel file.
        Parameters:
        - filename: str, Excel filename.
        - destinationfolder: str, where to save the file.
        - overwrite: bool, overwrite existing file.
        - t_range: tuple (t_min, t_max), optional
            The time range for interpolation (default: min & max of all stored results).
        - num_points: int, optional
            Number of points in the interpolated time grid (default: 100).
        - time_list: list or array, optional
            Explicit list of time points for interpolation (overrides t_range & num_points).
        """
        if not self.curves:
            print("No data to export.")
            return
        df = self.to_dataframe(t_range=t_range, num_points=num_points, time_list=time_list)
        filepath = os.path.join(destinationfolder, filename)
        if not overwrite and os.path.exists(filepath):
            print(f"File {filepath} already exists. Use overwrite=True to replace it.")
            return

        df.to_excel(filepath, index=False)
        print(f"Saved Excel file: {filepath}")


    def save_as_csv(self, filename="CF_data.csv", destinationfolder=os.getcwd(), overwrite=False,
                    t_range=None, num_points=200, time_list=None):
        """
        Save stored CF data to an Excel file.
        Parameters:
        - filename: str, Excel filename.
        - destinationfolder: str, where to save the file.
        - overwrite: bool, overwrite existing file.
        - t_range: tuple (t_min, t_max), optional
            The time range for interpolation (default: min & max of all stored results).
        - num_points: int, optional
            Number of points in the interpolated time grid (default: 100).
        - time_list: list or array, optional
            Explicit list of time points for interpolation (overrides t_range & num_points).
        """
        if not self.curves:
            print("No data to export.")
            return
        df = self.to_dataframe(t_range=t_range, num_points=num_points, time_list=time_list)
        filepath = os.path.join(destinationfolder, filename)
        if not overwrite and os.path.exists(filepath):
            print(f"File {filepath} already exists. Use overwrite=True to replace it.")
            return
        df.to_csv(filepath, index=False)
        print(f"Saved CSV file: {filepath}")


    def rgb(self):
        """Displays a categorized color chart with properly aligned headers."""
        plt.rc('text', usetex=False) # Enable LaTeX formatting for Matplotlib
        rgb()


# restartfile
class restartfile:
    """
    A container class for storing simulation restart data.

    This class facilitates storing and restoring simulation parameters and results,
    allowing simulations to be resumed or analyzed after computation.

    Methods:
    --------
    copy(what)
        Creates a deep copy of various data types to ensure safety in storage.

    Example:
    --------
    ```python
    restart = restartfile()
    copy_data = restart.copy([1, 2, 3])
    ```
    """
    @classmethod
    def copy(cls, what):
        """Safely copy a parameter that can be a float, str, dict, or a NumPy array"""
        if isinstance(what, (int, float, str, tuple,bool)):  # Immutable types (direct copy)
            return what
        elif isinstance(what, np.ndarray):  # NumPy array (ensure a separate copy)
            return np.copy(what)
        elif isinstance(what, dict):  # Dictionary (deep copy)
            return duplicate(what)
        elif what is None:
            return None
        else:  # Fallback for other complex types
            return duplicate(what)

# specific restartfile for senspatankar
class restartfile_senspantakar(restartfile):
    """
    Specialized restart file container for the `senspatankar` migration solver.

    This class stores the simulation inputs and computed results, enabling
    the resumption of a simulation from a saved state.

    Attributes:
    -----------
    inputs : dict
        Stores all initial simulation inputs.
    t : float or None
        Simulation time at the stored state.
    CF : float or None
        Concentration in food at the stored state.
    Cprofile : Cprofile or None
        Concentration profile at the stored state.

    Methods:
    --------
    freezeCF(t, CF)
        Saves the food concentration `CF` at time `t`.
    freezeCx(x, Cx)
        Saves the concentration profile `Cx` over `x`.

    Example:
    --------
    ```python
    restart = restartfile_senspatankar(multilayer, medium, name, description, ...)
    restart.freezeCF(t=1000, CF=0.05)
    ```
    """
    def __init__(self,multilayer,medium,name,description,
                 t,autotime,timescale,Cxprevious,
                 ntimes,RelTol,AbsTol,deepcopy=True):
        """constructor to be called at the intialization"""

        # eventual user override (for memory management)
        # deepcopy=False (unsecure copy) must never be overriden (if not the kernel will crash !!!)
        deepcopy = useroverride("deepcopy",deepcopy,expected_type=bool) if deepcopy else deepcopy

        if deepcopy:
            inputs = {
                "multilayer":multilayer.copy(),
                "medium":medium.copy(),
                "name":restartfile.copy(name),
                "description":restartfile.copy(description),
                "t":restartfile.copy(t), # t is a duration not absolute time (it should not be reused)
                "autotime":restartfile.copy(autotime),
                "timescale":restartfile.copy(timescale),
                "Cxprevious":Cxprevious,
                "ntimes":restartfile.copy(ntimes),
                "RelTol":restartfile.copy(RelTol),
                "AbsTol":restartfile.copy(AbsTol)
                }
        else:
            inputs = {
                "multilayer":multilayer,
                "medium":medium,
                "name":name,
                "description":description,
                "t":t, # t is a duration not absolute time (it should not be reused)
                "autotime":autotime,
                "timescale":timescale,
                "Cxprevious":Cxprevious,
                "ntimes":ntimes,
                "RelTol":RelTol,
                "AbsTol":AbsTol
                }
        # inputs
        self.inputs = inputs
        # outputs
        self.t = None # no result yet
        self.CF = None # no result yet
        self.Cprofile = None # no result yet

    def freezeCF(self,t,CF):
        """Freeze the CF solution CF(t)"""
        self.t = t
        self.CF = CF

    def freezeCx(self,x,Cx):
        """Freeze the Cx solution Cx(x)"""
        self.Cprofile = Cprofile(x,Cx)

    def __repr__(self):
        """representation of the restart object"""
        if self.t is None:
            print("Restart file with no result")
        else:
            print(f"Restart file at t={self.t} with CF={self.CF}")
            print("Details of the profile:")
            repr(self.Cprofile)
        return str(self)

    def __str__(self):
        """Formatted representation of the restart object"""
        res = "no result" if self.t is None else f"solution at t={self.t}"
        return f"<{self.__class__.__name__}: {res}"


# %% Core function
def senspatankar(multilayer=None, medium=None,
                 name=f"senspatantkar:{autoname(6)}", description="",
                 t=None, autotime=True, timescale="sqrt", Cxprevious=None,
                 ntimes=1000, RelTol=1e-6, AbsTol=1e-6,
                 container=None):
    """
    Simulates in 1D the mass transfer of a substance initially distributed in a multilayer
    packaging structure into a food medium (or liquid medium). This solver uses a finite-volume
    method adapted from Patankar to handle partition coefficients between all layers, and
    between the food and the contact layer.

    Two typical configurations are implemented:

    Configuration (PBC=False)
        - Robin (third-kind boundary condition) on the left (in contact with food)
        - Impervious boundary condition on the right (in contact with surrounding)

    Configuration (PBC=true)
        - periodic boundary condition between left and right to simulate infinite stacking or setoff

    The configuration nofood is a variant of PBC=False with h=Bi=0 (impervious boundary condition on the left).

    The behavior of the solver is decided by medium attributes (see food.py module).
    The property medium.PBC will determine whether periodic boundary conditions are used or not.


    Parameters
    ----------
    multilayer : layer
        A ``layer`` (or combined layers) object describing the packaging.
    medium : foodlayer or foodphysics
        A ``foodlayer`` object describing the food (or liquid) medium in contact.
    name : str, optional
        Simulation name, default = f"senspatantkar:{autoname(6)}" where autoname(6)
        is a random sequence of characters a-z A-Z 0-9
    description : str, optional
        Simulation description
    t : float or array_like, optional
        If a float is provided, it is taken as the total contact duration in seconds.
        If an array is provided, it is assumed to be time points where the solution
        will be evaluated. If None, it defaults to the contact time from the medium.
    autotime : bool, optional
        If True (default), an automatic time discretization is generated internally
        (linear or sqrt-based) between 0 and tmax (the maximum time). If False, the
        times in ``t`` are used directly.
    timescale : {"sqrt", "linear"}, optional
        Type of automatic time discretization if ``autotime=True``.
        "sqrt" (default) refines the early times more (useful for capturing rapid changes).
        "linear" uses a regular spacing.
    Cxprevious : Cprofile, optional (default=None)
        Concentration profile (from a previous simulation).
    ntimes : int, optional
        Number of time points in the automatically generated time vector if ``autotime=True``.
        The default is 1e3.
    RelTol : float, optional
        Relative tolerance for the ODE solver (``solve_ivp``). Default is 1e-4.
    AbsTol : float, optional
        Absolute tolerance for the ODE solver (``solve_ivp``). Default is 1e-4.

    Raises
    ------
    TypeError
        If ``multilayer`` is not a ``layer`` instance or ``medium`` is not a ``foodlayer`` instance,
        or if ``timescale`` is not a string.
    ValueError
        If an invalid ``timescale`` is given (not one of {"sqrt", "linear"}).

    Returns
    -------
    SensPatankarResult
        An object containing the time vector, concentration histories, fluxes, and
        spatial concentration profiles suitable for plotting and analysis.

    Notes
    -----
    - The geometry is assumed 1D: Food is on the left boundary, with a mass transfer coefficient
      `h = medium.h`, partition ratio `k0 = medium.k0`, and the packaging layers are to the right
      up to an impervious boundary.
    - Results are normalized internally using a reference layer (``iref``) specified in ``multilayer``.
      The reference layer is used to define dimensionless time (Fourier number Fo).
    - The dimensionless solution is solved by the Patankar approach with partition coefficients.

    Example
    -------
    .. code-block:: python

        from patankar.food import ethanol
        from patankar.layer import layer
        medium = ethanol()
        A = layer(layername="layer A")
        B = layer(layername="layer B")
        multilayer = A + B

        sol = senspatankar(multilayer, medium, autotime=True, timescale="sqrt")
        sol.plotCF()
        sol.plotC()
    """
    # Check arguments
    if not isinstance(multilayer, layer):
        raise TypeError(f"the input multilayer must be of class layer, not {type(multilayer).__name__}")
    if not isinstance(medium, (foodlayer,foodphysics)):
        raise TypeError(f"the input medium must be of class foodlayer, not {type(medium).__name__}")
    if not isinstance(timescale, str):
        raise TypeError(f"timescale must be a string, not {type(timescale).__name__}")

    # Apply User overrides if any
    timescale = useroverride("timescale",timescale,valuelist=("linear","sqrt"))
    ntimes = useroverride("ntimes",ntimes,valuemin=10,valuemax=20000)
    RelTol = useroverride("RelTol",RelTol,valuemin=1e-9,valuemax=1e-3)
    AbsTol = useroverride("AbsTol",AbsTol,valuemin=1e-9,valuemax=1e-3)

    # Refresh the physics of medium for parameters tunned by the end-user
    medium.refresh()

    # extract the PBC flag (True for setoff)
    PBC = medium.PBC

    # Restart file initialization (all parameters are saved - and cannot be changed)
    restart = restartfile_senspantakar(multilayer, medium, name,
            description, t, autotime, timescale, Cxprevious, ntimes, RelTol, AbsTol,deepcopy=True)
    # Restart file (unsecure version without deepcoy)
    restart_unsecure = restartfile_senspantakar(multilayer, medium, name,
            description, t, autotime, timescale, Cxprevious, ntimes, RelTol, AbsTol,deepcopy=False)

    # Contact medium properties
    CF0 = medium.get_param("CF0",0) # instead of medium.CF0 to get a fallback mechanism with nofood and setoff
    k0 = medium.get_param("k0",1)
    h = medium.get_param("h",0,acceptNone=False) # None will arise for PBC
    ttarget = medium.get_param("contacttime") # <-- ttarget is the time requested
    tmax = 2 * ttarget  # ensures at least up to 2*contacttime

    # Material properties
    k = multilayer.k / k0   # all k are normalized
    k0 = k0 / k0            # all k are normalized
    D = multilayer.D
    l = multilayer.l
    C0 = multilayer.C0

    # Validate/prepare time array
    if isinstance(t,tuple):
        t = check_units(t)[0]
    t = np.array(tmax if t is None else t, dtype=float) # <-- simulation time (longer than ttarget)
    if np.isscalar(t) or t.size == 1:
        t = np.array([0, t.item()],dtype=float)
    if t[0] != 0:
        t = np.insert(t, 0, 0)  # Ensure time starts at zero
    # Ensure t[-1] is greater than ttarget
    if t[-1] < ttarget.item():  # Convert ttarget to scalar before comparison
        t = np.append(t, [ttarget, 1.05*ttarget, 1.1*ttarget, 1.2*ttarget])  # Extend time array to cover requested time

    # Reference layer for dimensionless transformations
    iref = multilayer.referencelayer
    l_ref = l[iref]
    D_ref = D[iref]

    # Normalize lengths and diffusivities
    l_normalized = l / l_ref
    D_normalized = D / D_ref

    # Dimensionless time (Fourier number)
    timebase = l_ref**2 / D_ref
    Fo = t / timebase

    # Automatic time discretization if requested
    if autotime:
        if timescale.lower() == "linear":
            Fo_int = np.linspace(np.min(Fo), np.max(Fo), int(ntimes))
        elif timescale.lower() == "sqrt":
            Fo_int = np.linspace(np.sqrt(np.min(Fo)), np.sqrt(np.max(Fo)), int(ntimes))**2
        else:
            raise ValueError('timescale can be "sqrt" or "linear"')
        t = Fo_int * timebase
    else:
        Fo_int = Fo

    # L: dimensionless ratio of packaging to food volumes (scaled by reference layer thickness)
    A = medium.get_param("surfacearea",0)
    l_sum = multilayer.thickness
    VP = A * l_sum
    VF = medium.get_param("volume",1)
    LPF = VP / VF
    L = LPF * l_ref / l_sum

    # Bi: dimensionless mass transfer coefficient
    Bi = h * l_ref / D_ref

    # Compute equilibrium concentration factor
    sum_lL_C0 = np.sum(l_normalized * L * C0)
    sum_terms = np.sum((1 / k) * l_normalized * L)
    C0eq = (CF0 + sum_lL_C0) / (1 + sum_terms)
    if C0eq == 0:
        C0eq = 1.0

    # Normalize initial concentrations
    C0_normalized = C0 / C0eq
    CF0_normalized = CF0 / C0eq

    # Generate mesh (add offset x0 and concatenate them)
    meshes = multilayer.mesh()
    x0 = 0
    for i,mesh in enumerate((meshes)):
        mesh.xmesh += x0
        x0 += mesh.l
    xmesh = np.concatenate([m.xmesh for m in meshes])
    total_nodes = len(xmesh)

    # Positions of the interfaces (East and West)
    dw = np.concatenate([m.dw for m in meshes])
    de = np.concatenate([m.de for m in meshes])

    # Attach properties to nodes (flat interpolant)
    D_mesh = np.concatenate([D_normalized[m.index] for m in meshes])
    k_mesh = np.concatenate([k[m.index] for m in meshes])
    C0_mesh = np.concatenate([C0_normalized[m.index] for m in meshes])

    # Interpolate the initial solution if Cxprevious is supplied
    if Cxprevious is not None:
        if not isinstance(Cxprevious,Cprofile):
            raise TypeError(f"Cxprevisous should be a Cprofile object not a {type(Cxprevious).__name__}")
        C0_mesh = Cxprevious.interp(xmesh*l_ref) / C0eq # dimensionless

    # Conductances between the node and the next interface
    # item() is forced to avoid the (1,) Shape Issue (since NumPy 1.25)
    hw = np.zeros(total_nodes)
    he = np.zeros(total_nodes)
    if PBC:
        for i in range(total_nodes):
            prev = total_nodes-1 if i==0 else i-1
            hw[i] = (1 / ((de[prev] / D_mesh[prev] * k_mesh[prev] / k_mesh[i]) + dw[i] / D_mesh[i])).item()
    else:
        if Bi.item()==0:
            hw[0] = 0 # to prevent RuntimeWarning: divide by zero encountered in divide
        else:
            hw[0] = (1 / ((1 / k_mesh[0]) / Bi + dw[0] / D_mesh[0])).item()
        for i in range(1, total_nodes):
            hw[i] = (1 / ((de[i - 1] / D_mesh[i - 1] * k_mesh[i - 1] / k_mesh[i]) + dw[i] / D_mesh[i])).item()
    he[:-1] = hw[1:] # nodes are the center of FV elements: he = np.roll(hw, -1)
    he[-1]=hw[0] if PBC else 0.0 # we connect (PBC) or we enforce impervious (note that he was initialized to 0 already)

    if PBC: # periodic boundary condition

        # Assemble sparse matrix using COO format for efficient construction
        rows = np.zeros(3 * total_nodes, dtype=int) # row indices
        cols = np.zeros_like(rows) # col indices
        data = np.zeros_like(rows, dtype=np.float64) # values
        idx = 0
        for i in range(total_nodes):
            current = i
            west = (i-1) % total_nodes
            east = (i+1) % total_nodes
            denominator = dw[current] + de[current]
            k_current = k_mesh[current]
            k_west = k_mesh[west]
            k_east = k_mesh[east]
            # West neighbor
            rows[idx] = current
            cols[idx] = west
            data[idx] = hw[current] * k_west / k_current / denominator
            idx +=1
            # Diagonal
            rows[idx] = current
            cols[idx] = current
            data[idx] = (-hw[current] - he[current] * k_current/k_east) / denominator
            idx +=1
            # East neighbor
            rows[idx] = current
            cols[idx] = east
            data[idx] = he[current] / denominator
            idx +=1
        A = coo_matrix((data[:idx], (rows[:idx], cols[:idx])),
                     shape=(total_nodes, total_nodes)).tocsr()
        C_initial =  C0_mesh

    else: # Robin (left) + impervious (right) --> triband matrix

        # Assemble the tri-band matrix A as sparse for efficiency
        size = total_nodes + 1  # +1 for the food node
        main_diag = np.zeros(size)
        upper_diag = np.zeros(size - 1)
        lower_diag = np.zeros(size - 1)
        # Food node (index 0)
        main_diag[0] = (-L * hw[0] * (1 / k_mesh[0])).item()
        upper_diag[0] = (L * hw[0]).item()
        # Layer nodes
        for i in range(total_nodes):
            denom = dw[i] + de[i]
            if i == 0:
                main_diag[1] = (-hw[0] - he[0] * k_mesh[0] / k_mesh[1]) / denom
                upper_diag[1] = he[0] / denom
                lower_diag[0] = (hw[0] * (1 / k_mesh[0])) / denom
            elif i == total_nodes - 1:
                main_diag[i + 1] = (-hw[i]) / denom
                lower_diag[i] = (hw[i] * k_mesh[i - 1] / k_mesh[i]) / denom
            else:
                main_diag[i + 1] = (-hw[i] - he[i] * k_mesh[i] / k_mesh[i + 1]) / denom
                upper_diag[i + 1] = he[i] / denom
                lower_diag[i] = (hw[i] * k_mesh[i - 1] / k_mesh[i]) / denom
        A = diags([main_diag, upper_diag, lower_diag], [0, 1, -1], shape=(size, size), format='csr')
        C_initial = np.concatenate([CF0_normalized, C0_mesh])

    # ODE system: dC/dFo = A * C
    def odesys(_, C):
        return A.dot(C)

    sol = solve_ivp(   # <-- generic solver
        odesys,        # <-- our system (efficient sparse matrices)
        [Fo_int[0], Fo_int[-1]], # <-- integration range on Fourier scale
        C_initial,     # <-- initial solution
        t_eval=Fo_int, # <-- the solution is retrieved at these Fo values
        method='BDF',  # <-- backward differences are absolutely stable
        rtol=RelTol,   # <-- relative and absolute tolerances
        atol=AbsTol
    )

    # Check solution
    if not sol.success:
        print("Solver failed:", sol.message)

    # Extract solution
    if PBC:
        CF_dimless = np.full((sol.y.shape[1],), CF0 / C0eq)
        C_dimless = sol.y
        f = np.zeros_like(CF_dimless)
    else:
        CF_dimless = sol.y[0, :]
        C_dimless = sol.y[1:, :]
        # Robin flux
        f = hw[0] * (k0 * CF_dimless - C_dimless[0, :]) * C0eq

    # Compute cumulative flux
    fc = cumulative_trapezoid(f, t, initial=0)

    if PBC:
        # Build full (dimensionless) profile for plotting across each sub-node
        xfull, Cfull_dimless = compute_fc_profile_PBC(C_dimless, Fo_int, de, dw, he, hw, k_mesh, D_mesh, xmesh, xreltol=0)
        # Build full (dimensionless) profile for interpolation across each sub-node
        xfulli, Cfull_dimlessi = compute_fc_profile_PBC(C_dimless, Fo_int, de, dw, he, hw, k_mesh, D_mesh, xmesh, xreltol=1e-4)
    else:
        # Build full (dimensionless) profile for plotting across each sub-node
        xfull, Cfull_dimless = compute_fv_profile(xmesh, dw, de,C_dimless, k_mesh, D_mesh, hw, he, CF_dimless, k0, Fo_int, xreltol=0)
        # Build full (dimensionless) profile for interpolation across each sub-node
        xfulli, Cfull_dimlessi = compute_fv_profile(xmesh, dw, de,C_dimless, k_mesh, D_mesh, hw, he, CF_dimless, k0, Fo_int, xreltol=1e-4)


    # revert to dimensional concentrations
    CF = CF_dimless * C0eq
    Cx = Cfull_dimless * C0eq

    return SensPatankarResult(
        name=name,
        description=description,
        ttarget = ttarget,             # target time
        t=t,     # time where concentrations are calculated
        C= np.trapz(Cfull_dimless, xfull, axis=1)*C0eq,
        CF=CF,
        fc=fc,
        f=f,
        x=xfull * l_ref,           # revert to dimensional lengths
        Cx=Cx,
        tC=sol.t,
        C0eq=C0eq,
        timebase=timebase,
        restart=restart, # <--- restart info (inputs only)
        restart_unsecure=restart_unsecure,
        xi=xfulli*l_ref, # for restart only
        Cxi=Cfull_dimlessi*C0eq, # for restart only
        createcontainer = True,
        container=container
    )


# Exact FV interpolant (with Robin BC)
def compute_fv_profile(xmesh, dw, de, C_dimless, k_mesh, D_mesh, hw, he, CF_dimless, k0, Fo_int, xreltol=0):
    """
    Compute the full finite-volume concentration profile, including node values and interface values.
    (this function is not nested inside senspantar for better readability)

    Parameters:
        xmesh (np.ndarray): Node positions.
        dw (np.ndarray): Distance to west interfaces.
        de (np.ndarray): Distance to east interfaces.
        C_dimless (np.ndarray): Concentration at nodes.
        k_mesh (np.ndarray): Henri-like coefficient at nodes.
        D_mesh (np.ndarray): Diffusion coefficient at nodes.
        hw (np.ndarray): Conductance to west interface.
        he (np.ndarray): Conductance to east interface.
        CF_dimless (np.ndarray): Far-field (Food) concentration values.
        k0 (float): Partition coefficient at the boundary.
        Fo_int (np.ndarray): Time steps.
        xreltol (float, optional): Relative perturbation factor for interpolation accuracy. Defaults to 0.

    Returns:
        xfull (np.ndarray): Full spatial positions including nodes and interfaces.
        Cfull_dimless (np.ndarray): Full concentration profile.
    """
    num_nodes, num_timesteps = C_dimless.shape  # Extract shape

    # Compute xtol based on minimum interface distances
    xtol = np.min([np.min(de), np.min(dw)]) * xreltol

    # Adjust west and east interface positions
    xw = xmesh - dw + xtol  # Shift west interface
    xe = xmesh + de - xtol  # Shift east interface

    # Build full spatial profile
    xfull = np.empty(3 * num_nodes,dtype=np.float64)
    xfull[::3] = xw      # Every 3rd position is xw
    xfull[1::3] = xmesh  # Every 3rd position (offset by 1) is xmesh
    xfull[2::3] = xe     # Every 3rd position (offset by 2) is xe

    # Initialize concentration at interfaces
    Ce = np.zeros_like(C_dimless)  # East interfaces
    Cw = np.zeros_like(C_dimless)  # West interfaces

    # Compute Ce (east interface) for all timesteps at once
    Ce[:-1, :] = C_dimless[:-1, :] - (
        (de[:-1, None] * he[:-1, None] *
        ((k_mesh[:-1, None] / k_mesh[1:, None]) * C_dimless[:-1, :] - C_dimless[1:, :]))
        / D_mesh[:-1, None]
    )
    Ce[-1, :] = C_dimless[-1, :]  # Last node follows boundary condition

    # Compute Cw (west interface) for all timesteps at once
    Cw[1:, :] = C_dimless[1:, :] + (
        (dw[1:, None] * hw[1:, None] *
        ((k_mesh[:-1, None] / k_mesh[1:, None]) * C_dimless[:-1, :] - C_dimless[1:, :]))
        / D_mesh[1:, None]
    )

    # Compute Cw[0, :] separately to handle boundary condition
    Cw[0, :] = (C_dimless[0, :] + (
        dw[0] * hw[0] *
        (k0 / k_mesh[0] * CF_dimless - C_dimless[0, :])
        / D_mesh[0]
    )).flatten()  # Ensure correct shape

    # Interleave concentration values instead of using np.hstack and reshape
    Cfull_dimless = np.empty((num_timesteps, 3 * num_nodes),dtype=np.float64)
    Cfull_dimless[:, ::3] = Cw.T      # Every 3rd column is Cw
    Cfull_dimless[:, 1::3] = C_dimless.T  # Every 3rd column (offset by 1) is C
    Cfull_dimless[:, 2::3] = Ce.T      # Every 3rd column (offset by 2) is Ce

    return xfull, Cfull_dimless


def compute_fc_profile_PBC(C, t, de, dw, he, hw, k, D, xmesh, xreltol=0):
    """
    Computes the full concentration profile, including interface concentrations,
    for a system with periodic boundary conditions (PBC).

    This function calculates the concentrations at the east (`Ce`) and west (`Cw`)
    interfaces of each finite volume node, ensuring periodicity in the domain.

    Parameters
    ----------
    C : np.ndarray, shape (num_nodes, num_timesteps)
        Concentration values at each node for all time steps.
    t : np.ndarray, shape (num_timesteps,)
        Time points at which concentration profiles are computed.
    de : np.ndarray, shape (num_nodes,)
        Eastward diffusion lengths at each node.
    dw : np.ndarray, shape (num_nodes,)
        Westward diffusion lengths at each node.
    he : np.ndarray, shape (num_nodes,)
        Eastward mass transfer coefficients.
    hw : np.ndarray, shape (num_nodes,)
        Westward mass transfer coefficients.
    k : np.ndarray, shape (num_nodes,)
        Partition coefficients at each node.
    D : np.ndarray, shape (num_nodes,)
        Diffusion coefficients at each node.
    xmesh : np.ndarray, shape (num_nodes,)
        Spatial positions of the mesh points.
    xreltol : float, optional, default=0
        Relative tolerance applied to spatial positions to adjust interface locations.

    Returns
    -------
    xfull : np.ndarray, shape (3 * num_nodes,)
        Full spatial positions including center nodes and interface positions.
    Cfull : np.ndarray, shape (num_timesteps, 3 * num_nodes)
        Full concentration profiles, including node and interface values.

    Notes
    -----
    - This function enforces periodic boundary conditions by shifting indices in `C` and `k`.
    - Concentrations at the interfaces (`Ce` and `Cw`) are computed using the finite volume approach.
    - The result `Cfull` contains interleaved values: [Cw, C, Ce] for each node.

    Example
    -------
    ```python
    xfull, Cfull = compute_fc_profile_PBC(C, t, de, dw, he, hw, k, D, xmesh)
    ```
    """

    num_nodes, num_timesteps = C.shape  # Extract dimensions

    # Pre-calculate shifted indices for periodic BC
    east_shift = np.roll(np.arange(num_nodes), -1)  # Shift left (next node)
    west_shift = np.roll(np.arange(num_nodes), 1)   # Shift right (previous node)

    # Get shifted concentrations and diffusion coefficients
    C_east = C[east_shift, :]  # Shape (num_nodes, num_timesteps)
    C_west = C[west_shift, :]  # Shape (num_nodes, num_timesteps)
    k_east = k[east_shift][:, None]  # Make it broadcastable (num_nodes, 1)
    k_west = k[west_shift][:, None]  # Make it broadcastable (num_nodes, 1)

    # Eastern interface concentrations (vectorized)
    Ce = C - (de[:, None] * he[:, None] * ((k[:, None] / k_east) * C - C_east) / D[:, None])

    # Western interface concentrations (vectorized)
    Cw = C + (dw[:, None] * hw[:, None] * ((k_west / k[:, None]) * C_west - C) / D[:, None])

    # Create full concentration matrix with interfaces
    Cfull = np.empty((num_timesteps, 3*num_nodes),dtype=np.float64)

    # Compute positional tolerances
    xtol = np.min([np.min(de), np.min(dw)]) * xreltol
    xw = xmesh - dw + xtol  # Shifted west positions
    xe = xmesh + de - xtol  # Shifted east positions

    # Interleave values: West, Center, East
    Cfull[:, ::3] = Cw.T  # Ensure correct alignment
    Cfull[:, 1::3] = C.T
    Cfull[:, 2::3] = Ce.T

    # Create full position vector
    xfull = np.empty(3*num_nodes,dtype=np.float64)
    xfull[::3] = xw
    xfull[1::3] = xmesh
    xfull[2::3] = xe

    return xfull, Cfull

# %% test and debug
# -------------------------------------------------------------------
# Example usage (for debugging / standalone tests)
# -------------------------------------------------------------------
if __name__ == '__main__':
    from patankar.loadpubchem import migrant
    from patankar.food import nofood
    from patankar.layer import material
    from patankar.geometry import Packaging3D
    m = migrant("benzophenone")
    B =Packaging3D("bottle",body_radius=(20, "mm"),body_height=(15.5, "cm"),neck_radius=(1.2, "cm"),neck_height=(2,"cm"))
    HDPE = material('HDPE')(l=(1,'mm'),C0=0) # we assume 0 concentration in the bottle
    PVC = material('PVC cling film')(l=(150,'um'),C0=1000) # empirical concentral of 1000 ppm in the material
    P = HDPE+PVC
    storage = nofood(contacttime=(3,"months"),contacttemperature=(30,'degC'))
    S = m % storage << B >> P >> storage


    from patankar.food import ethanol, setoff, nofood
    from patankar.layer import PP
    from patankar.useroverride import useroverride
    useroverride.update(
        ntimes = 100,     # instead of 1000 (number of simulation times kept)
        nmesh = 300,       # instead of 600 (Finite-Volume resolution, number of FV nodes)
        tunit = "weeks",    # time units (can be any value s,min,days,weeks,months,years)
        lunit = "µm",      # length units (can be any value, nm, µm or um, mm,cm or even in)
        Cunit = "mg/kg",  # set concentration units instead of a.u.
        )


    medium = ethanol()
    medium.CF0 = 100 # works
    medium.update(CF0=100) # works also
    A = layer(layername="layer A",k=2,C0=0,D=1e-16)
    B = layer(layername="layer B")
    multilayer = A + B
    sol1 = senspatankar(multilayer, medium,t=(25,"days"))
    sol1.plotCF(t=np.array([3,10,14])*24*3600)
    sol1.plotCx()
    r=sol1.restart
    repr(r)

    # extend the solution for 40 days
    sol2 = sol1.resume((40,"days"))
    sol2.plotCF()
    sol2.plotCx()

    # extend the solution for 60 days from sol2
    sol3 = sol2.resume((60,"days"))
    sol3.update(name="sol3")
    sol3.plotCF()
    sol3.plotCx()

    # merge the previous solutions 1+2
    # extend the solution for 60 days from sol12=sol1+sol2
    sol12 = sol1+sol2
    sol123a = sol12.resume((60,"days"))
    sol123a.update(name="sol123a")
    sol123a.plotCF()
    sol123a.plotCx()

    # concat
    sol123a_ = sol12 + sol123a
    sol123a_.update(name="sol123a_ (full): sol12 + sol123a")
    sol123a_.plotCF()

    # compare with sol1+sol2+sol3
    sol123_ = sol1+sol2+sol3
    sol123_.update(name="sol123_ (full): sol1+sol2+sol3")
    sol123_.plotCF()
    sol123_.plotCx()

    # simulation of setoff
    packstorage = setoff(contacttime=(100,"days"))
    A = PP(l=(500,"um"),C0=0)
    B = PP(l=(300,"um"),C0=5000)
    AB = A+B
    print(medium)
    solAB = senspatankar(AB,packstorage)
    solAB.plotCx()

    # we extend the previous solution by putting medium in contact
    solABext = solAB.resume(medium=medium)
    solABext.plotCF()
    solABext.plotCx()
