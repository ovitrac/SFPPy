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
from migration import senspatankar
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

    from food import ethanol
    from layer import layer

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

@version: 1.0
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2022-01-17
@rev: 2025-02-24

"""
# Dependencies
import os
import random
import re
from datetime import datetime
from copy import deepcopy as duplicate
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, coo_matrix
from scipy.interpolate import interp1d
from scipy.integrate import simpson, cumulative_trapezoid
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
import pandas as pd

# Local dependencies
from layer import layer, check_units
from food import foodphysics,foodlayer

# Plot configuration (preferred units)
plotconfig = {
    "tscale": 24 * 3600, # days used as time scale
    "tunit": "days",
    "lscale": 1e-6, # µm
    "lunit": "µm",
    "Cscale": 1,
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
    return isinstance(fig, Figure) and plt.fignum_exists(fig.number)

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


# Define PrintableFigure class
class PrintableFigure(Figure):
    """Custom Figure class with print methods."""

    def print(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_figure(self, filename, destinationfolder, overwrite, dpi)

    def print_png(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_png(self, filename, destinationfolder, overwrite, dpi)

    def print_pdf(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_pdf(self, filename, destinationfolder, overwrite, dpi)

# ✅ Override `plt.figure()` and `plt.subplots()` to always use PrintableFigure
original_plt_figure = plt.figure
original_plt_subplots = plt.subplots

def custom_plt_figure(*args, **kwargs):
    """Ensure all figures are PrintableFigure."""
    kwargs.setdefault("FigureClass", PrintableFigure)
    return original_plt_figure(*args, **kwargs)

def custom_plt_subplots(*args, **kwargs):
    """Ensure plt.subplots() returns a PrintableFigure."""
    kwargs.setdefault("FigureClass", PrintableFigure)
    fig, ax = original_plt_subplots(*args, **kwargs)
    return fig, ax

# Apply overrides
plt.figure = custom_plt_figure
plt.subplots = custom_plt_subplots
plt.rcParams['figure.figsize'] = (8, 6)  # Optional default size


# %% Generic Classes to manipulate results
class Cprofile:
    """
    Class to store and interpolate a concentration profile (C(x)).
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


    """

    def __init__(self, name, description, ttarget, t, C, CF, fc, f, x, Cx, tC, C0eq, timebase,restart,xi,Cxi,_plotconfig=None):
        """constructor using positional arguments"""
        # xi and Cxi are close to x and Cx but they can be interpolated
        # their values are not saved but used to built Cprofile
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
        # Interpolated CF at ttarget
        self.interp_CF = interp1d(t, CF, kind="linear", fill_value="extrapolate")
        self.CFtarget = self.interp_CF(ttarget)
        # Interpolated concentration profile at ttarget
        self.interp_Cx = interp1d(t, Cx.T, kind="linear", axis=1, fill_value="extrapolate")
        self.Cxtarget = self.interp_Cx(ttarget)
        # Restart information including restults at ttarget
        # xi and Cxi are available only from a fresh simulation
        # these data are missing from operation +, in this case we use restart as supplied
        if xi is not None and Cxi is not None:
            Cxi_interp = interp1d(t, Cxi.T, kind="linear", axis=1, fill_value="extrapolate")
            Cxi_at_t = Cxi_interp(ttarget) # this profile has inceasing xi using xreltol
            restart.freezeCF(ttarget,self.CFtarget)
            restart.freezeCx(xi,Cxi_at_t)
        self.restart = restart
        if _plotconfig is None:
            self._plotconfig = plotconfig # fresh simulation
        else:
            self._plotconfig = _plotconfig # if from an existing SensPatankarResult
        # for simulation chaining
        self.savestate(self.restart.inputs["multilayer"],self.restart.inputs["medium"])

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
        newsol = senspatankar(multilayer=inputs["multilayer"],
                              medium=newmedium,
                              name=kwargs.get("name",inputs["name"]),
                              description=kwargs.get("description",inputs["description"]),
                              t=t,
                              autotime=kwargs.get("autotime",inputs["autotime"]),
                              timescale=kwargs.get("timescale",inputs["timescale"]),
                              Cxprevious=newCx0,
                              ntimes=kwargs.get("ntimes",inputs["ntimes"]),
                              RelTol=kwargs.get("RelTol",inputs["RelTol"]),
                              AbsTol=kwargs.get("AbsTol",inputs["AbsTol"]))
        return newsol

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
            xi=None,  # xi and Cxi values are available
            Cxi=None  # only from a fresh simulation
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
            f'CF({(self.ttarget / plotconfig["tscale"]).item():.4g} [{plotconfig["tunit"]}]) = '
            f'{(self.CFtarget / plotconfig["Cscale"]).item():.4g} [{plotconfig["Cunit"]}]>')



    def plotCF(self, t=None, trange=None):
        """
        Plot the concentration in the food (CF) as a function of time and highlight the target time(s).

        Parameters
        ----------
        t : float, list, or None, optional
            Specific time(s) for which the concentration should be highlighted.
            If None, defaults to `ttarget`.
        trange : None, float, or list [t_min, t_max], optional
            If None, the full profile is shown.
            If a float, it is treated as an upper bound (lower bound assumed 0).
            If a list `[t_min, t_max]`, the profile is interpolated between these values.
        """

        # extract plotconfig
        plotconfig = self._plotconfig

        # Ensure t is a list (even if a single value is given)
        if t is None:
            t_values = [self.ttarget]
        elif isinstance(t, (int, float)):
            t_values = [t]
        elif isinstance(t,np.ndarray):
            t_values = t.flatten()
        elif isinstance(t,tuple):
            t_values = check_units(t)[0]
        else:
            t_values = np.array(t)  # Convert to array

        # Interpolate CF values at given times
        CF_t_values = self.interp_CF(t_values)

        # Handle trange interpolation
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

            # Generate interpolated time values
            t_plot = np.linspace(t_min, t_max, 500)
            CF_plot = self.interp_CF(t_plot)  # Interpolated CF values

        # Set up colormap for multiple t values
        cmap = plt.get_cmap('viridis', len(t_values))
        norm = mcolors.Normalize(vmin=min(t_values), vmax=max(t_values))

        # Create the figure
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot CF curve (either original or interpolated)
        ax.plot(t_plot / plotconfig["tscale"], CF_plot / plotconfig["Cscale"],
                label='Concentration in Food', color='b')

        # Highlight each target time
        for i, tC in enumerate(t_values):
            color = tooclear(cmap(norm(tC))) if len(t_values) > 1 else 'r'  # Use color map only if multiple t values

            # Vertical and horizontal lines
            ax.axvline(tC / plotconfig["tscale"], color=color, linestyle='--', linewidth=1)
            ax.axhline(CF_t_values[i] / plotconfig["Cscale"], color=color, linestyle='--', linewidth=1)

            # Intersection point
            ax.scatter(tC / plotconfig["tscale"], CF_t_values[i] / plotconfig["Cscale"],
                       color=color, zorder=3)

            # Annotate time
            ax.text(tC / plotconfig["tscale"], min(CF_plot) / plotconfig["Cscale"],
                    f'{(tC / plotconfig["tscale"]).item():.2f} {plotconfig["tunit"]}',
                    verticalalignment='bottom', horizontalalignment='right', rotation=90, fontsize=10, color=color)

            # Annotate concentration
            ax.text(min(t_plot) / plotconfig["tscale"], CF_t_values[i] / plotconfig["Cscale"],
                    f'{(CF_t_values[i] / plotconfig["Cscale"]).item():.2f} {plotconfig["Cunit"]}',
                    verticalalignment='bottom', horizontalalignment='left', fontsize=10, color=color)

        # Labels and title
        ax.set_xlabel(f'Time [{plotconfig["tunit"]}]')
        ax.set_ylabel(f'Concentration in Food [{plotconfig["Cunit"]}]')
        title_main = "Concentration in Food vs. Time"
        title_sub = rf"$\bf{{{self.name}}}$" + (f": {self.description}" if self.description else "")
        ax.set_title(f"{title_main}\n{title_sub}", fontsize=10)
        ax.text(0.5, 1.05, title_sub, fontsize=8, ha="center", va="bottom", transform=ax.transAxes)
        ax.set_title(title_main)
        ax.legend()
        ax.grid(True)
        plt.show()
        # store metadata
        setattr(fig,_fig_metadata_atrr_,f"pltCF_{self.name}")
        return fig



    def plotCx(self, t=None, nmax=15):
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
        """

        # extract plotconfig
        plotconfig = self._plotconfig


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
        norm = mcolors.Normalize(vmin=t_values.min(), vmax=t_values.max()) if len(t_values) > 1 else mcolors.Normalize(vmin=self.t.min(), vmax=self.t.max())
        cmap = plt.get_cmap('viridis', nmax)  # 'viridis' is similar to Parula

        fig, ax = plt.subplots(figsize=(8, 6))  # Explicitly create a figure and axis

        # Plot all valid concentration profiles with time-based colormap
        for tC in t_values:
            C = self.interp_Cx(tC)
            color = tooclear(cmap(norm(tC)))  # Get color from colormap
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
        title_main = "Concentration Profiles in Packaging vs. Position"
        title_sub = rf"$\bf{{{self.name}}}$" + (f": {self.description}" if self.description else "")
        ax.set_title(f"{title_main}\n{title_sub}", fontsize=10)
        ax.text(0.5, 1.05, title_sub, fontsize=8, ha="center", va="bottom", transform=ax.transAxes)
        ax.set_title(title_main)
        ax.grid(True)
        ax.legend()
        plt.show()
        # store metadata
        setattr(fig,_fig_metadata_atrr_,f"pltCx_{self.name}")
        return fig


# Container for multiple simulations
class CFSimulationContainer:
    """
    Container to store and compare multiple CF results from different simulations.

    Attributes
    ----------
    curves : dict
        Stores CF results with unique keys. Each entry contains:
        - 'label': Label used for legend.
        - 'tmin', 'tmax': Time range of the simulation.
        - 'interpolant': Interpolated CF function.
        - 'color': Assigned color for plotting.
        - 'linestyle': Line style (default is '-').
        - 'linewidth': Line width (default is 2).
    """

    def __init__(self,name="",description=""):
        """Initialize an empty container for CF results."""
        self.curves = {}
        self._name = name
        self._description = description
        self._plotconfig = plotconfig

    @property
    def name(self):
        return self._name or autoname(6)

    @property
    def description(self):
        return self._description or f"comparison of {len(self.curves)} curves"


    def add(self, simulation_result, label=None, color=None, linestyle="-", linewidth=2):
        """
        Add a new CF result to the container.

        Parameters
        ----------
        simulation_result : SensPatankarResult
            The simulation result to store.
        label : str, optional
            Label for the curve (used in legend). Defaults to 'plot1', 'plot2', etc.
        color : str or tuple, optional
            Color for the curve. Defaults to automatic colormap.
        linestyle : str, optional
            Line style for the plot (e.g., '-', '--', '-.', ':'). Default is '-'.
        linewidth : float, optional
            Line width for the plot. Default is 2.
        """

        # Check input
        if not isinstance(simulation_result,SensPatankarResult):
            raise TypeError(f"simulation_result should be a SensPatankarResult object not a {type(simulation_result).__name__}")

        # Auto-generate a label if not provided
        label = label or self.name
        label = label or f"plot{len(self.curves) + 1}"

        # Generate a unique key (first 40 letters of label)
        key = label[:40] # 40 max

        # Create the interpolant function
        interpolant = interp1d(simulation_result.t, simulation_result.CF,
                               kind="linear", fill_value="extrapolate", bounds_error=False)

        # Assign color from a colormap if not provided
        if color is None:
            cmap = cm.get_cmap("tab10", len(self.curves) + 1)
            color = cmap(len(self.curves) % 10)  # Cycle through colors

        # Store the result in the dictionary (replacing existing key if necessary)
        self.curves[key] = {
            "label": label,
             "name": simulation_result.name,
      "description": simulation_result.description,
             "tmin": simulation_result.t.min(),
             "tmax": simulation_result.t.max(),
      "interpolant": interpolant,
            "color": color,
        "linestyle": linestyle,
        "linewidth": linewidth
        }

    def delete(self, label):
        """
        Remove a stored curve by its label.

        Parameters
        ----------
        label : str
            Label of the curve to delete.
        """
        key = label[:10]
        if key in self.curves:
            del self.curves[key]
            print(f"Deleted curve '{label}'")
        else:
            print(f"No curve found with label '{label}'")

    def __repr__(self):
        """Return a summary of stored CF curves."""
        if not self.curves:
            return "<CFSimulationContainer: No stored curves>"

        repr_str = "<CFSimulationContainer: Stored CF Curves>\n"
        repr_str += "--------------------------------------------------\n"
        for key, data in self.curves.items():
            repr_str += (f"[{data['label']}] "
                         f"Name: {data['name']} | "
                         f"Descr: {data['description']} | "
                         f"Time: [{data['tmin']:.2e}, {data['tmax']:.2e}] s | "
                         f"Color: {data['color']} | "
                         f"Style: {data['linestyle']} | "
                         f"Width: {data['linewidth']}\n")
        return repr_str


    def plotCF(self, t_range=None):
        """
        Plot all stored CF curves in a single figure.

        Parameters
        ----------
        t_range : tuple (t_min, t_max), optional
            Time range for plotting. If None, uses each curve's own range.
        plotconfig : dict, optional
            Dictionary with plotting configuration, containing:
            - "tunit": Time unit label (e.g., 's').
            - "Cunit": Concentration unit label (e.g., 'mg/L').
            - "tscale": Time scaling factor.
            - "Cscale": Concentration scaling factor.
        """

        # extract plotconfig
        plotconfig = self._plotconfig

        if not self.curves:
            print("No curves to plot.")
            return

        fig, ax = plt.subplots(figsize=(8, 6))

        for data in self.curves.values():
            # Determine the time range
            t_min, t_max = data["tmin"], data["tmax"]
            if t_range:
                t_min, t_max = max(t_min, t_range[0]), min(t_max, t_range[1])

            # Generate interpolated values
            t_plot = np.linspace(t_min, t_max, 500)
            CF_plot = data["interpolant"](t_plot)

            # Apply unit scaling if plotconfig is provided
            if plotconfig:
                t_plot /= plotconfig.get("tscale", 1)
                CF_plot /= plotconfig.get("Cscale", 1)

            # Plot the curve
            ax.plot(t_plot, CF_plot, label=data["label"],
                    color=data["color"], linestyle=data["linestyle"], linewidth=data["linewidth"])

        # Configure the plot
        ax.set_xlabel(f'Time [{plotconfig["tunit"]}]' if plotconfig else "Time")
        ax.set_ylabel(f'Concentration in Food [{plotconfig["Cunit"]}]' if plotconfig else "CF")
        title_main = "Concentration in Food vs. Time"
        title_sub = rf"$\bf{{{self.name}}}$" + (f": {self.description}" if self.description else "")
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
        rgb()


# restartfile
class restartfile:
    """
    Containter for the restartfile
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
    Containter for the restartfile
    """
    def __init__(self,multilayer,medium,name,description,
                 t,autotime,timescale,Cxprevious,
                 ntimes,RelTol,AbsTol):
        """constructor to be called at the intialization"""
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
                 ntimes=1e3, RelTol=1e-6, AbsTol=1e-6):
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

        from food import ethanol
        from layer import layer
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

    # Refresh the physics of medium for parameters tunned by the end-user
    medium.refresh()

    # extract the PBC flag (True for setoff)
    PBC = medium.PBC

    # Restart file initialization (all parameters are saved)
    restart = restartfile_senspantakar(multilayer, medium, name,
            description, t, autotime, timescale, Cxprevious, ntimes, RelTol, AbsTol)

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
        xi=xfulli*l_ref, # for restart only
        Cxi=Cfull_dimlessi*C0eq # for restart only
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
    xfull = np.empty(3 * num_nodes)
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
    Cfull_dimless = np.empty((num_timesteps, 3 * num_nodes))
    Cfull_dimless[:, ::3] = Cw.T      # Every 3rd column is Cw
    Cfull_dimless[:, 1::3] = C_dimless.T  # Every 3rd column (offset by 1) is C
    Cfull_dimless[:, 2::3] = Ce.T      # Every 3rd column (offset by 2) is Ce

    return xfull, Cfull_dimless


def compute_fc_profile_PBC(C, t, de, dw, he, hw, k, D, xmesh, xreltol=0):
    """Calculate interface concentrations with periodic boundary conditions"""

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
    Cfull = np.empty((num_timesteps, 3*num_nodes))

    # Compute positional tolerances
    xtol = np.min([np.min(de), np.min(dw)]) * xreltol
    xw = xmesh - dw + xtol  # Shifted west positions
    xe = xmesh + de - xtol  # Shifted east positions

    # Interleave values: West, Center, East
    Cfull[:, ::3] = Cw.T  # Ensure correct alignment
    Cfull[:, 1::3] = C.T
    Cfull[:, 2::3] = Ce.T

    # Create full position vector
    xfull = np.empty(3*num_nodes)
    xfull[::3] = xw
    xfull[1::3] = xmesh
    xfull[2::3] = xe

    return xfull, Cfull

# %% test and debug
# -------------------------------------------------------------------
# Example usage (for debugging / standalone tests)
# -------------------------------------------------------------------
if __name__ == '__main__':
    from food import ethanol, setoff, nofood
    from layer import PP

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
