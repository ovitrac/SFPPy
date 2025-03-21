#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Geometry
===============================================================================
Provides a framework for defining 3D packaging geometries and their surface/volume properties.
Supports standalone shapes and composite structures for packaging simulations.

**Main Components:**
- **Shape3D** (Base class for 3D shapes)
    - Subclasses: `Cylinder`, `Cone`, `RectangularPrism`, `Sphere`, `SquarePyramid`, `Hemisphere`
    - Implements `_compute_volume()`, `_compute_surface_area()`, `_compute_connectors()`
- **CompositeShape** (Combines multiple shapes while adjusting for overlapping volumes and shared faces)
- **Packaging3D** (High-level interface for working with packaging geometries)
- **Connector** (Defines interfaces between shapes for composite construction)
- **SHAPE_REGISTRY** (Maps packaging terms like "can" or "bottle" to their geometric models)

**Integration with SFPPy Modules:**
- Uses `check_units()` from `layer.py` to convert dimensions into SI units.
- Supports `migration.py` by providing accurate volume and surface area computations for mass transfer models.

Example:
```python
from patankar.geometry import Packaging3D
pkg = Packaging3D('bottle', body_radius=(5, 'cm'), body_height=(20, 'cm'))
vol, area = pkg.get_volume_and_area()
```


===============================================================================
Details
===============================================================================
Purpose:
  This module provides a framework for defining and combining various
  three-dimensional packaging shapes—both simple (e.g., cylinders, cones,
  rectangular prisms) and composite (e.g., 'bottle' = large cylinder + narrow
  cylinder). It also calculates each shape’s internal volume (in m³) and
  internal surface area (in m²).

Overview:
  1. **Shape3D and Subclasses**:
     - Each subclass (Cylinder, Cone, RectangularPrism, Sphere, SquarePyramid,
       Hemisphere) implements:
         * `_compute_volume()` and `_compute_surface_area()`.
         * `_compute_connectors()`, which returns a list of `Connector` objects
           representing the flat faces that can potentially connect to other
           shapes. Connectors have a face area and an axis (normal vector).
     - Connectors allow shapes to “snap” together in a `CompositeShape`.

  2. **CompositeShape**:
     - Manages multiple shapes added together, summing volumes and surface
       areas, minus overlaps along shared connector faces.
     - Uses `add_shape(...)` to join a new shape to any existing sub-shape
       via matching connector orientations. Overlapping face area is removed
       from the total surface area calculation.

  3. **Synonyms and Shape Registry**:
     - A dictionary `SHAPE_REGISTRY` maps real-world packaging names (like
       "can", "box", "glass") to their corresponding Shape3D classes.
     - Some "synonyms" map to more complex, composite constructs. For example,
       the name "bottle" creates a `CompositeShape` of two cylinders (body +
       neck).

  4. **Units**:
     - Dimensions can be given either as floats in meters or as `(value, "unit")`
       pairs. The helper `_to_m(...)` uses `check_units` (imported from
       `layer`) to convert numeric values to SI units (meters).

  5. **Packaging3D**:
     - A high-level interface for creating either a single shape or a
       composite shape by name and keyword arguments. It returns volume
       (in m³) and surface area (in m²) via `.get_volume_and_area()`.

Usage Example:
    from patankar.packaging import Packaging3D

    # Create a 'bottle' (two stacked cylinders) by specifying body and neck dims
    pkg = Packaging3D(
        'bottle',
        body_radius=(5, 'cm'),
        body_height=(20, 'cm'),
        neck_radius=(2, 'cm'),
        neck_height=(5, 'cm')
    )
    vol, area = pkg.get_volume_and_area()
    print("Volume (m^3):", vol)
    print("Surface Area (m^2):", area)

    # Create a single shape (e.g., 'can' which is a cylinder) with radius
    # and height specified in centimeters
    pkg2 = Packaging3D('can', radius=(4, 'cm'), height=(12, 'cm'))
    vol2, area2 = pkg2.get_volume_and_area()

About units:
    All lengths can be given:
        - without units: all lengths are assumed to be in meters (i.e., their SI unit)
        - with units by using a tupple (value,"unit"), where unit can be m,dm,cm,mm,um,nm...
    Input units can be heterogenerous, the result is always SI:
        - [m**2] for surface areas
        - [m**3] for volumes

Get help and syntax:
    help_geometry()

Notes:
  - This code is primarily illustrative. In a production system, you may expand
    the geometry classes, refine orientation logic for connectors, handle
    partial overlaps, or integrate with a 3D transform library for more
    sophisticated shape placement.
  - The overlap deduction for composite shapes is simplified. It subtracts
    `2 * (minimum overlapping face area)` from the total surface area, which
    assumes a perfect “face-to-face” join.



Dependencies:
  - Python 3.x
  - `check_units` function from the `layer` module


@version: 1.0
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-10-28, rev. 2025-03-13

===============================================================================
"""


# %% Dependencies
import math
import numpy as np
from collections import defaultdict
from patankar.layer import check_units

__all__ = ['Ball', 'Box', 'CompositeShape', 'Con', 'Cone', 'Connector', 'Cuboid', 'Cyl', 'Cylinder', 'HalfSphere', 'HemiSphere', 'Hemisphere', 'OpenBox', 'OpenBox2', 'OpenCon', 'OpenCone', 'OpenCyl', 'OpenCylinder1', 'OpenCylinder2', 'OpenPrism1', 'OpenPrism2', 'OpenSquare1', 'OpenSquare2', 'OpenSquareBox', 'OpenSquareBox2', 'Packaging3D', 'Pyramid', 'RectangularPrism', 'Shape3D', 'Sphere', 'SquarePyramid', 'Tube', 'check_units', 'create_packaging_widget', 'create_shape_by_name', 'get_all_shapes_info', 'get_geometries_and_synonyms', 'help_geometry']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.32"


# %% Widget
def create_packaging_widget(default_shape="cylinder", default_instance_name="shape1"):
    """
    Creates a widget interface to define a 3D packaging geometry.

    The user:
      - Picks a shape from a dropdown list (options include: 'cylinder', 'cone',
        'rectangular_prism', 'sphere', 'square_pyramid', 'hemisphere', 'cube', 'bottle').
      - Once a shape is selected, the interface displays input fields for each required
        dimension. For each parameter the user enters a value (via a FloatText) and selects
        a unit (from: µm, mm, cm, dm, in). For example, a Cylinder requires 'radius' and 'height'.
      - A short description (documentation and synonyms) for the selected shape is also shown.
      - The user enters an instance name.
      - When the "Create Packaging" button is clicked, the widget collects the input values,
        constructs a dictionary of dimensions (each as a tuple (value, unit)), and instantiates
        Packaging3D(selected_shape, **dimensions). The new instance is stored in a global
        dictionary (mypackaging) using the provided instance name.

    Special cases:
      - For "bottle", the required dimensions are: body_radius, body_height, neck_radius, neck_height.
      - For "cube", the only required dimension is "side" (which is used for all three dimensions).

    Returns:
      An ipywidgets.VBox instance containing the UI.

    Usage Example (in a Jupyter cell):
      from IPython.display import display
      pkg_widget = create_packaging_widget()
      display(pkg_widget)
      # Later, access the created packaging instances via the global dict 'mypackaging'
      import builtins
      print(builtins.mypackaging)
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError("ipywidgets and IPython are required for this interface.") from e

    import builtins
    if not hasattr(builtins, "mypackaging"):
        builtins.mypackaging = {}
    global mypackaging
    mypackaging = builtins.mypackaging

    # flag for preheated gui interface (widgets should be initialized manually, instead of being empty)
    _preheatedGUI_ = hasattr(builtins, "_PREHEATED_") and getattr(builtins, "_PREHEATED_") is True

    # Define the list of available shapes.
    # (You can later expand this list if needed.)
    shape_options = ["cylinder", "cone", "rectangular_prism", "sphere",
                     "square_pyramid", "hemisphere", "cube", "bottle"]
    shape_dropdown = widgets.Dropdown(
         options=shape_options,
         value=default_shape if default_shape in shape_options else shape_options[0],
         description="Shape:",
         layout=widgets.Layout(width="50%")
    )

    # Container for the dynamic parameter input fields.
    params_box = widgets.VBox([])

    # Label for shape documentation.
    doc_label = widgets.HTML(value="")

    # Text input for instance name.
    instance_name_input = widgets.Text(
         value=default_instance_name,
         description="Instance Name:",
         layout=widgets.Layout(width="50%")
    )

    # Button to create the Packaging3D instance.
    create_btn = widgets.Button(
         description="Create Packaging",
         button_style="success",
         tooltip="Click to create Packaging3D instance"
    )

    # Output area to show messages.
    output = widgets.Output()

    # Mapping to get canonical class names from the dropdown value.
    # For non-special shapes, we assume the first letter capitalized.
    canonical_map = {
         "cylinder": "Cylinder",
         "cone": "Cone",
         "rectangular_prism": "RectangularPrism",
         "sphere": "Sphere",
         "square_pyramid": "SquarePyramid",
         "hemisphere": "Hemisphere",
         "cube": "cube",    # special case
         "bottle": "bottle" # special composite case
    }

    # Function to update the parameter fields when shape selection changes.
    def update_params(*args):
        params_children = []
        selected_shape = shape_dropdown.value.lower()
        # For special cases, define required parameters explicitly.
        if selected_shape == "bottle":
            required_params = ["body_radius", "body_height", "neck_radius", "neck_height"]
            doc_text = ("<b>Bottle</b>: A composite shape built from two cylinders (body and neck).<br>"
                        "Required parameters: body_radius, body_height, neck_radius, neck_height.")
        elif selected_shape == "cube":
            required_params = ["side"]
            doc_text = ("<b>Cube</b>: A cube with equal dimensions (side).")
        else:
            # For basic shapes, try to get info from get_all_shapes_info().
            canonical = canonical_map.get(selected_shape, selected_shape.capitalize())
            info = get_all_shapes_info().get(canonical, None)
            if info is not None:
                required_params = info.get("required_params", [])
                doc_text = ("<b>" + canonical + "</b>: " +
                            info.get("doc", "No documentation available."))
            else:
                # Fallback: assume no required parameters.
                required_params = []
                doc_text = "No documentation available for this shape."
        doc_label.value = doc_text
        # Create dynamic widgets for each required parameter.
        param_widgets = []
        for param in required_params:
            # For each parameter, we create a FloatText for value and a Dropdown for unit.
            value_widget = widgets.FloatText(value=5.0, description=param + ":",
                                               layout=widgets.Layout(width="40%"))
            unit_widget = widgets.Dropdown(
                options=["µm", "mm", "cm", "dm", "in"],
                value="cm",
                description="Unit:",
                layout=widgets.Layout(width="20%")
            )
            param_widgets.append(widgets.HBox([value_widget, unit_widget]))
            # We attach the parameter name as an attribute for later retrieval.
            value_widget.param_name = param
            unit_widget.param_name = param
        params_box.children = param_widgets

    # Update parameters when shape selection changes.
    shape_dropdown.observe(update_params, names="value")
    update_params()

    # Callback when "Create Packaging" is clicked.
    def create_packaging_instance(b):
        with output:
            output.clear_output()
            selected_shape = shape_dropdown.value.lower()
            dims = {}
            # Collect parameter values from the params_box.
            for row in params_box.children:
                # row is an HBox with two children: a FloatText and a Dropdown.
                if len(row.children) < 2:
                    continue
                value_widget = row.children[0]
                unit_widget = row.children[1]
                param = value_widget.param_name
                dims[param] = (value_widget.value, unit_widget.value)
            # Create the Packaging3D instance.
            try:
                pkg = Packaging3D(selected_shape, **dims)
            except Exception as e:
                print("Error creating packaging shape:", e)
                return
            inst_name = instance_name_input.value.strip()
            if not inst_name:
                print("Please provide an instance name.")
                return
            builtins.mypackaging[inst_name] = pkg
            print(f"Packaging instance '{inst_name}' created:")
            print(pkg)
            # Optionally, also show volume and surface area.
            vol, area = pkg.get_volume_and_area()
            print(f"Volume: {vol.item():.4g} m³, Surface Area: {area.item():.4g} m²")
            print("\nCurrent packaging instances:", list(builtins.mypackaging.keys()))

    create_btn.on_click(create_packaging_instance)

    if _preheatedGUI_:
        create_packaging_instance(None) # we instantiate manually

    # Arrange the UI elements.
    ui = widgets.VBox([
         shape_dropdown,
         doc_label,
         params_box,
         instance_name_input,
         create_btn,
         output
    ])

    return ui


# %% Helper functions

# Convert lengths to SI units [m]
def _to_m(value):
    """
    Convert a dimension value to meters using check_units if it's a tuple.
    Otherwise assume the value is already in meters.
    """
    if isinstance(value, tuple):
        val_in_m, _ = check_units(value)  # check_units returns (value_in_SI, "m")
        return val_in_m
    else:
        return value

# %% Base Geometry Classes
class Connector:
    """
    Represents a 'connection face' on a shape:
      - area: the connectable area (m^2)
      - axis: a unit vector (tuple) indicating the orientation of the connector
      - name: optionally label the connector (e.g. 'top', 'bottom', etc.)
    """
    def __init__(self, area, axis=(0, 0, 1), name=""):
        self.area = area
        # Normalize axis for safety
        mag = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        if mag > 0:
            self.axis = (axis[0]/mag, axis[1]/mag, axis[2]/mag)
        else:
            self.axis = axis
        self.name = name

    def __repr__(self):
        """String representation of the Connector object."""
        axis_str = f"({self.axis[0]:.2f}, {self.axis[1]:.2f}, {self.axis[2]:.2f})"
        name_str = f"'{self.name}'" if self.name else "(unnamed)"
        print(f"Connector(name={name_str}, area={self.area.item():.4g} m², axis={axis_str})")
        return str(self)

    def __str__(self):
        """Formatted representation of the connector"""
        return f"<{self.__class__.__name__}: {self.name}>"


class Shape3D:
    """
    Base class for a 3D shape.
    Subclasses must implement:
      - _compute_volume()
      - _compute_surface_area()
      - _compute_connectors() -> list of Connector objects

    Additionally, each subclass should define the class attribute
      expected_dimensions: dict[str, list[str]]
    mapping the canonical dimension names to a list of acceptable synonyms.
    """
    # In the base class we keep an empty dictionary. Each subclass must override it.
    expected_dimensions = {}

    def __init__(self, **dimensions):
        # Resolve synonyms and check that every mandatory dimension is provided.
        resolved_dims = {}
        for canonical, synonyms in self.expected_dimensions.items():
            provided = None
            # Check all provided keywords (case-insensitive)
            for syn in synonyms:
                for key, value in dimensions.items():
                    if key.lower() == syn.lower():
                        provided = value
                        break
                if provided is not None:
                    break
            if provided is None:
                # Form a clear message indicating what is missing and the allowed options.
                raise ValueError(
                    f"Missing dimension for '{canonical}'. Expected one of: {synonyms}"
                )
            resolved_dims[canonical] = _to_m(provided)
        self.dimensions = resolved_dims

    def volume(self):
        return self._compute_volume()

    def surface_area(self):
        return self._compute_surface_area()

    def connectors(self):
        return self._compute_connectors()

    def _compute_volume(self):
        raise NotImplementedError

    def _compute_surface_area(self):
        raise NotImplementedError

    def _compute_connectors(self):
        return []

    @classmethod
    def list_expected_keywords(cls):
        msg = f"Expected keywords for {cls.__name__}:"
        for canonical, synonyms in cls.expected_dimensions.items():
            msg += f"\n  - {canonical}: {synonyms}"
        return msg

    def __repr__(self):
        class_name = self.__class__.__name__
        dimensions_str = ", ".join(
            f"{k}={v.item():.4g} m" if isinstance(v, np.ndarray) else f"{k}={v:.4f} m"
            for k, v in self.dimensions.items()
        )
        vol = self.volume()
        surf = self.surface_area()
        connectors = self.connectors()
        connector_str = (
            "\n  - ".join(repr(c) for c in connectors) if connectors else "None"
        )
        print(
            f"{class_name}(\n"
            f"  Dimensions: {dimensions_str}\n"
            f"  Volume: {vol.item():.4g} m³\n"
            f"  Surface Area: {surf.item():.4g} m²\n"
            f"  Connectors:\n  - {connector_str}\n"
            f")"
        )
        return str(self)

    def __str__(self):
        n = len(self.connectors())
        return f"<{self.__class__.__name__} with {n} connector{'s' if n>1 else ''}>"

# ----------------------------------------------------------------------------
# Basic shapes with synonyms defined
# ----------------------------------------------------------------------------

# For Cylinder, Cone, OpenCylinder1, OpenCylinder2, OpenCone:
common_rad_height = {
    "radius": ["radius", "r"],
    "height": ["height", "length", "h", "l"]
}

class Cylinder(Shape3D):
    """
    A cylinder with radius=r and height=h.
    Has two connectors (top and bottom).
    """
    expected_dimensions = common_rad_height

    def _compute_volume(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        return math.pi * r**2 * h

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        # Full cylinder: side + 2 ends
        return 2.0 * math.pi * r * h + 2.0 * math.pi * r**2

    def _compute_connectors(self):
        """
        Two circular faces: top (normal +z), bottom (normal -z).
        """
        r = self.dimensions['radius']
        area_face = math.pi * r**2
        c_top = Connector(area=area_face, axis=(0,0,1), name="cylinder_top")
        c_bottom = Connector(area=area_face, axis=(0,0,-1), name="cylinder_bottom")
        return [c_top, c_bottom]


class Cone(Shape3D):
    """
    A cone with radius=r, height=h.
    Typically only 1 connectable face: the circular base (normal -z).
    """
    expected_dimensions = common_rad_height

    def _compute_volume(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        return (1.0/3.0) * math.pi * r**2 * h

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        slant = math.sqrt(r**2 + h**2)
        base_area = math.pi * r**2
        lateral_area = math.pi * r * slant
        return base_area + lateral_area

    def _compute_connectors(self):
        r = self.dimensions['radius']
        area_face = math.pi * r**2
        # We'll define the base as normal -z
        return [Connector(area=area_face, axis=(0,0,-1), name="cone_base")]


class OpenCone(Shape3D):
    """
    A cone with the base removed, leaving a single open circular face.

    Required dimensions:
      radius, height
    Volume:
      Same as a full cone => (1/3)*π*r^2*h
    Surface area:
      Only the lateral surface => π*r*sqrt(r^2 + h^2)
      (No base area since it's open.)
    Connectors:
      One connector at the base (the open circle).
    """
    expected_dimensions = common_rad_height

    def _compute_volume(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        return (1.0/3.0) * math.pi * (r**2) * h

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        slant = math.sqrt(r**2 + h**2)
        # Lateral area only
        return math.pi * r * slant

    def _compute_connectors(self):
        r = self.dimensions['radius']
        area_face = math.pi * r**2
        # The open face is the base, normal -z
        return [Connector(area=area_face, axis=(0,0,-1), name="open_cone_base")]


class OpenCylinder1(Shape3D):
    """
    An open cylinder with exactly one open end (like a glass, pot, or jar).

    Volume:
      π * r^2 * h
    Surface area:
      Lateral area (2πrh) + base area (πr^2) => 2πrh + πr^2
    Connectors:
      Only one at the bottom (circular face).
    """
    expected_dimensions = common_rad_height

    def _compute_volume(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        return math.pi * r**2 * h

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        lateral_area = 2.0 * math.pi * r * h
        bottom_area = math.pi * r**2
        return lateral_area + bottom_area

    def _compute_connectors(self):
        r = self.dimensions['radius']
        bottom_area = math.pi * r**2
        return [Connector(area=bottom_area, axis=(0, 0, -1), name="open_cylinder1_bottom")]


class OpenCylinder2(Shape3D):
    """
    An open cylinder with two open ends (like a straw or tube).

    Volume:
      π * r^2 * h
    Surface area:
      Only lateral area => 2πrh
      (No top or bottom disk, since both ends are open.)
    Connectors:
      Two (top and bottom), each with area πr^2.
    """
    expected_dimensions = common_rad_height

    def _compute_volume(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        return math.pi * r**2 * h

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        h = self.dimensions['height']
        # No bases since both ends are open
        return 2.0 * math.pi * r * h

    def _compute_connectors(self):
        r = self.dimensions['radius']
        area_face = math.pi * r**2
        # top face (normal +z) and bottom face (normal -z)
        c_top = Connector(area=area_face, axis=(0,0, 1), name="open_cylinder2_top")
        c_bottom = Connector(area=area_face, axis=(0,0,-1), name="open_cylinder2_bottom")
        return [c_top, c_bottom]


# For RectangularPrism, OpenPrism1, OpenPrism2:
common_prism = {
    "length": ["length", "l"],
    "width": ["width", "w", "depth", "d"],
    "height": ["height", "h"]
}

class RectangularPrism(Shape3D):
    """
    A rectangular prism with length=l, width=w, height=h.
    Has 6 connectors for each face.
    """
    expected_dimensions = common_prism

    def _compute_volume(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        return l * w * h

    def _compute_surface_area(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        return 2.0 * (l*w + w*h + h*l)

    def _compute_connectors(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']

        # areas
        area_lw = l * w
        area_wh = w * h
        area_hl = h * l

        # Each face axis.
        # We'll define +z, -z, +y, -y, +x, -x as possible "connectors".
        return [
            Connector(area=area_lw, axis=(0,0, 1), name="top_face"),
            Connector(area=area_lw, axis=(0,0,-1), name="bottom_face"),
            Connector(area=area_wh, axis=(0, 1,0), name="front_face"),
            Connector(area=area_wh, axis=(0,-1,0), name="back_face"),
            Connector(area=area_hl, axis=( 1,0,0), name="right_face"),
            Connector(area=area_hl, axis=(-1,0,0), name="left_face")
        ]


class OpenPrism1(Shape3D):
    """
    A rectangular prism with ONE open face.

    Required dimensions:
      length, width, height
    Volume:
      length * width * height
    Surface area:
      (Sum of all faces) - area of the open face
      i.e. 2*(lw + lh + wh) - lw (assuming top is open).
      So total = lw + 2*(lh + wh).
    Connectors:
      One connector at the open face.

    By convention, let's treat the "top" (normal +z) as open.
    That means the bottom is length x width, and sides are intact.
    """
    expected_dimensions = common_prism

    def _compute_volume(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        return l * w * h

    def _compute_surface_area(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        total_closed = 2.0*(l*w + w*h + h*l)
        # Open face is the top (area = l*w).
        return total_closed - (l*w)

    def _compute_connectors(self):
        """
        The open face is the top: area = l*w, normal +z.
        """
        l = self.dimensions['length']
        w = self.dimensions['width']
        top_area = l * w
        return [Connector(area=top_area, axis=(0,0,1), name="open_prism1_top")]


class OpenPrism2(Shape3D):
    """
    A rectangular prism with TWO open faces (no top, no bottom).

    Required dimensions:
      length, width, height
    Volume:
      length * width * height
    Surface area:
      (Sum of all faces) - 2*(area of top + bottom)
      i.e. 2*(l*w + w*h + h*l) - 2*(l*w)
      = 2*(w*h + h*l)
    Connectors:
      Two connectors: top (+z), bottom (-z).
    """
    expected_dimensions = common_prism

    def _compute_volume(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        return l * w * h

    def _compute_surface_area(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        h = self.dimensions['height']
        # Full closed prism area: 2*(lw + lh + wh)
        # Remove top (lw) and bottom (lw), total of 2*lw
        return 2.0*(l*h + w*h)  # Just the 4 side faces

    def _compute_connectors(self):
        l = self.dimensions['length']
        w = self.dimensions['width']
        face_area = l * w
        top_connector = Connector(area=face_area, axis=(0,0,1),  name="open_prism2_top")
        bot_connector = Connector(area=face_area, axis=(0,0,-1), name="open_prism2_bottom")
        return [top_connector, bot_connector]


# For Sphere and Hemisphere:
common_sphere = {
    "radius": ["radius", "r"]
}

class Sphere(Shape3D):
    """
    A sphere with radius=r.
    In a strict sense, no perfectly flat 'connector' faces exist.
    So we typically return [] for connectors.
    """
    expected_dimensions = common_sphere

    def _compute_volume(self):
        r = self.dimensions['radius']
        return (4.0/3.0)*math.pi*(r**3)

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        return 4.0*math.pi*(r**2)

    def _compute_connectors(self):
        # Spheres have no truly flat face to connect.
        return []


class Hemisphere(Shape3D):
    """
    Hemisphere with radius=r.
    One connector (the flat circular base).
    """
    expected_dimensions = common_sphere

    def _compute_volume(self):
        r = self.dimensions['radius']
        return (2.0/3.0)*math.pi*(r**3)

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        # Curved surface area = 2πr^2
        # The flat cross-section area (open) = πr^2
        # If it's closed, we might add that, but typically "hemisphere" is open.
        # So total "internal" area might be 3πr^2 if we consider the open face.
        return 3.0*math.pi*(r**2)

    def _compute_connectors(self):
        r = self.dimensions['radius']
        return [Connector(area=math.pi*r**2, axis=(0,0,-1), name="hemisphere_flat")]


# For SquarePyramid:
square_pyramid = {
    "side": ["side", "length", "width", "a", "s"],
    "height": ["height", "h"]
}

class SquarePyramid(Shape3D):
    """
    Square-based pyramid with side=a and height=h.
    Has 1 connectable face (square base) with normal -z (assuming apex up).
    """
    expected_dimensions = square_pyramid

    def _compute_volume(self):
        a = self.dimensions['side']
        h = self.dimensions['height']
        return (a**2 * h) / 3.0

    def _compute_surface_area(self):
        a = self.dimensions['side']
        h = self.dimensions['height']
        base_area = a**2
        # Slant height
        slant = math.sqrt((a/2.0)**2 + h**2)
        # Four triangular faces
        lateral_area = a * slant * 2.0  # Because each triangle is (a*slant)/2, times 4 => 2*a*slant
        return base_area + lateral_area

    def _compute_connectors(self):
        a = self.dimensions['side']
        # The base area is a^2
        return [Connector(area=a**2, axis=(0,0,-1), name="pyramid_base")]

# For OpenSquare1 and OpenSquare2:
open_square = {
    "side": ["side", "width", "depth", "w", "d"],
    "height": ["height", "length", "h", "l"]
}


class OpenSquare1(Shape3D):
    """
    A square-based box with ONE open face (like an open-top box).

    Required dimensions:
      side (the length of each side of the square base)
      height
    Volume:
      side^2 * height
    Surface area:
      4 * side * height + (bottom face area)
      = (4 * side * height) + (side^2)
    Connectors:
      One connector at the open face (the top).
        - The bottom is closed, so no connector there.
    """
    expected_dimensions = open_square

    def _compute_volume(self):
        s = self.dimensions['side']
        h = self.dimensions['height']
        return s * s * h

    def _compute_surface_area(self):
        s = self.dimensions['side']
        h = self.dimensions['height']
        # Side walls: 4 * s * h
        # Bottom: s^2
        return (4.0 * s * h) + (s**2)

    def _compute_connectors(self):
        """
        The open face is the top: area = side^2, normal +z
        """
        s = self.dimensions['side']
        top_area = s**2
        return [Connector(area=top_area, axis=(0,0,1), name="open_square1_top")]


class OpenSquare2(Shape3D):
    """
    A square-based box with TWO open faces (no top, no bottom).

    Required dimensions:
      side
      height
    Volume:
      side^2 * height
    Surface area:
      Only the 4 vertical walls => 4 * side * height
    Connectors:
      Two connectors: top (+z) and bottom (-z).
    """
    expected_dimensions = open_square

    def _compute_volume(self):
        s = self.dimensions['side']
        h = self.dimensions['height']
        return s * s * h

    def _compute_surface_area(self):
        s = self.dimensions['side']
        h = self.dimensions['height']
        # No top or bottom => 4 side walls
        return 4.0 * s * h

    def _compute_connectors(self):
        s = self.dimensions['side']
        area_face = s**2
        top_connector = Connector(area=area_face, axis=(0,0,1),  name="open_square2_top")
        bot_connector = Connector(area=area_face, axis=(0,0,-1), name="open_square2_bottom")
        return [top_connector, bot_connector]

# %% Practical aliases for base shapes

# Aliases for cylindrical shapes:
Cyl = Cylinder         # Short alias for Cylinder
OpenCyl = OpenCylinder1  # Alias for OpenCylinder1 (one open end)
Tube = OpenCylinder2    # Alias for OpenCylinder2 (two open ends)

# Aliases for conical shapes:
Con = Cone             # Short alias for Cone
OpenCon = OpenCone     # Alias for OpenCone (base removed)

# Aliases for prismatic shapes:
Box = RectangularPrism  # Alias for RectangularPrism (a box or cuboid)
Cuboid = RectangularPrism  # Alternate alias for RectangularPrism
OpenBox = OpenPrism1      # Alias for OpenPrism1 (one open face)
OpenBox2 = OpenPrism2     # Alias for OpenPrism2 (two open faces)

# Aliases for spherical shapes:
Ball = Sphere         # Alias for Sphere
HemiSphere = Hemisphere  # Common alias for Hemisphere (also HemiSphere)
HalfSphere = Hemisphere  # Another alias for Hemisphere

# Aliases for pyramidal shapes:
Pyramid = SquarePyramid  # Alias for SquarePyramid (a square-based pyramid)
OpenSquareBox = OpenSquare1  # Alias for OpenSquare1 (one open face)
OpenSquareBox2 = OpenSquare2 # Alias for OpenSquare2 (two open faces)


# %% Shape registry and synonyms

# ----------------------------------------------------------------------------
# Shape Registry (synonyms):
# ----------------------------------------------------------------------------

SHAPE_REGISTRY = {
    # Existing geometry classes
    "cylinder": Cylinder,
    "cone": Cone,
    "rectangular_prism": RectangularPrism,
    "sphere": Sphere,
    "square_pyramid": SquarePyramid,
    "hemisphere": Hemisphere,
    "cube": RectangularPrism,   # special case => length=width=height
    "box": RectangularPrism,
    "prism": RectangularPrism,
    "can": Cylinder,
    "pyramid": SquarePyramid,
    "cuboid": RectangularPrism,
    "ball": Sphere,

    # Open geometries
    "open_box": OpenPrism1,
    "openbox": OpenPrism1,
    "open_box2": OpenPrism2,
    "openbox2": OpenPrism2,
    "open_cylinder_1": OpenCylinder1,
    "open_cylinder_2": OpenCylinder2,
    "open_square1": OpenSquare1,
    "open_square2": OpenSquare2,
    "open_prism1": OpenPrism1,
    "open_prism2": OpenPrism2,
    "open_cone": OpenCone,


    # Synonyms for open containers
    "box_container": OpenPrism1,
    "bowl": Hemisphere,
    "halfsphere": Hemisphere,


    # Synonyms for an open cylinder with one open end:
    "glass": OpenCylinder1,
    "pot": OpenCylinder1,
    "jar": OpenCylinder1,

    # Synonym for an open cylinder with two open ends:
    "straw": OpenCylinder2,
    "tube": OpenCylinder2,

}

# ----------------------------------------------------------------------------
# Shape Metadata:
# ----------------------------------------------------------------------------

SHAPE_PARAMETER_SPEC = {
    "Cylinder": {
        "required": ["radius", "height"],
        "doc": (
            "A standard cylinder with top and bottom faces.\n"
            "Volume = π r² h. Surface area includes top and bottom disks."
        ),
    },
    "OpenCylinder1": {
        "required": ["radius", "height"],
        "doc": (
            "A cylinder with exactly one open end (like a glass).\n"
            "Volume = π r² h. Surface area = 2πrh + πr²."
        ),
    },
    "OpenCylinder2": {
        "required": ["radius", "height"],
        "doc": (
            "A cylinder with two open ends (like a straw).\n"
            "Volume = π r² h. Surface area = 2πrh (no top or bottom)."
        ),
    },
    "Cone": {
        "required": ["radius", "height"],
        "doc": (
            "A full cone with closed circular base.\n"
            "Volume = (1/3) π r² h. Surface area = base + lateral area."
        ),
    },
    "OpenCone": {
        "required": ["radius", "height"],
        "doc": (
            "A cone with its base removed, leaving a single open circular face.\n"
            "Volume = (1/3) π r² h. Surface area = π r * slant (no base)."
        ),
    },
    "RectangularPrism": {
        "required": ["length", "width", "height"],
        "doc": (
            "A rectangular prism with all faces closed.\n"
            "Volume = l * w * h. Surface area = 2(lw + lh + wh)."
        ),
    },
    "SquarePyramid": {
        "required": ["side", "height"],
        "doc": (
            "A square-based pyramid.\n"
            "Volume = (side² * height) / 3. Surface area = base + 4 triangles."
        ),
    },
    "Hemisphere": {
        "required": ["radius"],
        "doc": (
            "A hemisphere (half a sphere) typically open at the flat side.\n"
            "Volume = (2/3) π r³. Surface area = 3π r² (2πr² curved + πr² open)."
        ),
    },
    "Sphere": {
        "required": ["radius"],
        "doc": (
            "A full sphere.\n"
            "Volume = (4/3) π r³. Surface area = 4π r²."
        ),
    },
    "OpenSquare1": {
        "required": ["side", "height"],
        "doc": (
            "A square-based box with ONE open face (like an open-top box).\n"
            "Volume = side² * height.\n"
            "Surface area = bottom + 4 walls = side² + 4 side * height."
        ),
    },
    "OpenSquare2": {
        "required": ["side", "height"],
        "doc": (
            "A square-based box with TWO open faces (no top, no bottom).\n"
            "Volume = side² * height.\n"
            "Surface area = 4 side * height."
        ),
    },
    "OpenPrism1": {
        "required": ["length", "width", "height"],
        "doc": (
            "A rectangular prism with ONE open face (e.g. open top).\n"
            "Volume = l * w * h.\n"
            "Surface area = 2(lw + lh + wh) - lw (remove top)."
        ),
    },
    "OpenPrism2": {
        "required": ["length", "width", "height"],
        "doc": (
            "A rectangular prism with TWO open faces (no top, no bottom).\n"
            "Volume = l * w * h.\n"
            "Surface area = 2(lw + lh + wh) - 2(lw)."
        ),
    },
}

# ----------------------------------------------------------------------------
# Helper Functions for Documentation
# ----------------------------------------------------------------------------

def get_geometries_and_synonyms():
    """
    Returns a dictionary mapping each shape class name
    to a sorted list of all registry keys (synonyms) that point to it.

    Example return:
    {
      "Cylinder": ["can", "cylinder"],
      "OpenCylinder1": ["glass", "jar", "open_cylinder_1", "pot"],
      ...
    }
    """
    class_to_names = defaultdict(list)
    for shape_name, shape_cls in SHAPE_REGISTRY.items():
        class_to_names[shape_cls.__name__].append(shape_name)
    result = {}
    for cls_name, synonyms in class_to_names.items():
        result[cls_name] = sorted(synonyms)
    return result

def get_all_shapes_info():
    """
    Returns a dictionary that combines synonyms, required parameters,
    documentation, and expected keywords (derived from the shape class)
    for each shape class.

    Example structure:
    {
      'Cylinder': {
          'synonyms': ['can', 'cylinder'],
          'required_params': ['radius', 'height'],
          'doc': '...',
          'expected_keywords': 'Expected keywords for Cylinder: ...'
      },
      ...
    }
    """
    shape_synonyms_map = get_geometries_and_synonyms()  # {class_name -> [synonyms]}
    # Build a mapping from class name to the shape class object (one instance per class)
    class_to_obj = {}
    for shape_name, shape_cls in SHAPE_REGISTRY.items():
        class_to_obj[shape_cls.__name__] = shape_cls

    all_info = {}
    for cls_name, synonyms in shape_synonyms_map.items():
        param_spec = SHAPE_PARAMETER_SPEC.get(cls_name, {})
        required = param_spec.get("required", [])
        doc_str = param_spec.get("doc", "No documentation available.")
        shape_cls = class_to_obj[cls_name]
        expected_keywords = shape_cls.list_expected_keywords() if hasattr(shape_cls, "list_expected_keywords") else "N/A"
        all_info[cls_name] = {
            "synonyms": synonyms,
            "required_params": required,
            "doc": doc_str,
            "expected_keywords": expected_keywords,
        }
    return all_info


def help_geometry():
    """
    Returns a pretty-formatted string showing all shape classes,
    their synonyms, required parameters, expected keywords (synonym mapping),
    and documentation.

    Example usage:
      help_geometry()
    """
    info = get_all_shapes_info()
    lines = []
    lines.append("=== List of Implemented Geometries & Synonyms ===\n")
    for cls_name in sorted(info.keys()):
        synonyms = info[cls_name]["synonyms"]
        required_params = info[cls_name]["required_params"]
        doc_text = info[cls_name]["doc"]
        expected_keywords = info[cls_name]["expected_keywords"]

        lines.append(f"Shape Class: {cls_name}")
        lines.append(f"  Synonyms       : {', '.join(synonyms)}")
        lines.append(f"  Required Params: {', '.join(required_params) if required_params else 'None'}")
        lines.append("  Expected Keywords:")
        for line in expected_keywords.splitlines():
            lines.append("    " + line)
        lines.append("  Documentation:")
        for line in doc_text.splitlines():
            lines.append("    " + line)
        lines.append("-" * 60)
    prettytxt = "\n".join(lines)
    print(prettytxt)

# ----------------------------------------------------------------------------
# Main Factory Function:
# ----------------------------------------------------------------------------

def create_shape_by_name(name, **dimensions):
    """
    Factory function to create either a single shape or a known composite shape.

    For a direct shape, we find it in SHAPE_REGISTRY.
    For a composite shape (like 'bottle'), we build it from simpler shapes.
    """
    lower_name = name.lower()
    if lower_name == "bottle":
        # Example composite: a bottle with a body and neck.
        body_radius = _to_m(dimensions["body_radius"])
        body_height = _to_m(dimensions["body_height"])
        neck_radius = _to_m(dimensions["neck_radius"])
        neck_height = _to_m(dimensions["neck_height"])
        body = Cylinder(radius=body_radius, height=body_height)
        neck = Cylinder(radius=neck_radius, height=neck_height)
        bottle_composite = CompositeShape()  # Assuming CompositeShape is defined elsewhere.
        bottle_composite.add_shape(body)
        bottle_composite.add_shape(neck, connect_axis=(0, 0, 1))
        return bottle_composite
    else:
        shape_class = SHAPE_REGISTRY.get(lower_name, None)
        if shape_class is None:
            raise ValueError(f"Unknown shape or composite name '{name}'.")
        if lower_name == "cube":
            side = _to_m(dimensions['side'])
            return RectangularPrism(length=side, width=side, height=side)
        return shape_class(**dimensions)


# %% Class for composite shapes
# ----------------------------------------------------------------------------
# Composite Geometry Class (valid approximation for 3D-->1D simulation)
# ----------------------------------------------------------------------------

class CompositeShape(Shape3D):
    """
    Represents a shape made by combining multiple sub-shapes.
    The total volume is the sum of sub-shapes' volumes.

    For surface area, we use the naive sum minus the overlapped face
    (twice the minimum connectable area).

    This version tracks connector usage so that free (external) connectors
    can be computed. This approach is generic: any Shape3D, composite or not,
    exposes its connection faces via the connectors() method.
    """
    def __init__(self):
        # Since a composite shape is built from sub-shapes,
        # no intrinsic dimensions exist.
        super().__init__()
        self.shapes = []
        # Now, we record connections as tuples:
        # (shapeA, shapeB, connector_from_A, connector_from_B, overlap_area)
        self.connections = []

    def add_shape(self, new_shape, connect_axis=None):
        """
        Add a new shape to this composite. If connect_axis is provided,
        we attempt to find a matching connector on 'new_shape' and an
        opposite connector on an existing shape.

        The matching pair's connectors are recorded in self.connections.
        """
        if not self.shapes:
            # First shape: nothing to connect to.
            self.shapes.append(new_shape)
            return

        if connect_axis is None:
            self.shapes.append(new_shape)
            return

        # Step 1: Gather new_shape connectors matching the provided axis.
        new_connectors = [
            c for c in new_shape.connectors()
            if _axes_almost_equal(c.axis, connect_axis)
        ]
        if not new_connectors:
            self.shapes.append(new_shape)
            return

        # Step 2: Determine the opposite axis.
        opposite_axis = tuple([-a for a in connect_axis])

        # Try to find an existing shape with a matching connector.
        for existing in self.shapes:
            existing_connectors = [
                c for c in existing.connectors()
                if _axes_almost_equal(c.axis, opposite_axis)
            ]
            if not existing_connectors:
                continue  # No matching connector on this shape.

            # Use the first matching pair.
            conn_new = new_connectors[0]
            conn_existing = existing_connectors[0]
            overlap_area = _compute_min_overlap(conn_new, conn_existing)
            self.connections.append((new_shape, existing, conn_new, conn_existing, overlap_area))
            self.shapes.append(new_shape)
            return

        # If no pairing was found, add the shape unconnected.
        self.shapes.append(new_shape)

    def _compute_volume(self):
        return sum(s.volume() for s in self.shapes)

    def _compute_surface_area(self):
        """
        Naively sum the surface areas and subtract twice the total overlap.
        """
        total_area = sum(s.surface_area() for s in self.shapes)
        total_overlap = sum(overlap for (_, _, _, _, overlap) in self.connections)
        return total_area - 2.0 * total_overlap

    def _compute_connectors(self):
        """
        Compute the free (external) connectors of the composite.

        Method:
          1. For each sub-shape, retrieve its connectors.
          2. For each connection recorded in self.connections,
             mark both the involved connectors as used.
          3. Return all connectors from sub-shapes that are not marked used.
        """
        # Dictionary mapping each shape to a list of its used connectors.
        used_connectors = {}
        for shapeA, shapeB, connA, connB, _ in self.connections:
            used_connectors.setdefault(shapeA, []).append(connA)
            used_connectors.setdefault(shapeB, []).append(connB)

        free_connectors = []
        for shape in self.shapes:
            for connector in shape.connectors():
                # If this connector was not used, consider it free.
                if shape in used_connectors and connector in used_connectors[shape]:
                    continue
                free_connectors.append(connector)
        return free_connectors

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

def _axes_almost_equal(axis1, axis2, tol=1e-5):
    """
    Check if two unit vectors are nearly the same (or exactly opposite).
    Because connectors are face normals, we consider "matching" to be
    an axis that is within tolerance of the negative direction or the same,
    depending on your design rules.

    In the code above, for matching we do EXACT direction or EXACT opposite.
    Adjust to your preference.
    """
    # We'll check if either they're almost the same or almost exact opposites
    dot = axis1[0]*axis2[0] + axis1[1]*axis2[1] + axis1[2]*axis2[2]
    # If dot ~ 1.0 => same direction, if dot ~ -1.0 => opposite direction
    return abs(abs(dot) - 1.0) < tol

def _compute_min_overlap(connector1, connector2):
    """
    The overlap area is the minimum of the two connectable faces,
    since you can't overlap more than the smaller face area.
    """
    return min(connector1.area, connector2.area)

# %% Packaging3D

# ----------------------------------------------------------------------------
# High Level "Packaging3D" class
# ----------------------------------------------------------------------------

class Packaging3D:
    """
    High-level interface that creates a shape/composite shape by name
    and provides volume & surface area in SI units.

    usage:
      pkg = Packaging3D('bottle',
                        body_radius=(5, 'cm'),
                        body_height=(20, 'cm'),
                        neck_radius=(1.5, 'cm'),
                        neck_height=(5, 'cm'))
      vol, area = pkg.get_volume_and_area()
    """
    def __init__(self, geometry_name, **dimensions):
        self.geometry_name = geometry_name
        self.shape = create_shape_by_name(geometry_name, **dimensions)

    def get_volume_and_area(self):
        """
        Returns: (volume_in_m3, surface_area_in_m2)
        """
        return (self.shape.volume(), self.shape.surface_area())

    def __repr__(self):
        """String representation of Packaging3D, including the nested shape."""
        print(f"Packaging3D(geometry_name='{self.geometry_name}', shape=\n{repr(self.shape)})")
        return str(self)

    def __str__(self):
        """Formatted string representation of the Packaging 3D"""
        return f"<{self.__class__.__name__}: {self.geometry_name}>"


    # --------------------------------------------------------------------
    # For convenience, several operators have been overloaded
    #   packaging >> medium      # sets the volume and the surfacearea
    # --------------------------------------------------------------------

    # method: medium._to(material) and its associated operator >>
    def _to(self,other=None):
        """Propagates volume and area to a food instance"""
        from patankar.food import foodphysics
        if not isinstance(other,foodphysics):
            raise TypeError(f"other must be a foodphysics instance not a {type(other).__name__}")
        other.volume,other.surfacearea = self.get_volume_and_area()
        # we record in other the properties inherited and then transferable
        other.acknowledge(what={"volume","surfacearea"},category="geometry")
        return other

    # the >> and @ operators are depreciated and should not be used anymore
    def __rshift__(self, other): # depreciated
        """
            Overloads >> to enable: packaging >> food --> food
            This use is depreciated, use food << packaging instead
        """
        self._to(other)
        return other

    def __matmul__(self, other): # depreciated
        """
            Overloads @ to enable: packaging@food --> food
            This use is depreciated, use food << packaging instead
        """
        self._to(other)
        return other

# %% Test
# ----------------------------------------------------------------------------
# USAGE EXAMPLES
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    # 1) A "bottle" composed of two cylinders
    bottle_pkg = Packaging3D(
        "bottle",
        body_radius=(50, "mm"), # 0.05 m
        body_height=(0.2, "m"), # 0.20 m
        neck_radius=(2, "cm"),  # 0.02 m
        neck_height=0.05        # 0.05 m
    )
    b_vol, b_area = bottle_pkg.get_volume_and_area()
    print("Bottle Volume (m**3):", b_vol)
    print("Bottle Surface (m**2):", b_area)

    # 2) A single shape: "can" is just a cylinder
    can_pkg = Packaging3D("can", radius=(4,"cm"), height=(12,"cm"))
    c_vol, c_area = can_pkg.get_volume_and_area()
    print("Can Volume (m**3):", c_vol)
    print("Can Surface (m**2):", c_area)

    # 3) A "cube" with side=10 cm
    cube_pkg = Packaging3D("cube", side=(10,"cm"))
    cu_vol, cu_area = cube_pkg.get_volume_and_area()
    print("Cube Volume (m**3):", cu_vol)
    print("Cube Surface (m**2):", cu_area)
