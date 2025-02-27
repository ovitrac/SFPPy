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
@Date: 2024-10-28, rev. 2025-02-22

===============================================================================
"""


# %% Dependencies
import math
import numpy as np
from collections import defaultdict
from patankar.layer import check_units

__all__ = ['CompositeShape', 'Cone', 'Connector', 'Cylinder', 'Hemisphere', 'OpenCone', 'OpenCylinder1', 'OpenCylinder2', 'OpenPrism1', 'OpenPrism2', 'OpenSquare1', 'OpenSquare2', 'Packaging3D', 'RectangularPrism', 'Shape3D', 'Sphere', 'SquarePyramid', 'check_units', 'create_shape_by_name', 'get_all_shapes_info', 'get_geometries_and_synonyms', 'help_geometry']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.2"
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

# %% Private Classes
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
    Base class for a 3D shape. Subclasses must implement:
      - _compute_volume()
      - _compute_surface_area()
      - _compute_connectors() -> list of Connector objects
    """
    def __init__(self, **dimensions):
        # Convert every dimension to meters
        self.dimensions = {k: _to_m(v) for k, v in dimensions.items()}

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
        """
        Return a list of Connector objects that represent the shape’s
        possible connections. By default, shapes with no flat faces return [].
        """
        return []

    def __repr__(self):
        """String representation of the Shape3D object."""
        class_name = self.__class__.__name__

        # Convert numpy arrays to scalars before formatting
        dimensions_str = ", ".join(f"{k}={v.item():.4g} m" if isinstance(v, np.ndarray) else f"{k}={v:.4f} m"
                                   for k, v in self.dimensions.items())

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
        """Formatted string representing the 3D shape"""
        n = len(self.connectors())
        return f"<{self.__class__.__name__} with {n} connector{'s' if n>1 else ''}>"


# ----------------------------------------------------------------------------
# Basic shapes
# ----------------------------------------------------------------------------

class Cylinder(Shape3D):
    """
    A cylinder with radius=r and height=h.
    Has two connectors (top and bottom).
    """
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

class RectangularPrism(Shape3D):
    """
    A rectangular prism with length=l, width=w, height=h.
    Has 6 connectors for each face.
    """
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

class Sphere(Shape3D):
    """
    A sphere with radius=r.
    In a strict sense, no perfectly flat 'connector' faces exist.
    So we typically return [] for connectors.
    """
    def _compute_volume(self):
        r = self.dimensions['radius']
        return (4.0/3.0)*math.pi*(r**3)

    def _compute_surface_area(self):
        r = self.dimensions['radius']
        return 4.0*math.pi*(r**2)

    def _compute_connectors(self):
        # Spheres have no truly flat face to connect.
        return []

class SquarePyramid(Shape3D):
    """
    Square-based pyramid with side=a and height=h.
    Has 1 connectable face (square base) with normal -z (assuming apex up).
    """
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

class Hemisphere(Shape3D):
    """
    Hemisphere with radius=r.
    One connector (the flat circular base).
    """
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


# %% Main code

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
    "bowl": Hemisphere,

    # Open geometries
    "open_cylinder_1": OpenCylinder1,
    "open_cylinder_2": OpenCylinder2,
    "open_square1": OpenSquare1,
    "open_square2": OpenSquare2,
    "open_prism1": OpenPrism1,
    "open_prism2": OpenPrism2,
    "open_cone": OpenCone,

    # Synonyms for open containers
    "box_container": OpenPrism1,

    # Synonyms for an open cylinder with one open end:
    "glass": OpenCylinder1,
    "pot": OpenCylinder1,
    "jar": OpenCylinder1,

    # Synonym for an open cylinder with two open ends:
    "straw": OpenCylinder2,
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
# autodoc function:
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
        # e.g. shape_cls.__name__ => 'Cylinder'
        class_to_names[shape_cls.__name__].append(shape_name)

    # Sort synonyms for a consistent presentation
    result = {}
    for cls_name, synonyms in class_to_names.items():
        result[cls_name] = sorted(synonyms)
    return result

def get_all_shapes_info():
    """
    Returns a dictionary that combines synonyms, required parameters,
    and doc strings for each shape class.

    Example structure:
    {
      'Cylinder': {
          'synonyms': ['can', 'cylinder'],
          'required_params': ['radius', 'height'],
          'doc': '...'
      },
      'OpenCylinder1': {
          'synonyms': ['glass', 'jar', 'open_cylinder_1', 'pot'],
          'required_params': ['radius', 'height'],
          'doc': '...'
      },
      ...
    }
    """
    shape_synonyms_map = get_geometries_and_synonyms()  # {class_name -> [synonyms]}
    all_info = {}
    for cls_name, synonyms in shape_synonyms_map.items():
        param_spec = SHAPE_PARAMETER_SPEC.get(cls_name, {})
        required = param_spec.get("required", [])
        doc_str = param_spec.get("doc", "No documentation available.")

        all_info[cls_name] = {
            "synonyms": synonyms,
            "required_params": required,
            "doc": doc_str,
        }
    return all_info

def help_geometry():
    """
    Returns a pretty-formatted string showing all shape classes,
    their synonyms, required parameters, and documentation.

    Example usage:
      help_geometry()
    """
    info = get_all_shapes_info()

    lines = []
    lines.append("=== List of Implemented Geometries & Synonyms ===\n")

    # Sort by class name for consistency
    for cls_name in sorted(info.keys()):
        synonyms = info[cls_name]["synonyms"]
        required_params = info[cls_name]["required_params"]
        doc_text = info[cls_name]["doc"]

        lines.append(f"Shape Class: {cls_name}")
        lines.append(f"  Synonyms       : {', '.join(synonyms)}")
        lines.append(f"  Required Params: {', '.join(required_params) if required_params else 'None'}")

        # Optionally indent doc lines
        doc_lines = doc_text.split("\n")
        for dl in doc_lines:
            lines.append(f"    {dl}")
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

    # Example of a special composite shape: 'bottle'
    # A "bottle" can be modeled as: a large cylinder (body) + smaller cylinder (neck)
    # joined along their circular faces.
    if lower_name == "bottle":
        # We expect: body_radius, body_height, neck_radius, neck_height
        body_radius = _to_m(dimensions["body_radius"])
        body_height = _to_m(dimensions["body_height"])
        neck_radius = _to_m(dimensions["neck_radius"])
        neck_height = _to_m(dimensions["neck_height"])

        # Create the big cylinder
        body = Cylinder(radius=body_radius, height=body_height)
        # Create the smaller cylinder for neck
        neck = Cylinder(radius=neck_radius, height=neck_height)

        # Combine them
        bottle_composite = CompositeShape()
        bottle_composite.add_shape(body)
        bottle_composite.add_shape(neck, connect_axis=(0,0,1))

        return bottle_composite
    else:
        # If it's a direct geometry name or a synonym:
        shape_class = SHAPE_REGISTRY.get(lower_name, None)
        if shape_class is None:
            raise ValueError(f"Unknown shape or composite name '{name}'.")

        # Special case for "cube": user might supply side=...
        # or length=... In normal usage for a "cube" we do side=...
        # We'll unify as rectangular prism with l=w=h=side
        if lower_name == "cube":
            # We expect 'side' => l=w=h
            side = _to_m(dimensions['side'])
            return RectangularPrism(length=side, width=side, height=side)

        # Otherwise just create the shape with the given dimensions
        return shape_class(**dimensions)

# ----------------------------------------------------------------------------
# Composite Geometry Class (valid approximation for 3D-->1D simulation)
# ----------------------------------------------------------------------------

class CompositeShape(Shape3D):
    """
    Represents a shape made by combining multiple sub-shapes.
    The total volume is the sum of sub-shapes' volumes.

    For surface area, we use the naive sum minus the overlapped face
    (twice the minimum connectable area).

    In a real system you might track the exact arrangement in 3D space,
    but here we keep it conceptual for demonstration:
      - add_shape(shape, connect_axis): we try to connect the new shape along
        a matching connector from an existing shape if axes align.
    """
    def __init__(self):
        super().__init__()  # empty base
        self.shapes = []
        self.connections = []  # List of (shapeA, shapeB, overlap_area)

    def add_shape(self, new_shape, connect_axis=None):
        """
        Add a new shape to this composite. If connect_axis is provided,
        we attempt to find a connector on 'new_shape' that matches
        a connector on an existing shape in self.shapes.

        For demonstration, we connect to the first available match.
        Overlap area is min( area1, area2 ).
        """
        if not self.shapes:
            # If this is the first shape in the composite, just add it.
            self.shapes.append(new_shape)
            return

        if connect_axis is None:
            # No connector logic needed, just add it unconnected
            self.shapes.append(new_shape)
            return

        # We'll search for a matching connector (by orientation) in the new_shape,
        # and see if we can pair it with an existing shape's connector.
        # If found, record the overlap.

        # Step 1: gather new_shape connectors that match the connect_axis
        new_connectors = [
            c for c in new_shape.connectors()
            if _axes_almost_equal(c.axis, connect_axis)
        ]
        if not new_connectors:
            # If we found no matching connectors, we just add shape
            self.shapes.append(new_shape)
            return

        # Step 2: gather existing shapes' connectors that match the opposite axis
        # i.e. connect_axis is (0,0,1), we might need existing axis to be (0,0,-1)
        # so that they can connect “face-to-face”.
        opposite_axis = tuple([-a for a in connect_axis])

        for existing in self.shapes:
            existing_connectors = [
                c for c in existing.connectors()
                if _axes_almost_equal(c.axis, opposite_axis)
            ]
            if not existing_connectors:
                continue  # no suitable connector on that shape

            # We’ll just connect the first pair of connectors we find
            # (In a real system, you’d define a better approach or prompt the user.)
            overlap_area = _compute_min_overlap(new_connectors[0], existing_connectors[0])
            self.connections.append((new_shape, existing, overlap_area))

            # Add the new shape to the composite
            self.shapes.append(new_shape)
            return

        # If we get here, no suitable pairing was found. Just add shape unconnected
        self.shapes.append(new_shape)

    def _compute_volume(self):
        return sum(s.volume() for s in self.shapes)

    def _compute_surface_area(self):
        """
        Sum of the sub-shapes’ surface areas minus
        2 * sum of each overlapping face area.
        """
        total_area = sum(s.surface_area() for s in self.shapes)

        overlap = 0.0
        for (shapeA, shapeB, overlap_area) in self.connections:
            overlap += overlap_area

        # We remove twice the overlap area for each connection
        # (once from shape A, once from shape B).
        return total_area - 2.0 * overlap

    def _compute_connectors(self):
        """
        As a composite, its external connectors might be complicated.
        Here we return an empty list, or you could gather connectors
        that are not overlapped.
        """
        return []

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

    def __rshift__(self, other):
        """Overloads >> to propagate to other."""
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
