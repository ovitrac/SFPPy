#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Food Layer
===============================================================================
Defines **food materials** for migration simulations. Models food as a **0D layer** with:
- **Mass transfer resistance (`h`)**
- **Partitioning (`k0`)**
- **Contact time & temperature**

**Main Components:**
- **Base Class: `foodphysics`** (Stores all food-related parameters)
    - Defines mass transfer properties (`h`, `k0`)
    - Implements property propagation (`food >> layer`)
- **Subclasses:**
    - `foodlayer`: General food layer model
    - `setoff`: Periodic boundary conditions (e.g., stacked packaging)
    - `nofood`: Impervious boundary (no mass transfer)
    - `realcontact` & `testcontact`: Standardized storage and testing conditions

**Integration with SFPPy Modules:**
- Works with `migration.py` as the **left-side boundary** for simulations.
- Can inherit properties from `layer.py` for **contact temperature propagation**.
- Used in `geometry.py` when defining food-contacting packaging.

Example:
```python
from patankar.food import foodlayer
medium = foodlayer(name="ethanol", contacttemperature=(40, "degC"))
```


@version: 1.2
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2023-01-25
@rev: 2025-02-14

"""

# Dependencies
import sys
import inspect
import textwrap
import numpy as np
from copy import deepcopy as duplicate

from patankar.layer import check_units, NoUnits, layer # to convert units to SI

__all__ = ['ambient', 'aqueous', 'boiling', 'check_units', 'chemicalaffinity', 'chilled', 'ethanol', 'ethanol50', 'fat', 'foodlayer', 'foodphysics', 'foodproperty', 'frozen', 'get_defined_init_params', 'help_food', 'hotfilled', 'intermediate', 'is_valid_classname', 'layer', 'liquid', 'list_food_classes', 'nofood', 'oven', 'pasteurization', 'perfectlymixed', 'realcontact', 'realfood', 'rolled', 'semisolid', 'setoff', 'simulant', 'solid', 'stacked', 'sterilization', 'tenax', 'testcontact', 'texture', 'water', 'wrap_text', 'yogurt']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.2"
#%% Private Properties and functions

# List of the default SI units used by physical quantity
parametersWithUnits = {"volume":"m**3",
                       "surfacearea":"m**2",
                       "density":"kg/m**3",
                       "contacttemperature":"degC",
                       "h":"m/s",
                       "k0":NoUnits,
                       "CF0":NoUnits,
                       "contacttime":"s"
                       }
# corresponding protperty names                       }
paramaterNamesWithUnits = [p+"Units" for p in parametersWithUnits.keys()]

# List parameters not used with nofood, noPBC
parametersWithUnits_andfallback = [key for key in parametersWithUnits if key != "contacttime"]

LEVEL_ORDER = {"base": 0, "root": 1, "property":2, "contact":3, "user": 4}  # Priority order for sorting

def wrap_text(text, width=20):
    """Wraps text within a specified width and returns a list of wrapped lines."""
    if not isinstance(text, str):
        return [str(text)]
    return textwrap.wrap(text, width) or [""]  # Ensure at least one line

def get_defined_init_params(instance):
    """Returns which parameters from parametersWithUnits are defined in the instance."""
    return [param for param in parametersWithUnits.keys() if hasattr(instance, param)]

def is_valid_classname(name):
    """Returns True if class name is valid (not private/internal)."""
    return name.isidentifier() and not name.startswith("_")  # Exclude _10, __, etc.

def list_food_classes():
    """
    Lists all classes in the 'food' module with:
    - name and description
    - level (class attribute)
    - Inheritance details
    - Parameters from parametersWithUnits that are set in the instance
    """
    subclasses_info = []
    current_module = sys.modules[__name__]  # Reference to the food module

    for name, obj in inspect.getmembers(current_module, inspect.isclass):
        if obj.__module__ == current_module.__name__ and is_valid_classname(name):  # Ensure valid class name
            try:
                instance = obj()  # Try to instantiate
                init_params = get_defined_init_params(instance)
                level = getattr(obj, "level", "other")  # Default to "other" if no level is set

                class_info = {
                    "Class Name": wrap_text(name),
                    "Name": wrap_text(getattr(instance, "name", "N/A")),
                    "Description": wrap_text(getattr(instance, "description", "N/A")),
                    "Level": wrap_text(level),
                    "Inheritance": wrap_text(", ".join(base.__name__ for base in obj.__bases__)),
                    "Init Params": wrap_text(", ".join(init_params) if init_params else ""),
                    "Level Sorting": LEVEL_ORDER.get(level, 3)  # Used for sorting, not for table output
                }
                subclasses_info.append(class_info)
            except TypeError:
                class_info = {
                    "Class Name": wrap_text(name),
                    "Name": ["N/A"],
                    "Description": ["N/A"],
                    "Level": wrap_text(getattr(obj, "level", "other")),
                    "Inheritance": wrap_text(", ".join(base.__name__ for base in obj.__bases__)),
                    "Init Params": wrap_text("⚠️ Cannot instantiate"),
                    "Level Sorting": LEVEL_ORDER.get(getattr(obj, "level", "other"), 3)
                }
                subclasses_info.append(class_info)

    # **Sort first by level priority, then alphabetically within each level**
    subclasses_info.sort(key=lambda x: (x["Level Sorting"], x["Class Name"]))

    return subclasses_info

def help_food():
    """
    Prints all food-related classes with relevant attributes in a **formatted Markdown table**.
    """
    derived = list_food_classes()

    # Define table headers (excluding "Level Sorting" because it's only used for sorting)
    headers = ["Class Name", "Name", "Description", "Level", "Inheritance", "Init Params"]

    # Find the maximum number of lines in any wrapped column (excluding "Level Sorting")
    max_lines_per_row = [
        max(len(value) for key, value in row.items() if key != "Level Sorting")
        for row in derived
    ]

    # Convert dictionary entries to lists and ensure they all have the same number of lines
    formatted_rows = []
    for row, max_lines in zip(derived, max_lines_per_row):
        wrapped_row = {
            key: (value if isinstance(value, list) else [value]) + [""] * (max_lines - len(value))
            for key, value in row.items() if key != "Level Sorting"  # Exclude "Level Sorting"
        }
        for i in range(max_lines):  # Transpose wrapped lines into multiple rows
            formatted_rows.append([wrapped_row[key][i] for key in headers])

    # Compute column widths dynamically
    col_widths = [max(len(str(cell)) for cell in col) for col in zip(headers, *formatted_rows)]

    # Create a formatting row template
    row_format = "| " + " | ".join(f"{{:<{w}}}" for w in col_widths) + " |"

    # Print the table header
    print(row_format.format(*headers))
    print("|-" + "-|-".join("-" * w for w in col_widths) + "-|")

    # Print all table rows
    for row in formatted_rows:
        print(row_format.format(*row))


#%% Base physics class
# -------------------------------------------------------------------
# Base Class to convert class defaults to instance attributes
# -------------------------------------------------------------------
class foodphysics:
    """
    Base class that automatically assigns instance attributes from class defaults,
    except for the 'description' attribute.
    Check the physical meaning of quantities with units.

    Implemented methods include:
        - refresh() validates all quantities before a simulation
        - update(name="new name", description="new description",parameter1=value)
          assigns new values to physical parameters and attributes
        - getparam() returns physical parameters even if they undefined
    Available properties:
        PBC returns True in periodic boundary conditions are enforced (setoff)
        impervious returns True if impervious boundary condition is appled (no food)
    """

    # General descriptors
    description = "Root physics class used to implement food and mass transfer physics"  # Remains as class attribute
    name = "food physics"
    level = "base"

    # ------------------------------------------------------
    # Transfer rules for food1 >> food2 and food1 >> result
    # ------------------------------------------------------

    # Mapping of properties to their respective categories
    _list_categories = {
        "contacttemperature": "contact",
        "contacttime": "contact",
        "surfacearea": "geometry",
        "volume": "geometry"
    }

    # Rules for property transfer based on object type
    _transferable_properties = {
        "contacttemperature": {
            "foodphysics": {
                "onlyifinherited": True,
                "checkNumPy": False,
                "as": "",
                "prototype": None,
                "category": "contact"
            },
            "layer": {
                "onlyifinherited": False,
                "checkNumPy": True,
                "as": "T",
                "prototype": None
            }
        },
        "contacttime": {
            "foodphysics": {
                "onlyifinherited": True,
                "checkNumPy": True,
                "as": "",
                "prototype": None
            },
            "SensPatankarResult": {
                "onlyifinherited": False,
                "checkNumPy": True,
                "as": "t",
                "prototype": None
            }
        },
        "surfacearea": {
            "foodphysics": {
                "onlyifinherited": False,
                "checkNumPy": False,
                "as": "surfacearea",
                "prototype": None
            }
        },
        "volume": {
            "foodphysics": {
                "onlyifinherited": False,
                "checkNumPy": True,
                "as": "",
                "prototype": None
            }
        }
    }


    def __init__(self, **kwargs):
        """general constructor"""

        # local import
        from patankar.migration import SensPatankarResult

        # numeric validator
        def numvalidator(key,value):
            if key in parametersWithUnits:          # the parameter is a physical quantity
                if isinstance(value,tuple):         # the supplied value as unit
                    value,_ = check_units(value)    # we convert to SI, we drop the units
                if not isinstance(value,np.ndarray):
                    value = np.array([value])       # we force NumPy class
            return value

        # Iterate through the MRO (excluding foodphysics and object)
        for cls in reversed(self.__class__.__mro__):
            if cls in (foodphysics, object):
                continue
            # For each attribute defined at the class level,
            # if it is not 'description', not callable, and not a dunder, set it as an instance attribute.
            for key, value in cls.__dict__.items(): # we loop on class attributes
                if key in ("description","level") or key.startswith("__") or callable(value):
                    continue
                if key not in kwargs:
                    setattr(self, key, numvalidator(key,value))
        # Now update/override with any keyword arguments provided at instantiation.
        for key, value in kwargs.items():
            value = numvalidator(key,value)
            if key not in paramaterNamesWithUnits: # we protect the values of units (they are SI, they cannot be changed)
                setattr(self, key, value)
        # we initialize the acknowlegment process for future property propagation
        self._hasbeeninherited = {}
        # we initialize the _simstate storing the last simulation result available
        self._simstate = None # simulation results
        self._inpstate = None # their inputs
        # For cooperative multiple inheritance, call the next __init__ if it exists.
        super().__init__()
        # Define actual class references to avoid circular dependency issues
        if self.__class__._transferable_properties["contacttemperature"]["foodphysics"]["prototype"] is None:
            self.__class__._transferable_properties["contacttemperature"]["foodphysics"]["prototype"] = foodphysics
            self.__class__._transferable_properties["contacttemperature"]["layer"]["prototype"] = layer
            self.__class__._transferable_properties["contacttime"]["foodphysics"]["prototype"] = foodphysics
            self.__class__._transferable_properties["contacttime"]["SensPatankarResult"]["prototype"] = SensPatankarResult
            self.__class__._transferable_properties["surfacearea"]["foodphysics"]["prototype"] = foodphysics
            self.__class__._transferable_properties["volume"]["foodphysics"]["prototype"] = foodphysics

    # ------- [properties to access/modify simstate] --------
    @property
    def lastinput(self):
        """Getter for last layer input."""
        return self._inpstate

    @lastinput.setter
    def lastinput(self, value):
        """Setter for last layer input."""
        self._inpstate = value

    @property
    def lastsimulation(self):
        """Getter for last simulation results."""
        return self._simstate

    @lastsimulation.setter
    def lastsimulation(self, value):
        """Setter for last simulation results."""
        self._simstate = value

    @property
    def hassimulation(self):
        """Returns True if a simulation exists"""
        return self.lastsimulation is not None


    # ------- [inheritance registration mechanism] --------
    def acknowledge(self, what=None, category=None):
        """
        Register inherited properties under a given category.

        Parameters:
        -----------
        what : str or list of str or a set
            The properties or attributes that have been inherited.
        category : str
            The category under which the properties are grouped.

        Example:
        --------
        >>> b = B()
        >>> b.acknowledge(what="volume", category="geometry")
        >>> b.acknowledge(what=["surfacearea", "diameter"], category="geometry")
        >>> print(b._hasbeeninherited)
        {'geometry': {'volume', 'surfacearea', 'diameter'}}
        """
        if category is None or what is None:
            raise ValueError("Both 'what' and 'category' must be provided.")
        if isinstance(what, str):
            what = {what}  # Convert string to a set
        elif isinstance(what, list):
            what = set(what)  # Convert list to a set for uniqueness
        elif not isinstance(what,set):
            raise TypeError("'what' must be a string, a list, or a set of strings.")
        if category not in self._hasbeeninherited:
            self._hasbeeninherited[category] = set()
        self._hasbeeninherited[category].update(what)


    def refresh(self):
        """refresh all physcal paramaters after instantiation"""
        for key, value in self.__dict__.items():    # we loop on instance attributes
            if key in parametersWithUnits:          # the parameter is a physical quantity
                if isinstance(value,tuple):         # the supplied value as unit
                    value = check_units(value)[0]   # we convert to SI, we drop the units
                    setattr(self,key,value)
                if not isinstance(value,np.ndarray):
                    value = np.array([value])      # we force NumPy class
                    setattr(self,key,value)

    def update(self, **kwargs):
        """
        Update modifiable parameters of the foodphysics object.

        Modifiable Parameters:
            - name (str): New name for the object.
            - description (str): New description.
            - volume (float or tuple): Volume (can be tuple like (1, "L")).
            - surfacearea (float or tuple): Surface area (can be tuple like (1, "cm^2")).
            - density (float or tuple): Density (can be tuple like (1000, "kg/m^3")).
            - CF0 (float or tuple): Initial concentration in the food.
            - contacttime (float or tuple): Contact time (can be tuple like (1, "h")).
            - contacttemperature (float or tuple): Temperature (can be tuple like (25, "degC")).
            - h (float or tuple): Mass transfer coefficient (can be tuple like (1e-6,"m/s")).
            - k0 (float or tuple): Henri-like coefficient for the food (can be tuple like (1,"a.u.")).

        """
        if not kwargs:  # shortcut
            return self # for chaining
        def checkunits(value):
            """Helper function to convert physical quantities to SI."""
            if isinstance(value, tuple) and len(value) == 2:
                scale = check_units(value)[0]  # Convert to SI, drop unit
                return np.array([scale], dtype=float)  # Ensure NumPy array
            elif isinstance(value, (int, float, np.ndarray)):
                return np.array([value], dtype=float)  # Ensure NumPy array
            else:
                raise ValueError(f"Invalid value for physical quantity: {value}")
        # Update `name` and `description` if provided
        if "name" in kwargs:
            self.name = str(kwargs["name"])
        if "description" in kwargs:
            self.description = str(kwargs["description"])
        # Update physical properties
        for key in parametersWithUnits.keys():
            if key in kwargs:
                value = kwargs[key]
                setattr(self, key, checkunits(value))  # Ensure NumPy array in SI
        return self  # Return self for method chaining if needed

    def get_param(self, key, default=None, acceptNone=True):
        """Retrieve instance attribute with a default fallback if enabled."""
        paramdefaultvalue = 1
        if isinstance(self,(setoff,nofood)):
            if key in parametersWithUnits_andfallback:
                value =  self.__dict__.get(key, paramdefaultvalue) if default is None else self.__dict__.get(key, default)
                if isinstance(value,np.ndarray):
                    value = value.item()
                if value is None and not acceptNone:
                    value = paramdefaultvalue if default is None else default
                return np.array([value])
            if key in paramaterNamesWithUnits:
                return self.__dict__.get(key, parametersWithUnits[key]) if default is None else self.__dict__.get(key, default)
        if key in parametersWithUnits:
            if hasattr(self, key):
                return getattr(self,key)
            else:
                raise KeyError(
                    f"Missing property: '{key}' in instance of class '{self.__class__.__name__}'.\n"
                    f"To define it, use one of the following methods:\n"
                    f"  - Direct assignment:   object.{key} = value\n"
                    f"  - Using update method: object.update({key}=value)\n"
                    f"Note: The value can also be provided as a tuple (value, 'unit')."
                )
        elif key in paramaterNamesWithUnits:
            return self.__dict__.get(key, paramaterNamesWithUnits[key]) if default is None else self.__dict__.get(key, default)
        raise KeyError(f'Use getattr("{key}") to retrieve the value of {key}')

    def __repr__(self):
        """Formatted string representation of the foodphysics object."""
        # refresh all definitions
        self.refresh()
        # Header with name and description
        repr_str = f'Food object "{self.name}" ({self.description}) with properties:\n'

        # Helper function to extract a numerical value safely
        def format_value(value):
            """Ensure the value is a float or a single-item NumPy array."""
            if isinstance(value, np.ndarray):
                return value.item() if value.size == 1 else value[0]  # Ensure scalar representation
            elif value is None:
                return value
            return float(value)
        # Loop through parameters that should be printed
        for key, unit in parametersWithUnits.items():
            if hasattr(self, key):  # Print only defined parameters
                value = format_value(getattr(self, key))
                unit_str = self.get_param(key+"Units", parametersWithUnits[key])  # Retrieve unit safely
                if value is not None:
                    repr_str += f"{key:15s}: {value:0.8g} [{unit_str}]\n"
        print(repr_str.strip())  # Remove trailing newline
        return str(self)


    def __str__(self):
        """Formatted string representation of the property"""
        simstr = ' [simulated]' if self.hassimulation else ""
        return f"<{self.__class__.__name__}: {self.name}>{simstr}"

    def copy(self,**kwargs):
        """Creates a deep copy of the current food instance."""
        return duplicate(self).update(**kwargs)


    @property
    def PBC(self):
        """
        Returns true if h is not defined or None
            This property is used to identified periodic boundary condition also called setoff mass transfer.

        """
        if not hasattr(self,"h"):
            return None
        htmp = getattr(self,"h")
        if isinstance(htmp,np.ndarray):
            htmp = htmp.item()
        return htmp is None


    # --------------------------------------------------------------------
    # For convenience, several operators have been overloaded
    #   medium >> packaging      # sets the volume and the surfacearea
    #   medium >> material       # propgates the contact temperature from the medium to the material
    #   sol = medium << material # simulate migration from the material to the medium
    # --------------------------------------------------------------------

    # method: medium._to(material) and its associated operator >>
    def _to(self, other = None):
        """
        Transfers inherited properties to another object based on predefined rules.

        Parameters:
        -----------
        other : object
            The recipient object that will receive the transferred properties.

        Notes:
        ------
        - Only properties listed in `_transferable_properties` are transferred.
        - A property can only be transferred if `other` matches the expected class.
        - The property may have a different name in `other` as defined in `as`.
        - If `onlyifinherited` is True, the property must have been inherited by `self`.
        - If `checkNumPy` is True, ensures NumPy array compatibility.
        - Updates `other`'s `_hasbeeninherited` tracking.
        """
        for prop, classes in self._transferable_properties.items():
            if prop not in self._list_categories:
                continue  # Skip properties not categorized

            category = self._list_categories[prop]

            for class_name, rules in classes.items():

                if not isinstance(other, rules["prototype"]):
                    continue  # Skip if other is not an instance of the expected prototype class

                if rules["onlyifinherited"] and category not in self._hasbeeninherited:
                    continue  # Skip if property must be inherited but is not

                if rules["onlyifinherited"] and prop not in self._hasbeeninherited[category]:
                    continue  # Skip if the specific property has not been inherited

                if not hasattr(self, prop):
                    continue  # Skip if the property does not exist on self

                # Determine the target attribute name in other
                target_attr = rules["as"] if rules["as"] else prop

                # Retrieve the property value
                value = getattr(self, prop)

                # Handle NumPy array check
                if rules["checkNumPy"] and hasattr(other, target_attr):
                    existing_value = getattr(other, target_attr)
                    if isinstance(existing_value, np.ndarray):
                        value = np.full(existing_value.shape, value)

                # Assign the value to other
                setattr(other, target_attr, value)

                # Register the transfer in other’s inheritance tracking
                other.acknowledge(what=target_attr, category=category)

                # to chain >>
                return other

    def __rshift__(self, other):
        """Overloads >> to propagate to other."""
        return self._to(other)

    # migration method
    def migration(self,material,**kwargs):
        from patankar.migration import senspatankar
        self._to(material) # propagate contact conditions first
        return senspatankar(material,self,**kwargs)

    def contact(self,material,**kwargs):
        return self.migration(self,material,**kwargs)


# %% Root classes
# -------------------------------------------------------------------
# ROOT CLASSES
#   - The foodlayer class represents physically the food
#   - The chemicalaffinity class represents the polarity of the medium (with respect to the substance)
#   - The texture class represents the mass transfer reistance between the food and the material in contact
#   - The nofood class enforces an impervious boundary condition on the food side preventing any transfer.
#     This class is useful to simulate mass transfer within the packaging layer in the absence of food.
#   - The setoff class enforces periodic conditions such as when packaging are stacked together.
# -------------------------------------------------------------------

class foodlayer(foodphysics):
    """
    Foodlayer class is a generic class to define food as a 1D layer in a symmetric manner with layer class
    applicable to materials in contact.
    Since mass transfer are much faster in the food than in the materials in contact, food is represented
    as an almost 0D layer. Only a mass transfer resistance is applied at the food-material interface
    controlled by the mass transfer coefficient h. A Henri-like coefficient k0 controls the eventual
    partitioning of the substance between the food and the layer of the materials.

    Food are geometrically defined by their volume and surface area in contact with the material.

    Contact time (contacttime) and contact temperature (contacttemperature) are defined via foodlayer.

    """
    level = "root"
    description = "root food class"  # Remains as class attribute
    name = "generic food layer"
    volume,volumeUnits = check_units((1,"dm**3"))
    surfacearea,surfaceareaUnits = check_units((6,"dm**2"))
    density,densityUnits = check_units((1000,"kg/m**3"))
    CF0,CF0units = check_units((0,NoUnits))  # initial concentration (arbitrary units)
    contacttime, contacttime_units = check_units((10,"days"))
    contactemperature,contactemperatureUnits = check_units((40,"degC"),ExpectedUnits="degC") # temperature ALWAYS in °C


class texture(foodphysics):
    """Parent food texture class"""
    description = "default class texture"
    name = "undefined"
    level = "root"
    h = 1e-3

class chemicalaffinity(foodphysics):
    """Parent chemical affinity class"""
    description = "default chemical affinity"
    name = "undefined"
    level = "root"
    k0 = 1

class nofood(foodphysics):
    """Impervious boundary condition"""
    description = "impervious boundary condition"
    name = "undefined"
    level = "root"
    h = 0

class setoff(foodphysics):
    """periodic boundary conditions"""
    description = "periodic boundary conditions"
    name = "setoff"
    level = "root"
    h = None

class realcontact(foodphysics):
    """real contact conditions"""
    description = "real storage conditions"
    name = "contact conditions"
    level = "root"
    [contacttime,contacttimeUnits] = check_units((200,"days"))
    [contacttemperature,contacttemperatureUnits] = check_units((25,"degC"))

class testcontact(foodphysics):
    """conditions of migration testing"""
    description = "migration testing conditions"
    name = "migration testing"
    level = "root"
    [contacttime,contacttimeUnits] = check_units((10,"days"))
    [contacttemperature,contacttemperatureUnits] = check_units((40,"degC"))

# %% Property classes
# -------------------------------------------------------------------
# SECOND LEVEL CLASSES
# This classes are used as keyword to define new food with a combination of properties.
# -------------------------------------------------------------------

# Food/chemical properties
class foodproperty(foodlayer):
    """Class wrapper of food properties"""
    level="property"

class realfood(foodproperty):
    """Core real food class (second level)"""
    description = "real food class"

class simulant(foodproperty):
    """Core food simulant class (second level)"""
    name = "generic food simulant"
    description = "food simulant"

class solid(foodproperty):
    """Solid food texture"""
    name = "solid food"
    description = "solid food products"
    [h,hUnits] = check_units((1e-8,"m/s"))

class semisolid(texture):
    """Semi-solid food texture"""
    name = "solid food"
    description = "solid food products"
    [h,hUnits] = check_units((1e-7,"m/s"))

class liquid(texture):
    """Liquid food texture"""
    name = "liquid food"
    description = "liquid food products"
    [h,hUnits] = check_units((1e-6,"m/s"))

class perfectlymixed(texture):
    """Perfectly mixed liquid (texture)"""
    name = "perfectly mixed liquid"
    description = "maximize mixing, minimize the mass transfer boundary layer"
    [h,hUnits] = check_units((1e-4,"m/s"))

class fat(chemicalaffinity):
    """Fat contact"""
    name = "fat contact"
    description = "maximize mass transfer"
    [k0,k0Units] = check_units((1,NoUnits))

class aqueous(chemicalaffinity):
    """Aqueous food contact"""
    name = "aqueous contact"
    description = "minimize mass transfer"
    [k0,k0Units] = check_units((1000,NoUnits))

class intermediate(chemicalaffinity):
    """Intermediate chemical affinity"""
    name = "intermediate"
    description = "intermediate chemical affinity"
    [k0,k0Units] = check_units((10,NoUnits))

# Contact conditions

class chilled(realcontact):
    """real contact conditions"""
    description = "ambient storage conditions"
    name = "ambient"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((30,"days"))
    [contacttemperature,contacttemperatureUnits] = check_units((4,"degC"))

class frozen(realcontact):
    """real contact conditions"""
    description = "freezing storage conditions"
    name = "frrozen"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((6,"months"))
    [contacttemperature,contacttemperatureUnits] = check_units((-20,"degC"))

class ambient(realcontact):
    """real contact conditions"""
    description = "ambient storage conditions"
    name = "ambient"
    level = "contact"
    #[contacttime,contacttimeUnits] = check_units((200,"days"))
    #[contacttemperature,contacttemperatureUnits] = check_units((25,"degC"))

class hotfilled(realcontact):
    """real contact conditions"""
    description = "hot-filling conditions"
    name = "hotfilled"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((20,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((80,"degC"))

class boiling(realcontact):
    """real contact conditions"""
    description = "boiling conditions"
    name = "boiling"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((30,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((100,"degC"))

class pasteurization(realcontact):
    """real contact conditions"""
    description = "pasteurization conditions"
    name = "pasteurization"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((20,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((100,"degC"))

class sterilization(realcontact):
    """real contact conditions"""
    description = "sterilization conditions"
    name = "sterilization"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((20,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((121,"degC"))

class oven(realcontact):
    """real contact conditions"""
    description = "oven conditions"
    name = "oven"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((1,"hour"))
    [contacttemperature,contacttemperatureUnits] = check_units((180,"degC"))


# %% End-User classes
# -------------------------------------------------------------------
# THIRD LEVEL CLASSES
# Theses classes correspond to real cases and can be hybridized to
# derive new classes, for instance, for a specific brand of yoghurt.
# -------------------------------------------------------------------
class stacked(setoff):
    """stacked storage"""
    name = "stacked"
    description = "storage in stacks"
    level = "user"

class rolled(setoff):
    """rolled storage"""
    name = "rolled"
    description = "storage in rolls"
    level = "user"

class ethanol(simulant, perfectlymixed, fat):
    """Ethanol food simulant"""
    name = "ethanol"
    description = "ethanol = from pure ethanol down to ethanol 95%"
    level = "user"

class ethanol50(simulant, perfectlymixed, intermediate):
    """Ethanol 50 food simulant"""
    name = "ethanol 50"
    description = "ethanol 50, food simulant of dairy products"
    level = "user"

class water(simulant, perfectlymixed, aqueous):
    """Water food simulant"""
    name = "water"
    description = "water food simulant"
    level = "user"

class tenax(simulant, solid, fat):
    """Tenax(r) food simulant"""
    name = "Tenax"
    description = "simulant of dry food products"
    level = "user"

class yogurt(realfood, semisolid, ethanol50):
    """Yogurt as an example of real food"""
    description = "yogurt"
    level = "user"
    [k0,k0Units] = check_units((1,NoUnits))
    volume,volumeUnits = check_units((125,"mL"))

    # def __init__(self, name="no brand", volume=None, **kwargs):
    #     # Prepare a parameters dict: if a value is provided (e.g. volume), use it;
    #     # otherwise, the default (from class) is used.
    #     params = {}
    #     if volume is not None:
    #         params['volume'] = volume
    #     params['name'] = name
    #     params.update(kwargs)
    #     super().__init__(**params)

# -------------------------------------------------------------------
# Example usage (for debugging)
# -------------------------------------------------------------------
if __name__ == '__main__':
    F = foodlayer()
    E95 = ethanol()
    Y = yogurt()
    YF = yogurt(name="danone", volume=(150,"mL"))
    YF.description = "yogurt with fruits"  # You can still update the description on the instance if needed

    print("\n",repr(F),"\n"*2)
    print("\n",repr(E95),"\n"*2)
    print("\n",repr(Y),"\n"*2)
    print("\n",repr(YF),"\n"*2)

    # How to define a new food easily:
    class sandwich(realfood, solid, fat):
        name = "sandwich"
    S = sandwich()
    print("\n", repr(S))

    help_food()
