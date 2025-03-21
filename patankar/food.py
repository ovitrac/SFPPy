#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Food Layer
===============================================================================
Defines **food materials** for migration simulations. Models food as a **0D layer** with:
- **Mass transfer resistance (`h`)**
- **Partitioning (`k`)**
- **Contact time & temperature**

**Main Components:**
- **Base Class: `foodphysics`** (Stores all food-related parameters)
    - Defines mass transfer properties (`h`, `k`)
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

**Advanced Operators**
[substance] in [food] in [packaging] | [layer] >> [condition1] >> [condition2]


@version: 1.22
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2023-01-25
@rev: 2025-03-13

"""

# Dependencies
import sys
import inspect
import textwrap
import numpy as np
from copy import deepcopy as duplicate

from patankar.layer import check_units, NoUnits, layer # to convert units to SI
from patankar.loadpubchem import migrant

__all__ = ['acetonitrile', 'ambient', 'aqueous', 'boiling', 'check_units', 'chemicalaffinity', 'chilled', 'create_food_tree_widget', 'ethanol', 'ethanol50', 'ethanol95', 'fat', 'foodlayer', 'foodphysics', 'foodproperty', 'frozen', 'frying', 'get_defined_init_params', 'help_food', 'hotambient', 'hotfilled', 'hotoven', 'intermediate', 'is_valid_classname', 'isooctane', 'layer', 'liquid', 'list_food_classes', 'methanol', 'microwave', 'migrant', 'nofood', 'oil', 'oliveoil', 'oven', 'panfrying', 'pasteurization', 'perfectlymixed', 'realcontact', 'realfood', 'rolled', 'semisolid', 'setoff', 'simulant', 'solid', 'stacked', 'sterilization', 'tenax', 'testcontact', 'texture', 'transportation', 'update_class_list', 'water', 'water3aceticacid', 'wrap_text', 'yogurt']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.32"
#%% Private Properties and functions

# List of the default SI units used by physical quantity
parametersWithUnits = {"volume":"m**3",
                       "surfacearea":"m**2",
                       "density":"kg/m**3",
                       "contacttemperature":"degC",
                       "h":"m/s",
                       "k":NoUnits,  # user (preferred)
                       "k0":NoUnits, # alias (can be used after instantiation)
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

# class checked for the registry and other applications
def update_class_list(clslist, clsnew):
    """
    Updates clslist with clsnew if clsnew is a subclass of an existing class in clslist.
    If clsnew is a subclass of clslist[i], clslist[i] is replaced with clsnew.
    If no subclass relationship is found, clsnew is appended to clslist.

    Parameters:
        clslist (list): A list of classes.
        clsnew (type): A new class to check and possibly add.
    """
    for i, cls in enumerate(clslist):
        if issubclass(clsnew, cls):  # Check if clsnew is a subclass of clslist[i]
            clslist[i] = clsnew  # Replace the parent class with the new subclass
            return  # Exit immediately after updating

    clslist.append(clsnew)  # If no subclass relationship is found, append clsnew


# %% help food
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

# %% Widget
def create_food_tree_widget():
    """
    Creates a widget interface for food/contact conditions that uses a hierarchical
    tree to collect one class name per level. A custom class is then built via multiple
    inheritance from the classes selected at each level. This custom class is instantiated
    with overridden contact time and temperature values and stored in builtins.mycontacts.

    The tree is organized as follows (each branch has an associated 'class' and a 'children'
    dictionary):

    realfood
      ├─ liquid
      │    ├─ fat
      │    │     ├─ frozen
      │    │     ├─ chilled
      │    │     ├─ boiling
      │    │     ├─ pasteurization
      │    │     ├─ sterilization
      │    │     ├─ oven
      │    │     ├─ frying
      │    │     └─ hotoven
      │    ├─ aqueous
      │    │     ├─ frozen
      │    │     ├─ chilled
      │    │     ├─ boiling
      │    │     ├─ pasteurization
      │    │     └─ sterilization
      │    └─ intermediate
      │          ├─ frozen
      │          ├─ chilled
      │          ├─ boiling
      │          ├─ pasteurization
      │          ├─ sterilization
      │          ├─ oven
      │          ├─ frying
      │          └─ hotoven
      ├─ semisolid
      │      ... (similar structure)
      └─ solid
             ... (similar structure)

    simulant
      └─ (various branches; each branch has two leaves: realcontact and testcontact)

    setoff, stacked, rolled, nofood, yogurt, sandwich: with one or two levels

    Returns:
      An ipywidgets.VBox instance containing the full UI.
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError("ipywidgets and IPython are required for this interface.") from e

    import builtins
    # Ensure the global dictionary for storing instantiated food conditions exists.
    if not hasattr(builtins, "mycontacts"):
        builtins.mycontacts = {}
    global mycontacts
    mycontacts = builtins.mycontacts
    # flag for preheated gui interface (widgets should be initialized manually, instead of being empty)
    _preheatedGUI_ = hasattr(builtins, "_PREHEATED_") and getattr(builtins, "_PREHEATED_") is True

    # ---- Build the tree registry as a nested dict.
    # Each node is a dict with keys:
    #   "class": the class associated with that branch,
    #   "children": a dict of further branches (if any).
    food_tree = {
        "realfood": {
            "class": realfood,
            "children": {
                "liquid": {
                    "class": liquid,
                    "children": {
                        "fat": {
                            "class": fat,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        },
                        "aqueous": {
                            "class": aqueous,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization}
                            }
                        },
                        "intermediate": {
                            "class": intermediate,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        }
                    }
                },
                "semisolid": {
                    "class": semisolid,
                    "children": {
                        "fat": {
                            "class": fat,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        },
                        "aqueous": {
                            "class": aqueous,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven}
                            }
                        },
                        "intermediate": {
                            "class": intermediate,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        }
                    }
                },
                "solid": {
                    "class": solid,
                    "children": {
                        "fat": {
                            "class": fat,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        },
                        "aqueous": {
                            "class": aqueous,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        },
                        "intermediate": {
                            "class": intermediate,
                            "children": {
                                "ambient": {"class": ambient},
                                "frozen": {"class": frozen},
                                "chilled": {"class": chilled},
                                "boiling": {"class": boiling},
                                "pasteurization": {"class": pasteurization},
                                "sterilization": {"class": sterilization},
                                "oven": {"class": oven},
                                "frying": {"class": frying},
                                "hotoven": {"class": hotoven}
                            }
                        }
                    }
                }
            }
        },
        "simulant": {
            "class": simulant,
            "children": {
#                "fat": {"class": fat, "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
#                "aqueous": {"class": aqueous, },# "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
#                "intermediate": {"class": intermediate,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "oliveoil": {"class": oliveoil,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "ethanol95": {"class": ethanol95,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "ethanol50": {"class": ethanol50,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "acetonitrile": {"class": acetonitrile,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "methanol": {"class": methanol,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "water": {"class": water,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "water3aceticacid": {"class": water3aceticacid,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}},
                "isooctane": {"class": isooctane,}, # "children": {"realcontact": {"class": realcontact}, "testcontact": {"class": testcontact}}}
            }
        },
        "setoff": {
            "class": setoff,
            "children": {"ambient": {"class": ambient}, "transportation": {"class": transportation}, "hotambient": {"class": hotambient}}
        },
        "stacked": {
            "class": stacked,
            "children": {"ambient": {"class": ambient}, "transportation": {"class": transportation}, "hotambient": {"class": hotambient}}
        },
        "rolled": {
            "class": rolled,
            "children": {"ambient": {"class": ambient}, "transportation": {"class": transportation}, "hotambient": {"class": hotambient}}
        },
        "nofood": {
            "class": nofood,
            "children": {"ambient": {"class": ambient}, "transportation": {"class": transportation}, "hotambient": {"class": hotambient}}
        },
        "yogurt": {
            "class": yogurt,
            "children": {"chilled": {"class": chilled}, "ambient": {"class": ambient}}
        },
    }


    # ---- Create cascading dropdowns for up to 4 levels.
    dd1 = widgets.Dropdown(options=list(food_tree.keys()), description="Level 1:")
    dd2 = widgets.Dropdown(options=[], description="Level 2:")
    dd3 = widgets.Dropdown(options=[], description="Level 3:")
    dd4 = widgets.Dropdown(options=[], description="Level 4:")

    def update_dd2(*args):
        selected = dd1.value
        #print("dd2:",selected)
        opts = list(food_tree[selected]["children"].keys())
        dd2.options = opts
        dd2.value = opts[0] if opts else None
        update_dd3()
    def update_dd3(*args):
        if dd2.value is None:
            dd3.options = []
            dd3.value = None
        else:
            selected1 = dd1.value
            selected2 = dd2.value
            #print("dd3:",selected1,selected2) # debug
            subtree = food_tree[selected1]["children"][selected2]
            if "children" in subtree:
                opts = list(subtree["children"].keys())
                dd3.options = opts
                dd3.value = opts[0] if opts else None
            else:
                dd3.options = []
                dd3.value = None
        update_dd4()
    def update_dd4(*args):
        if dd3.value is None:
            dd4.options = []
            dd4.value = None
        else:
            selected1 = dd1.value
            selected2 = dd2.value
            selected3 = dd3.value
            #print("dd4:",selected1,selected2,selected3) # debug
            subtree = food_tree[selected1]["children"][selected2]["children"][selected3]
            if "children" in subtree:
                opts = list(subtree["children"].keys())
                dd4.options = opts
                dd4.value = opts[0] if opts else None
            else:
                dd4.options = []
                dd4.value = None
    dd1.observe(update_dd2, names='value')
    dd2.observe(update_dd3, names='value')
    dd3.observe(update_dd4, names='value')
    update_dd2()  # Initialize

    # ---- Additional inputs for overrides and instance name ----
    time_val = widgets.FloatText(value=60, description="Contact Time:")
    time_unit = widgets.Dropdown(options=["s", "min", "h", "days", "months", "years"],
                                 value="days", description="Time Unit:")
    temp_val = widgets.FloatText(value=25, description="Contact Temp:")
    temp_unit = widgets.Dropdown(options=["degC", "K"],
                                 value="degC", description="Temp Unit:")
    inst_name = widgets.Text(value="contact1", description="Instance Name:")

    # ---- Button to instantiate the custom food/contact condition ----
    btn = widgets.Button(description="Instantiate Custom Condition", button_style="success")
    out = widgets.Output()

    def instantiate_custom(b):
        with out:
            out.clear_output()
            selected_classes = []
            # Collect the class from Level 1
            try:
                cls1 = food_tree[dd1.value]["class"]
                selected_classes.append(cls1)
            except Exception as e:
                print("Error at Level 1:", e)
                return
            # Level 2 (if available)
            if dd2.value is not None:
                try:
                    cls2 = food_tree[dd1.value]["children"][dd2.value]["class"]
                    update_class_list(selected_classes,cls2) # avoid class conflicts (children replace parents)
                except Exception as e:
                    print("Error at Level 2:", e)
                    return
            # Level 3 (if available)
            if dd3.value is not None:
                try:
                    cls3 = food_tree[dd1.value]["children"][dd2.value]["children"][dd3.value]["class"]
                    update_class_list(selected_classes,cls3) # avoid class conflicts (children replace parents)
                except Exception as e:
                    print("Error at Level 3:", e)
                    return
            # Level 4 (if available)
            if dd4.options and dd4.value is not None:
                try:
                    cls4 = food_tree[dd1.value]["children"][dd2.value]["children"][dd3.value]["children"][dd4.value]["class"]
                    update_class_list(selected_classes,cls4) # avoid class conflicts (children replace parents)
                except Exception as e:
                    print("Error at Level 4:", e)
                    return
            # Dynamically create a custom class with multiple inheritance.
            try:
                CustomClass = type("CustomFood", tuple(selected_classes), {})
            except Exception as e:
                try:
                    # for debugging, we rely on level2 (situation fixed with update_class_list())
                    print("ERROR detected in foodtree:",e) # for debuging
                    [print(c.__name__) for c in selected_classes]
                    print("DEBUG:: keep only the last valid class:", selected_classes[-1])
                    CustomClass = type("CustomFood", tuple(selected_classes[-1:]), {})
                except Exception as e:
                    print("Error creating custom class:", e)
                    return
            # Instantiate the custom class with overrides.
            try:
                instance = CustomClass(contacttime=(time_val.value, time_unit.value),
                                         contacttemperature=(temp_val.value, temp_unit.value))
            except Exception as e:
                print("Error instantiating custom class:", e)
                return
            name_val = inst_name.value.strip()
            if not name_val:
                print("Please provide an instance name.")
                return
            instance.update(contacttime=(time_val.value,time_unit.value),
                            contactemperature=(temp_val.value,temp_unit.value))
            builtins.mycontacts[name_val] = instance
            print(f"Instantiated custom condition '{name_val}':")
            print(instance)
            print("\nCurrent conditions:", list(builtins.mycontacts.keys()))
    btn.on_click(instantiate_custom)

    if _preheatedGUI_:
        instantiate_custom(None) # we instantiate manually

    ui = widgets.VBox([
         widgets.HBox([dd1, dd2, dd3, dd4]),
         widgets.HBox([time_val, time_unit, temp_val, temp_unit]),
         inst_name,
         btn,
         out
    ])

    return ui



#%% Base physics class
# -------------------------------------------------------------------
# Base Class to convert class defaults to instance attributes
# -------------------------------------------------------------------
class foodphysics:
    """
    ===============================================================================
    SFPPy Module: Food Physics (Base Class)
    ===============================================================================
    `foodphysics` serves as the **base class** for all food-related objects in mass
    transfer simulations. It defines key parameters for food interaction with packaging
    materials and implements dynamic property propagation for simulation models.

    ------------------------------------------------------------------------------
    **Core Functionality**
    ------------------------------------------------------------------------------
    - Defines **mass transfer properties**:
      - `h`: Mass transfer coefficient (m/s)
      - `k`: Partition coefficient (dimensionless)
    - Implements **contact conditions**:
      - `contacttime`: Duration of food-packaging contact
      - `contacttemperature`: Temperature of the contact interface
    - Supports **inheritance and property propagation** to layers.
    - Provides **physical state representation** (`solid`, `liquid`, `gas`).
    - Allows **customization of mass transfer coefficients** via `kmodel`.

    ------------------------------------------------------------------------------
    **Key Properties**
    ------------------------------------------------------------------------------
    - `h`: Mass transfer coefficient (m/s) defining resistance at the interface.
    - `k`: Henry-like partition coefficient between the food and the material.
    - `contacttime`: Time duration of the packaging-food interaction.
    - `contacttemperature`: Temperature at the packaging interface (°C).
    - `surfacearea`: Contact surface area between packaging and food (m²).
    - `volume`: Volume of the food medium (m³).
    - `density`: Density of the food medium (kg/m³).
    - `substance`: The migrating substance (e.g., a chemical compound).
    - `medium`: The food medium in contact with packaging.
    - `kmodel`: Custom partitioning model (can be overridden by the user).

    ------------------------------------------------------------------------------
    **Methods**
    ------------------------------------------------------------------------------
    - `__rshift__(self, other)`: Propagates food properties to a layer (`food >> layer`).
    - `__matmul__(self, other)`: Equivalent to `>>`, enables `food @ layer`.
    - `migration(self, material, **kwargs)`: Simulates migration into a packaging layer.
    - `contact(self, material, **kwargs)`: Alias for `migration()`.
    - `update(self, **kwargs)`: Dynamically updates food properties.
    - `get_param(self, key, default=None, acceptNone=True)`: Retrieves a parameter safely.
    - `refresh(self)`: Ensures all properties are validated before simulation.
    - `acknowledge(self, what, category)`: Tracks inherited properties.
    - `copy(self, **kwargs)`: Creates a deep copy of the food object.

    ------------------------------------------------------------------------------
    **Integration with SFPPy Modules**
    ------------------------------------------------------------------------------
    - Works with `migration.py` to define the **left-side boundary condition**.
    - Interfaces with `layer.py` to apply contact temperature propagation.
    - Connects with `geometry.py` for food-contacting packaging surfaces.

    ------------------------------------------------------------------------------
    **Usage Example**
    ------------------------------------------------------------------------------
    ```python
    from patankar.food import foodphysics
    from patankar.layer import layer

    medium = foodphysics(contacttemperature=(40, "degC"), h=(1e-6, "m/s"))
    packaging_layer = layer(D=1e-14, l=50e-6)

    # Propagate food properties to the layer
    medium >> packaging_layer

    # Simulate migration
    from patankar.migration import senspatankar
    solution = senspatankar(packaging_layer, medium)
    solution.plotCF()
    ```

    ------------------------------------------------------------------------------
    **Notes**
    ------------------------------------------------------------------------------
    - The `foodphysics` class is the parent of `foodlayer`, `nofood`, `setoff`,
      `realcontact`, and `testcontact`.
    - The `PBC` property identifies periodic boundary conditions (used in `setoff`).
    - This class provides **dynamic inheritance** for mass transfer properties.

    """

    # General descriptors
    description = "Root physics class used to implement food and mass transfer physics"  # Remains as class attribute
    name = "food physics"
    level = "base"

    # Low-level prediction properties (F=contact medium, i=solute/migrant)
    # these @properties are defined by foodlayer, they should be duplicated
    _lowLevelPredictionPropertyList = [
        "chemicalsubstance","simulant","polarityindex","ispolymer","issolid", # F: common with patankar.layer
        "physicalstate","chemicalclass", # phase F properties
        "substance","migrant","solute", # i properties with synonyms substance=migrant=solute
        # users use "k", but internally we use k0, keep _kmodel in the instance
        "k0","k0unit","kmodel","_compute_kmodel", # Henry-like coefficients returned as properties with possible user override with medium.k0model=None or a function
        # SML properties
        "SML","SMLunit"
        ]

    # ------------------------------------------------------
    # Transfer rules for food1 >> food2 and food1 >> result
    # ------------------------------------------------------

    # Mapping of properties to their respective categories
    _list_categories = {
        "contacttemperature": "contact",
        "contacttime": "contact",
        "surfacearea": "geometry",
        "volume": "geometry",
        "substance": "substance",
        "medium": "medium"
    }

    # Rules for property transfer wtih >> or @ based on object type
    # ["property name"]["name of the destination class"][attr]
    #   - if onlyifinherited, only inherited values are transferred
    #   - if checkNmPy, the value will be transferred as a np.ndarray
    #   - name is the name of the property in the destination class (use "" to keep the same name)
    #   - prototype is the class itself (available only after instantiation, keep None here)
    _transferable_properties = {
        "contacttemperature": {
            "foodphysics": {
                "onlyifinherited": True,
                "checkNumPy": False,
                "as": "",
                "prototype": None,
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
                "prototype": None,
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
        },
        "substance": {
            "foodlayer": {
                "onlyifinherited": False,
                "checkNumPy": False,
                "as": "",
                "prototype": None,
            },
            "layer": {
                "onlyifinherited": False,
                "checkNumPy": False,
                "as": "",
                "prototype": None
            }
        },
        "medium": {
            "layer": {
                "onlyifinherited": False,
                "checkNumPy": False,
                "as": "",
                "prototype": None
            }
        },
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
                if key in ("description","level") or key in self._lowLevelPredictionPropertyList or key.startswith("__") or key.startswith("_") or callable(value):
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
        # we initialize _kmodel if _compute_kmodel exists
        if hasattr(self,"_compute_kmodel"):
            self._kmodel = "default" # do not initialize at self._compute_kmodel (default forces refresh)
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
            self.__class__._transferable_properties["substance"]["foodlayer"]["prototype"] = migrant
            self.__class__._transferable_properties["substance"]["layer"]["prototype"] = layer
            self.__class__._transferable_properties["medium"]["layer"]["prototype"] = layer
        # Subtance propagation (if needed)
        if self.hassubstance:
            if self.substance.hasSML:
                self.SML = self.substance.SML
                self.SMLunit = self.substance.SMLunit

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
            - k (float or tuple): Henry-like coefficient for the food (can be tuple like (1,"a.u.")).

        """
        if not kwargs:  # shortcut
            return self # for chaining
        def checkunits(value):
            """Helper function to convert physical quantities to SI."""
            if isinstance(value, tuple) and len(value) == 2:
                scale = check_units(value)[0]  # Convert to SI, drop unit
                if isinstance(scale,np.ndarray):
                    return scale
                else:
                    return np.array([scale], dtype=float)  # Ensure NumPy array
            elif isinstance(value, (int, float)):
                return np.array([value], dtype=float)  # Ensure NumPy array
            elif isinstance(value,np.ndarray):
                return value
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
        # Update medium, migrant (they accept aliases)
        lex = {
            "substance": ("substance", "migrant", "chemical", "solute"),
            "medium": ("medium", "simulant", "food", "contact"),
        }
        used_aliases = {}
        def get_value(canonical_key):
            """Find the correct alias in kwargs and return its value, or None if not found."""
            found_key = None
            for alias in lex.get(canonical_key, ()):  # Get aliases, default to empty tuple
                if alias in kwargs:
                    if alias in used_aliases:
                        raise ValueError(f"Alias '{alias}' is used for multiple canonical keys!")
                    found_key = alias
                    used_aliases[alias] = canonical_key
                    break  # Stop at the first match
            return kwargs.get(found_key, None)  # Return value if found, else None
        # Assign values only if found in kwargs
        new_substance = get_value("substance")
        new_medium = get_value("medium")
        if new_substance is not None: self.substance = new_substance
        if new_medium is not None:self.medium = new_medium
        # return
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
        """Formatted string representation of the FOODlayer object."""
        # Refresh all definitions
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
        # Collect defined properties and their formatted values
        properties = []
        excluded = ("k") if self.haskmodel else ("k0")
        for key, unit in parametersWithUnits.items():
            if hasattr(self, key) and key not in excluded:  # Include only defined parameters
                value = format_value(getattr(self, key))
                unit_str = self.get_param(f"{key}Units", unit)  # Retrieve unit safely
                if value is not None:
                    properties.append((key, f"{value:0.8g}", unit_str))

        # Sort properties alphabetically
        properties.sort(key=lambda x: x[0])

        # Determine max width for right-aligned names
        max_key_length = max(len(key) for key, _, _ in properties) if properties else 0
        # Construct formatted property list
        for key, value, unit_str in properties:
            repr_str += f"{key.rjust(max_key_length)}: {value} [{unit_str}]\n"
            if key == "k0":
                extra_info = f"{self._substance.k.__name__}(<{self.chemicalsubstance}>,{self._substance})"
                repr_str += f"{' ' * (max_key_length)}= {extra_info}\n"
        print(repr_str.strip())  # Print formatted output
        return str(self)  # Simplified representation for repr()



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
        Returns True if h is not defined or None
        This property is used to identified periodic boundary condition also called setoff mass transfer.

        """
        if not hasattr(self,"h"):
            return False # None
        htmp = getattr(self,"h")
        if isinstance(htmp,np.ndarray):
            htmp = htmp.item()
        return htmp is None

    @property
    def hassubstance(self):
        """Returns True if substance is defined (class migrant)"""
        if not hasattr(self, "_substance"):
            return False
        return isinstance(self._substance,migrant)



    # --------------------------------------------------------------------
    # For convenience, several operators have been overloaded
    #   medium >> packaging      # sets the volume and the surfacearea
    #   medium >> material       # propgates the contact temperature from the medium to the material
    #   medium in packaging    # simulate migration from the material to the medium
    #   medium << packaging
    # --------------------------------------------------------------------

    # method: medium._from(packaging) and its associatd operators in or <<
    def _from(self, other = None):
        """Inherit packaging properties"""
        from patankar.geometry import Packaging3D
        if not isinstance(other,Packaging3D):
            raise TypeError(f"other must be a Packaging3d or migrant not a {type(other).__name__}")
        return other._to(self) # we use the reverse definition

    def __lshift__(self,other):
        """
            overload <<
            food << packaging --> food
            full example: substance in food << packaging >> layer --> layer
        """
        return self._from(other)

    def __rmod__(self,other):
        """
            overload % as: migrant % food --> food
            full example: migrant in food in packaging | layer >> condition >> condition
        """
        from patankar.loadpubchem import migrant
        if not isinstance(other,migrant):
            raise TypeError(f"other must a migrant and not a {type(other).__name__}")
        self.update(substance=other)
        return self


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
        """
            Overloads >> to propagate to other.
            e.g. food>>layer --> layer
        """
        # inherit substance/migrant from other if self.migrant is None
        if isinstance(other,(layer,foodlayer)):
            if isinstance(self,foodlayer):
                if self.substance is None and other.substance is not None:
                    self.substance = other.substance
        return self._to(other) # propagates other

        # The operators @ and | are depreciated and should not be used anymore
    def __matmul__(self, other): # depreciated @
        """
            Overload @: equivalent to >> if other is a layer.
            e.g. food@layer --> layer (depreciated)
        """
        if not isinstance(other, layer):
            raise TypeError(f"Right operand must be a layer not a {type(other).__name__}")
        return self._to(other)  # propagates other

    def __or__(self, other): # depreciated |
        """
            Overload @: equivalent to >> if other is a layer.
            e.g. food|layer --> layer (depreciated)
        """
        if not isinstance(other, layer):
            raise TypeError(f"Right operand must be a layer not a {type(other).__name__}")
        return self._to(other)  # propagates other

    # migration method
    def migration(self,material,**kwargs):
        """interface to simulation engine: senspantankar"""
        from patankar.migration import senspatankar
        self._to(material) # propagate contact conditions first
        sim = senspatankar(material,self,**kwargs)
        self.lastsimulation = sim # store the last simulation result in medium
        self.lastinput = material # store the last input (material)
        sim.savestate(material,self) # store store the inputs in sim for chaining
        return sim

    def contact(self,material,**kwargs):
        """alias to migration method"""
        return self.migration(self,material,**kwargs)

    @property
    def haskmodel(self):
        """Returns True if a kmodel has been defined"""
        if hasattr(self, "_compute_kmodel"):
            if self._compute_kmodel() is not None:
                return True
            elif callable(self.kmodel):
                return self.kmodel() is not None
        return False

    @property
    def hasSML(self):
        """Returns True if SML has been defined"""
        return hasattr(self,"_SML") and self._SML is not None
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
    ===============================================================================
    SFPPy Module: Food Layer
    ===============================================================================
    `foodlayer` models food as a **0D layer** in mass transfer simulations, serving
    as the primary class for defining the medium in contact with a packaging material.

    ------------------------------------------------------------------------------
    **Core Functionality**
    ------------------------------------------------------------------------------
    - Models food as a **zero-dimensional (0D) medium** with:
      - A **mass transfer resistance (`h`)** at the interface.
      - A **partitioning behavior (`k`)** between food and packaging.
      - **Contact time (`contacttime`) and temperature (`contacttemperature`)**.
    - Defines **food geometry**:
      - `surfacearea`: Contact area with the material (m²).
      - `volume`: Total volume of the food medium (m³).
    - Supports **impervious (`nofood`) and periodic (`setoff`) conditions**.

    ------------------------------------------------------------------------------
    **Key Properties**
    ------------------------------------------------------------------------------
    - `h`: Mass transfer coefficient (m/s) defining resistance at the interface.
    - `k`: Partition coefficient describing substance solubility in food.
    - `contacttime`: Time duration of the packaging-food interaction.
    - `contacttemperature`: Temperature at the packaging interface (°C).
    - `surfacearea`: Contact surface area between packaging and food (m²).
    - `volume`: Volume of the food medium (m³).
    - `density`: Density of the food medium (kg/m³).
    - `substance`: Migrant (chemical) diffusing into food.
    - `medium`: Food medium in contact with packaging.
    - `impervious`: `True` if no mass transfer occurs (`nofood` class).
    - `PBC`: `True` if periodic boundary conditions apply (`setoff` class).

    ------------------------------------------------------------------------------
    **Methods**
    ------------------------------------------------------------------------------
    - `__rshift__(self, other)`: Propagates food properties to a packaging layer (`food >> layer`).
    - `__matmul__(self, other)`: Equivalent to `>>`, enables `food @ layer`.
    - `migration(self, material, **kwargs)`: Simulates migration into a packaging layer.
    - `contact(self, material, **kwargs)`: Alias for `migration()`.
    - `update(self, **kwargs)`: Dynamically updates food properties.
    - `get_param(self, key, default=None, acceptNone=True)`: Retrieves a parameter safely.
    - `refresh(self)`: Ensures all properties are validated before simulation.
    - `acknowledge(self, what, category)`: Tracks inherited properties.
    - `copy(self, **kwargs)`: Creates a deep copy of the food object.

    ------------------------------------------------------------------------------
    **Integration with SFPPy Modules**
    ------------------------------------------------------------------------------
    - Used as the **left-side boundary** in `migration.py` simulations.
    - Interacts with `layer.py` to propagate temperature and partitioning effects.
    - Interfaces with `geometry.py` for food-contacting packaging simulations.

    ------------------------------------------------------------------------------
    **Usage Example**
    ------------------------------------------------------------------------------
    ```python
    from patankar.food import foodlayer
    medium = foodlayer(name="ethanol", contacttemperature=(40, "degC"))

    from patankar.layer import LDPE
    packaging = LDPE(l=50e-6, D=1e-14)

    # Propagate food properties to the packaging
    medium >> packaging

    # Simulate migration
    from patankar.migration import senspatankar
    solution = senspatankar(packaging, medium)
    solution.plotCF()
    ```

    ------------------------------------------------------------------------------
    **Notes**
    ------------------------------------------------------------------------------
    - The `foodlayer` class extends `foodphysics` and provides a physical
      representation of food in contact with packaging.
    - Subclasses include:
      - `setoff`: Periodic boundary conditions (stacked packaging).
      - `nofood`: Impervious boundary (no mass transfer).
      - `realcontact`, `testcontact`: Standardized food contact conditions.
    - The `h` parameter determines if the medium is **well-mixed** or **diffusion-limited**.

    """
    level = "root"
    description = "root food class"  # Remains as class attribute
    name = "generic food layer"
    # -----------------------------------------------------------------------------
    # Class attributes that can be overidden in instances.
    # Their default values are set in classes and overriden with similar
    # instance properties with @property.setter.
    # These values cannot be set during construction, but only after instantiation.
    # A common scale for polarity index for solvents is from 0 to 10:
    #     - 0-3: Non-polar solvents (e.g., hexane)
    #     - 4-6: Moderately polar solvents (e.g., acetone)
    #     - 7-10: Polar solvents (e.g., water)
    # -----------------------------------------------------------------------------
    # These properties are essential for model predictions, they cannot be customized
    # beyond the rules accepted by the model predictors (they are not metadata)
    # note: similar attributes exist for patanaker.layer objects (similar possible values)
    _physicalstate = "liquid"   # solid, liquid (default), gas, porous
    _chemicalclass = "other"    # polymer, other (default)
    _chemicalsubstance = None   # None (default), monomer for polymers
    _polarityindex = 0.0        # polarity index (roughly: 0=hexane, 10=water)
    # -----------------------------------------------------------------------------
    # Class attributes duplicated as instance parameters
    # -----------------------------------------------------------------------------
    volume,volumeUnits = check_units((1,"dm**3"))
    surfacearea,surfaceareaUnits = check_units((6,"dm**2"))
    density,densityUnits = check_units((1000,"kg/m**3"))
    CF0,CF0units = check_units((0,NoUnits))  # initial concentration (arbitrary units)
    contacttime, contacttime_units = check_units((10,"days"))
    contactemperature,contactemperatureUnits = check_units((40,"degC"),ExpectedUnits="degC") # temperature in °C
    _substance = None # substance container / similar construction in pantankar.layer = migrant
    _k0model = None
    # -----------------------------------------------------------------------------
    # Getter methods for class/instance properties: same definitions as in patankar.layer (mandatory)
    # medium properties
    # -----------------------------------------------------------------------------
    # PHASE PROPERTIES  (attention chemicalsubstance=F substance, substance=i substance)
    @property
    def physicalstate(self): return self._physicalstate
    @property
    def chemicalclass(self): return self._chemicalclass
    @property
    def chemicalsubstance(self): return self._chemicalsubstance
    @property
    def simulant(self): return self._chemicalsubstance # simulant is an alias of chemicalsubstance
    @property
    def polarityindex(self): return self._polarityindex
    @property
    def ispolymer(self): return self.physicalstate == "polymer"
    @property
    def issolid(self): return self.solid == "solid"
    # SUBSTANCE/SOLUTE/MIGRANT properties  (attention chemicalsubstance=F substance, substance=i substance)
    @property
    def substance(self): return self._substance # substance can be ambiguous
    @property
    def migrant(self): return self.substance    # synonym
    @property
    def solute(self): return self.substance     # synonym
    # SML
    @property
    def SML(self): return self._SML if self.hasSML else None
    @property
    def SMLunit(self): return self._SMLunit if self.hasSML else None

    # -----------------------------------------------------------------------------
    # Setter methods for class/instance properties: same definitions as in patankar.layer (mandatory)
    # -----------------------------------------------------------------------------
    # PHASE PROPERTIES  (attention chemicalsubstance=F substance, substance=i substance)
    @physicalstate.setter
    def physicalstate(self,value):
        if value not in ("solid","liquid","gas","supercritical"):
            raise ValueError(f"physicalstate must be solid/liduid/gas/supercritical and not {value}")
        self._physicalstate = value
    @chemicalclass.setter
    def chemicalclass(self,value):
        if value not in ("polymer","other"):
            raise ValueError(f"chemicalclass must be polymer/oher and not {value}")
        self._chemicalclass= value
    @chemicalsubstance.setter
    def chemicalsubstance(self,value):
        if not isinstance(value,str):
            raise ValueError("chemicalsubtance must be str not a {type(value).__name__}")
        self._chemicalsubstance= value
    @simulant.setter
    def simulant(self,value):
        self.chemicalsubstance = value # simulant is an alias of chemicalcalsubstance
    @polarityindex.setter
    def polarityindex(self,value):
        if not isinstance(value,(float,int)):
            raise ValueError("polarity index must be float not a {type(value).__name__}")
        # rescaled to match predictions - standard scale [0,10.2] - predicted scale [0,7.12]
        return self._polarityindex * migrant("water").polarityindex/10.2
    # SUBSTANCE/SOLUTE/MIGRANT properties  (attention chemicalsubstance=F substance, substance=i substance)
    @substance.setter
    def substance(self,value):
        if isinstance(value,str):
            value = migrant(value)
        if not isinstance(value,migrant):
            raise TypeError(f"substance/migrant/solute must be a migrant not a {type(value).__name__}")
        self._substance = value
        if value.hasSML:
            self.SML = value.SML
            self.SMLunit = value.SMLunit
    @migrant.setter
    def migrant(self,value):
        self.substance = value
        if self.substance.hasSML:
            self.SML = self.substance.SML
            self.SMLunit = self.substance.SMLunit
    @solute.setter
    def solute(self,value):
        self.substance = value
    @SML.setter
    def SML(self,value):
        if isinstance(value,np.ndarray):
            value = value.item(0)
        if not isinstance(value,(float,int)):
            raise TypeError(f"SML value must be a float or int not a {type(value).__name__}")
        self._SML = float(value)
        self._SMLunit = "mg/kg"
    @SMLunit.setter
    def SMLunit(self,units):
        units = "mg/kg" if units is None else units
        if not isinstance(units,str):
            raise TypeError(f"SMLunit value must be a str not a {type(units).__name__}")
        self._SMLunit = units

    # -----------------------------------------------------------------------------
    # Henry-like coefficient k and its alias k0 (internal use)
    # -----------------------------------------------------------------------------
    #   - k is the name of the Henry-like property for food (as set and seen by the user)
    #   - k0 is the property operated by migration
    #   - k0 = k except if kmodel (lambda function) does not returns None
    #   - kmodel returns None if _substance is not set (proper migrant)
    #   - kmodel = None will override any existing kmodel
    #   - kmodel must be intialized to "default" to refresh its definition with self
    # note: The implementation is almost symmetric with kmodel in patankar.layer.
    # The main difference are:
    #    - food classes are instantiated by foodphysics
    #    - k is used to store the value of k0 (not _k or _k0)
    # -----------------------------------------------------------------------------
    # layer substance (of class migrant or None)
    # k0 and k0units (k and kunits are user inputs)
    @property
    def k0(self):
        ktmp = None
        if self.kmodel == "default": # default behavior
            ktmp = self._compute_kmodel()
        elif callable(self.kmodel): # user override (not the default function)
            ktmp = self.kmodel()
        if ktmp:
            return np.full_like(self.k, ktmp,dtype=np.float64)
        return self.k
    @k0.setter
    def k0(self,value):
        if not isinstance(value,(int,float,np.ndarray)):
            TypeError("k0 must be int, float or np.ndarray")
        if isinstance(self.k,int): self.k = float(self.k)
        self.k = np.full_like(self.k,value,dtype=np.float64)
    @property
    def kmodel(self):
        return self._kmodel
    @kmodel.setter
    def kmodel(self,value):
        if value is None or callable(value):
            self._kmodel = value
        else:
            raise ValueError("kmodel must be None or a callable function")
    @property
    def _compute_kmodel(self):
        """Return a callable function that evaluates k with updated parameters."""
        if not isinstance(self._substance,migrant) or self._substance.keval() is None or self.chemicalsubstance is None:
            return lambda **kwargs: None  # Return a function that always returns None
        template = self._substance.ktemplate.copy()
        # add solute (i) properties: Pi and Vi have been set by loadpubchem already
        template.update(ispolymer = False)
        def func(**kwargs):
            if self.chemicalsubstance:
                if "+" in self.chemicalsubstance: # mixture (e.g.: water + ethanol)
                    ks = [migrant(s) for s in self.chemicalsubstance.split("+")] # several k: ks
                    Pk = np.mean([k.polarityindex for k in ks]) # we average Pk assuming 50:50 mixure
                    Vk = np.mean([k.molarvolumeMiller for k in ks]) # we average Vk assuming 50:50 mixure
                else: # pure simulant
                    simulant = migrant(self.chemicalsubstance)
                    Pk = simulant.polarityindex
                    Vk = simulant.molarvolumeMiller
                template.update(Pk = Pk, Vk = Vk)
                k = self._substance.k.evaluate(**dict(template, **kwargs))
                return k
            else:
                self.k
        return func # we return a callable function not a value


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
    k = 1.0

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
    _physicalstate = "solid"    # it will be enforced if solid is defined first (see obj.mro())
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
    [k,kUnits] = check_units((1,NoUnits))

class aqueous(chemicalaffinity):
    """Aqueous food contact"""
    name = "aqueous contact"
    description = "minimize mass transfer"
    [k,kUnits] = check_units((1000,NoUnits))

class intermediate(chemicalaffinity):
    """Intermediate chemical affinity"""
    name = "intermediate"
    description = "intermediate chemical affinity"
    [k,kUnits] = check_units((10,NoUnits))

# Contact conditions

class frozen(realcontact):
    """real contact conditions"""
    description = "freezing storage conditions"
    name = "frrozen"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((6,"months"))
    [contacttemperature,contacttemperatureUnits] = check_units((-20,"degC"))

class chilled(realcontact):
    """real contact conditions"""
    description = "ambient storage conditions"
    name = "ambient"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((30,"days"))
    [contacttemperature,contacttemperatureUnits] = check_units((4,"degC"))

class ambient(realcontact):
    """real contact conditions"""
    description = "ambient storage conditions"
    name = "ambient"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((200,"days"))
    [contacttemperature,contacttemperatureUnits] = check_units((25,"degC"))

class transportation(realcontact):
    """hot transportation contact conditions"""
    description = "hot transportation storage conditions"
    name = "hot transportation"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((1,"month"))
    [contacttemperature,contacttemperatureUnits] = check_units((40,"degC"))

class hotambient(realcontact):
    """real contact conditions"""
    description = "hot ambient storage conditions"
    name = "hot ambient"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((2,"months"))
    [contacttemperature,contacttemperatureUnits] = check_units((40,"degC"))

class hotfilled(realcontact):
    """real contact conditions"""
    description = "hot-filling conditions"
    name = "hotfilled"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((20,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((80,"degC"))

class microwave(realcontact):
    """real contact conditions"""
    description = "microwave-oven conditions"
    name = "microwave"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((10,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((100,"degC"))

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

class panfrying(realcontact):
    """real contact conditions"""
    description = "panfrying conditions"
    name = "panfrying"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((20,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((120,"degC"))


class frying(realcontact):
    """real contact conditions"""
    description = "frying conditions"
    name = "frying"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((10,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((160,"degC"))

class oven(realcontact):
    """real contact conditions"""
    description = "oven conditions"
    name = "oven"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((1,"hour"))
    [contacttemperature,contacttemperatureUnits] = check_units((180,"degC"))

class hotoven(realcontact):
    """real contact conditions"""
    description = "hot oven conditions"
    name = "hot oven"
    level = "contact"
    [contacttime,contacttimeUnits] = check_units((30,"min"))
    [contacttemperature,contacttemperatureUnits] = check_units((230,"degC"))


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

class isooctane(simulant, perfectlymixed, fat):
    """Isoactane food simulant"""
    _chemicalsubstance = "isooctane"
    _polarityindex = 1.0 # Very non-polar hydrocarbon. Dielectric constant ~1.9.
    name = "isooctane"
    description = "isooctane food simulant"
    level = "user"

class oliveoil(simulant, perfectlymixed, fat):
    """Isoactane food simulant"""
    _chemicalsubstance = "methyl stearate"
    _polarityindex = 1.0 # Primarily triacylglycerides; still quite non-polar, though it contains some polar headgroups (the glycerol backbone).
    name = "olive oil"
    description = "olive oil food simulant"
    level = "user"
class oil(oliveoil): pass # synonym of oliveoil

class ethanol(simulant, perfectlymixed, fat):
    """Ethanol food simulant"""
    _chemicalsubstance = "ethanol"
    _polarityindex = 5.0 # Polar protic solvent; dielectric constant ~24.5. Lower polarity than methanol.
    name = "ethanol"
    description = "ethanol = from pure ethanol down to ethanol 95%"
    level = "user"
class ethanol95(ethanol): pass # synonym of ethanol

class ethanol50(simulant, perfectlymixed, intermediate):
    """Ethanol 50% food simulant"""
    _chemicalsubstance = "water+ethanol"
    _polarityindex = 7.0 # Intermediate polarity between ethanol and water.
    name = "ethanol 50"
    description = "ethanol 50, food simulant of dairy products"
    level = "user"

class acetonitrile(simulant, perfectlymixed, aqueous):
    """Acetonitrile food simulant"""
    _chemicalsubstance = "acetonitrile"
    _polarityindex = 6.8 # Polar aprotic solvent; dielectric constant ~36. Comparable to methanol in some polarity rankings.
    name = "acetonitrile"
    description = "acetonitrile"
    level = "user"

class methanol(simulant, perfectlymixed, aqueous):
    """Methanol food simulant"""
    _chemicalsubstance = "methanol"
    _polarityindex = 8.1 # Polar protic, dielectric constant ~33. Highly capable of hydrogen bonding, but still less so than water.
    name = "methanol"
    description = "methanol"
    level = "user"

class water(simulant, perfectlymixed, aqueous):
    """Water food simulant"""
    _chemicalsubstance = "water"
    _polarityindex = 10.2
    name = "water"
    description = "water food simulant"
    level = "user"

class water3aceticacid(simulant, perfectlymixed, aqueous):
    """Water food simulant"""
    _chemicalsubstance = "water"
    _polarityindex = 10.0 # Essentially still dominated by water’s polarity; 3% acetic acid does not drastically lower overall polarity.
    name = "water 3% acetic acid"
    description = "water 3% acetic acid - simulant for acidic aqueous foods"
    level = "user"

class tenax(simulant, solid, fat):
    """Tenax(r) food simulant"""
    _physicalstate = "porous"    # it will be enforced if tenax is defined first (see obj.mro())
    name = "Tenax"
    description = "simulant of dry food products"
    level = "user"

class yogurt(realfood, semisolid, ethanol50):
    """Yogurt as an example of real food"""
    description = "yogurt"
    level = "user"
    [k,kUnits] = check_units((1,NoUnits))
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
    food_tree_widget = create_food_tree_widget()
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
