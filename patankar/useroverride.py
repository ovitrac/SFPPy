
__all__ = ['UserOverride', 'useroverride']
"""
This module defines the UserOverride class which provides a dynamic mechanismto override SFPPy parameters.

The legacy system was using a global variable named 'SetSFPPy'. User scripts could override parameters by simply doing:

    SetSFPPy["param"] = value # this mechanism is depreciated and supported anymore by default (no more code injection)

The updated system (for SFPPy 1.40 and above) uses:
    from patankar.useroverride import useroverride
    useroverride.param = new value

The SFPPy code can then check for overrides via the provided methods.


Practical usage:
-----------------
    1) In one of your notebook cells, import the module and the global instance. This instance (named useroverride) is shared across all SFPPy modules. It can be injected into the global namespace as SetSFPPy.

    from patankar.useroverride import useroverride

⚠️❗ note: the instance useroverride is imported NOT the class UserOverride

useroverride <ENTER> will give you the current definitions
useroverride.param = value defines a new value

From the side of the program, useroverride.check("param",defaultvalue, constraints) is used to check whether the overriden value is acceptable or not. If not, the default value sets in the program will be used instead.

The check method can be omitted and seroverride("param",defaultvalue,...) works also (syntax used in SFPPy).

Legacy usage:
---------------
    2) Call the inject() method so that SetSFPPy is defined globally. You can also pass a dictionary with initial parameters if desired.

    # Inject without an update dictionary; SetSFPPy will be an empty dict initially.
    useroverride.inject()

    # Alternatively, provide an initial update:
    # useroverride.inject({"exp0": 2.718})

    3) Once injected, you can access and modify SetSFPPy from any cell. For instance, you can override parameters or check them:

    # Override parameters via item assignment
    SetSFPPy["param"] = value


Synopsis:
The intended design is that SFPPy's code accesses a single, predominantly read-only instance (via SetSFPPy) while user-side scripts rarely inject new values to update it. This ensures that all parts of SFPPy consistently refer to one global configuration, avoiding the potential pitfalls of having multiple conflicting instances. The mechanism is meant to decouple parameter management from function arguments, enabling a clear, centralized override process.


@version: 1.40
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-03-10
@rev: 2025-03-26
"""


import builtins
import numpy as np
from collections.abc import MutableMapping


__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.40"


if not hasattr(builtins, 'SetSFPPy'):
    # Fallback: create an empty configuration dictionary
    builtins.SetSFPPy = {}

class UserOverride(MutableMapping):
    r"""
    A dynamic override container for SFPPy parameters.

    This class implements a dictionary-like object that can be injected into
    the global (builtins) namespace as `SetSFPPy` so that all SFPPy modules
    share the same override settings. It supports:

      - Dynamic attribute and item access (e.g., `useroverride.param` or
        `useroverride["param"]`).
      - An `inject()` method to publish the current override dictionary globally.
      - A custom `__repr__` that prints a nicely tabulated list of parameters.
      - A `__str__` method returning a short summary.
      - A `check()` method to validate parameter values with type and range checks.
      - An `update()` method accepting multiple key/value pairs.
      - A messaging mechanism so that errors and warnings can be recorded.

    Usage example:

        from useroverride import useroverride
        # Optionally update parameters from a dictionary:
        useroverride.inject({"alpha": 3.14})
        # Dynamic access:
        useroverride.beta = 42
        print(useroverride["alpha"])
        # Check and validate a parameter:
        alpha_val = useroverride.check("alpha", 1.0, expected_type=float, valuemin=0.0, valuemax=10.0)
        # Update several parameters:
        useroverride.update(gamma=100, delta=[1, 2, 3])

    Parameters are stored internally and always available via the global
    variable `SetSFPPy` (after injection).

    Author: Olivier Vitrac
    Email: olivier.vitrac@gmail.com
    Project: SFPPy, Version: 1.40
    """

    def __init__(self):
        # Internal dictionary to hold override parameters.
        super().__setattr__("_data", {})
        # List to hold messages (warnings, errors, info)
        super().__setattr__("_messages", [])

    def inject(self, update_dict=None):
        r"""
        Injects the current override container into the global namespace as 'SetSFPPy'.

        Optionally, an update dictionary can be provided which will update the
        internal parameters.

        Parameters:
            update_dict (dict, optional): A dictionary with additional parameters
                to update the override container.

        Returns:
            bool: True if a parameter named "param" exists in the overrides,
                  False otherwise.

        Example:
            >>> useroverride.inject({"param": 10})
            True
        """
        if update_dict:
            if not isinstance(update_dict, dict):
                self.add_message("inject: update_dict must be a dict.", level="error")
            else:
                self._data.update(update_dict)
        # Inject the instance into the builtins so that SFPPy modules can access it
        builtins.SetSFPPy = self
        return "param" in self._data

    def add_message(self, message, level="info"):
        r"""
        Adds a message to the internal message list for errors, warnings, or info.

        Parameters:
            message (str): The message text.
            level (str): The level of the message. Can be 'info', 'warning', or 'error'.

        Example:
            >>> useroverride.add_message("Parameter X missing.", level="warning")
        """
        self._messages.append((level, message))

    def check(self, key, default, expected_type=None, nparray=True,
              valuemin=None, valuemax=None, valuelist=[], acceptNone=False):
        r"""
        Validates and returns the parameter value.

        If the parameter `key` is not in the overrides, returns the `default`
        value. If it is overridden, the value is checked for type conformity,
        optionally converted to a NumPy array (if it is a numeric list and `nparray`
        is True), and tested against range and allowed values constraints.

        Parameters:
            key (str): The parameter name.
            default: The default value to use if the parameter is not overridden.
            expected_type (type or tuple of types): Expected type(s) (default: float).
            nparray (bool): If True, converts a numeric list to a NumPy array of the
                specified type (default: True).
            valuemin: Minimum acceptable value (or None if not used).
            valuemax: Maximum acceptable value (or None if not used).
            valuelist (list): List of acceptable values. If non-empty, the value
                must be in this list.
            acceptNone (bool): If True, a value of None is accepted (default: False).

        Returns:
            The validated parameter value if the override exists and passes all checks,
            otherwise the provided default value.

        Example:
            >>> useroverride["alpha"] = [1.0, 2.0, 3.0]
            >>> alpha_val = useroverride.check("alpha", 0.0, expected_type=float, nparray=True,
            ...                                valuemin=0.0, valuemax=10.0)
            >>> print(alpha_val)  # Will print a NumPy array of floats.
        """
        if key not in self._data:
            return default

        value = self._data[key]

        # Check for None
        if value is None and not acceptNone:
            self.add_message(f"Parameter '{key}' is None but None is not accepted. Using default value.", level="warning")
            return default
        if expected_type is None:
            expected_type = type(default)


        # Convert numeric list to numpy array if requested
        if nparray and isinstance(value, list) and all(isinstance(x, (int, float)) for x in value):
            try:
                value = np.array(value, dtype=expected_type)
                self._data[key] = value  # update stored value with converted array
            except Exception as e:
                self.add_message(f"Conversion of parameter '{key}' to numpy array failed: {e}. Using default value.", level="error")
                return default

        # Check type
        if not isinstance(value, expected_type):
            self.add_message(f"Parameter '{key}' is not of expected type {expected_type}. Using default value.", level="warning")
            return default

        # Range check (works for scalars or array-like objects)
        try:
            if valuemin is not None:
                if hasattr(value, "min"):
                    if value.min() < valuemin:
                        self.add_message(f"Parameter '{key}' has minimum value {value.min()} lower than allowed {valuemin}. Using default.", level="warning")
                        return default
                else:
                    if value < valuemin:
                        self.add_message(f"Parameter '{key}' with value {value} is less than allowed minimum {valuemin}. Using default.", level="warning")
                        return default
            if valuemax is not None:
                if hasattr(value, "max"):
                    if value.max() > valuemax:
                        self.add_message(f"Parameter '{key}' has maximum value {value.max()} higher than allowed {valuemax}. Using default.", level="warning")
                        return default
                else:
                    if value > valuemax:
                        self.add_message(f"Parameter '{key}' with value {value} is greater than allowed maximum {valuemax}. Using default.", level="warning")
                        return default
        except Exception as e:
            self.add_message(f"Range check for parameter '{key}' failed: {e}. Using default.", level="error")
            return default

        # Allowed values check
        if valuelist:
            if isinstance(value, (list, np.ndarray)):
                # For list/array, check each element
                values = value if isinstance(value, list) else value.tolist()
                if not all(v in valuelist for v in values):
                    self.add_message(f"One or more elements of parameter '{key}' are not in the allowed list {valuelist}. Using default.", level="warning")
                    return default
            else:
                if value not in valuelist:
                    self.add_message(f"Parameter '{key}' with value {value} is not in the allowed list {valuelist}. Using default.", level="warning")
                    return default

        return value


    def update(self, **kwargs):
        r"""
        Update the override parameters using keyword arguments.

        Accepts keyword arguments in the form:

            update(alpha=3.14, beta=42, ...)

        The method updates the internal parameters accordingly.

        Example:
            >>> useroverride.update(alpha=3.14, beta=42)
        """
        self._data.update(kwargs)


    # --- Shortcut ---
    def __call__(self, key, default, expected_type=None, nparray=True,
                 valuemin=None, valuemax=None, valuelist=[], acceptNone=False):
        """
        Allows the UserOverride instance to be called directly as a shortcut to the check() method.

        If expected_type is not provided, it is inferred from the type of the default value.

        Example:
            useroverride("toto", 12)  # Will perform the same check as useroverride.check("toto", 12, expected_type=type(12))
        """
        if expected_type is None:
            expected_type = type(default)
        return self.check(key, default, expected_type=expected_type, nparray=nparray,
                          valuemin=valuemin, valuemax=valuemax, valuelist=valuelist, acceptNone=acceptNone)


    # --- Dictionary-like interface methods ---
    def __getitem__(self, key):
        return self._data.get(key, None)

    def __setitem__(self, key, value):
        self._data[key] = value

    def __delitem__(self, key):
        if key in self._data:
            del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def keys(self):
        r"""Return a set-like object providing a view on the override parameter keys."""
        return self._data.keys()

    def values(self):
        r"""Return an object providing a view on the override parameter values."""
        return self._data.values()

    def todict(self):
        r"""
        Returns a shallow copy of the override parameters as a plain dictionary.

        Returns:
            dict: A copy of the parameters.
        """
        return self._data.copy()

    # --- Dynamic attribute access ---
    def __getattr__(self, name):
        """
        Enables dynamic attribute access.

        If an attribute is not found in the instance's __dict__,
        this method checks the internal parameter dictionary.
        """
        if name in self._data:
            return self._data[name]
        raise AttributeError(f"'UserOverride' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        """
        Enables dynamic attribute setting.

        Attributes that are not internal (i.e. not '_data' or '_messages')
        are stored in the internal parameter dictionary.
        """
        if name in {"_data", "_messages"} or name.startswith("__"):
            super().__setattr__(name, value)
        else:
            self._data[name] = value

    # --- String representations ---
    def __repr__(self):
        r"""
        Returns a tabulated representation of the override parameters.

        The keys are right aligned and the values are left aligned.
        For example:

                  param1: value1
                  param2: value2
            anotherparam: anothervalue
        """
        if not self._data:
            return "<UserOverride: No parameters defined>"
        max_len = max(len(str(k)) for k in self._data)
        lines = []
        for k, v in self._data.items():
            lines.append(f"{str(k).rjust(max_len)}: {v}")
        print("\n".join(lines))
        return str(self)

    def __str__(self):
        r"""
        Returns a short summary of the override container.

        For example:
            <UserOverride: 3 parameters defined>
        """
        return f"<{type(self).__name__}: {len(self._data)} parameters defined>"



# Create a global instance to be used throughout SFPPy.
useroverride = UserOverride()


# Here a list of useful overrides for patankar.migration
useroverride.plotconfig = {
    "tscale": 24.0 * 3600, # days used as time scale
    "tunit": "days",       # units
    "lscale": 1e-6,        # length scale (from SI), here µm
    "lunit": "µm",         # units
    "Cscale": 1.0,         # concentration conversion
    "Cunit": "a.u."        # use "mg/kg" if you are sure of your units
    }
useroverride.ntimes = 1000      # number of stored simulation times (max=20000)
useroverride.timescale = "sqrt" # best for the first step
useroverride.RelTol=1e-6        # relative tolerance for integration of PDE in time
useroverride.AbsTol=1e-6        # absolute tolerance for integration of PDE in time
useroverride.deepcopy = None    # forcing False will have side effects (keep None to have overrides)
useroverride.nmax = 15          # number of concentration profiles per profile
useroverride.plotSML = None # keep it to None if not it will override all plotSML values in plotCF()

# Here a list of useful overrides for patankar.layer
useroverride.nmeshmin = 20 # number of minimal FV volumes per layer
useroverride.nmesh = 600   # total number of FV volumes in the assembly (the result will be ntimes x nmesh)


# Optionally (legacy), one could automatically inject the override container into builtins (set SFPPy):
#useroverride.inject()
