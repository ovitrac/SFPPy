
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

    2) plotconfig configuration with a flexible syntax
    Setting the units alone will change automatically the scale (for t and l, not for C0)

    useroverride.plotconfig(tunit="month",lunit="mm",Cunit="ppm"); # note ";" to suppress output
    useroverride.plotconfig.tunit = "week"
    useroverride.plotconfig.tunit = (1,"week")
    useroverride.plotconfig.tscale = lref**2/Dref # diffusion time based on lref and Dref
    useroverride.plotconfig.tunit = "-"
    useroverride.plotconfig["tunit"] = "-" # works also

Legacy usage:
---------------
    3) Call the inject() method so that SetSFPPy is defined globally. You can also pass a dictionary with initial parameters if desired.

    # Inject without an update dictionary; SetSFPPy will be an empty dict initially.
    useroverride.inject()

    # Alternatively, provide an initial update:
    # useroverride.inject({"exp0": 2.718})

    3) Once injected, you can access and modify SetSFPPy from any cell. For instance, you can override parameters or check them:

    # Override parameters via item assignment
    SetSFPPy["param"] = value


Synopsis:
The intended design is that SFPPy's code accesses a single, predominantly read-only instance (via SetSFPPy) while user-side scripts rarely inject new values to update it. This ensures that all parts of SFPPy consistently refer to one global configuration, avoiding the potential pitfalls of having multiple conflicting instances. The mechanism is meant to decouple parameter management from function arguments, enabling a clear, centralized override process.


@version: 1.41
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-03-10
@rev: 2025-03-30
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

# plotcondfig keys
plotconfig_keys = {"tscale", "tunit", "lscale", "lunit", "Cscale", "Cunit"}

# plotconfig Class container
class _PlotConfigDescriptor:
    """
    A callable property-like descriptor for plotconfig that supports dual usage:
        1. Getter: returns always-valid config
            cfg = obj.plotconfig
        2. Callable: updates and returns config
            cfg = obj.plotconfig(tscale=(1, "hours"), Cunit="mg/kg")
    """
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return self._Wrapper(instance)

    class _Wrapper:
        def __init__(self, instance):
            object.__setattr__(self, "instance", instance)


        def __call__(self, **kwargs):
            validated = self.instance.plotconfig_validator(self.instance._plotconfig, **kwargs)
            self.instance._plotconfig = validated
            return validated

        def __getattr__(self, name):
            if self.instance._plotconfig is None:
                self()  # triggers default init
            try:
                return self.instance._plotconfig[name]
            except KeyError:
                raise AttributeError(f"'plotconfig' has no key '{name}'")

        def __getitem__(self, key):
            if self.instance._plotconfig is None:
                self()
            return self.instance._plotconfig[key]

        def __setitem__(self, key, value):
            if self.instance._plotconfig is None:
                self()
            validated = self.instance.plotconfig_validator(self.instance._plotconfig, **{key: value})
            self.instance._plotconfig = validated

        def __setattr__(self, name, value):
            if name == "instance":
                object.__setattr__(self, name, value)
            else:
                self.__setitem__(name, value)

        def __repr__(self):
            return repr(self())

        def as_dict(self):
            return self()

        def __bool__(self):
            return True



# main class UserOverride (do not call it directly, import its instance useroverride)
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

    Special case for plotconfig (a flexible interface with validation has been implemented):
        useroverride = UserOverride()
        useroverride.plotconfig.tscale = 3600          # Attribute-style
        useroverride.plotconfig.tunit = "seconds"
        useroverride.plotconfig["Cunit"] = "mg/kg"     # Dict-style
        useroverride.plotconfig(tscale=(2, "days"))    # Callable-style
        useroverride.plotconfig(tunit="min",lunit="nm")# Callable style that will set tscale and lscale
        print(useroverride.plotconfig.tscale)          # Access as attribute
        print(useroverride.plotconfig["Cunit"])        # Access as key

    Parameters are stored internally and always available via the global
    variable `SetSFPPy` (after injection).

    Author: Olivier Vitrac
    Email: olivier.vitrac@gmail.com
    Project: SFPPy, Version: 1.40
    """

    plotconfig = _PlotConfigDescriptor()

    def __init__(self):
        # Internal dictionary to hold override parameters.
        super().__setattr__("_data", {})
        super().__setattr__("_plotconfig", None) # special container for plotconfig
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

        Example:
            >>> useroverride.inject({"param": 10})

        """
        if update_dict:
            if not isinstance(update_dict, dict):
                self.add_message("inject: update_dict must be a dict.", level="error")
            else:
                self._data.update(update_dict)
        # Inject the instance into the builtins so that SFPPy modules can access it
        builtins.SetSFPPy = self # we provide now the full instance, not only self._data

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
        if key == "plotconfig":      # shortcut / plotconfig is well protected with _PlotConfigDescriptor
            return self.plotconfig

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

        Special handling:
        -----------------
        - If 'plotconfig' is provided, it is validated and stored in _plotconfig.
        - If any of the individual plotconfig keys (tscale, tunit, lscale, lunit, Cscale, Cunit)
          are provided, they are passed to self.plotconfig(...) and removed from _data.

        All other keys are stored in the internal override dictionary (_data).

        Example:
            >>> useroverride.update(alpha=3.14, tscale=(1, "days"))
        """
        # 1. Handle full plotconfig
        if "plotconfig" in kwargs:
            cfg = kwargs.pop("plotconfig")
            self._plotconfig = self.plotconfig_validator(cfg)

        # 2. Handle individual plotconfig keys
        pc_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in plotconfig_keys}
        if pc_kwargs:
            self.plotconfig(**pc_kwargs)

        # 3. Update remaining keys into _data
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
        if key == "plotconfig":
            self._plotconfig = self.plotconfig_validator(value)
        elif key in plotconfig_keys:
            self.plotconfig(**{key: value})
        else:
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
        # 🛠 FIX: allow access to descriptors like plotconfig
        cls_attr = getattr(type(self), name, None)
        if hasattr(cls_attr, "__get__"):
            return cls_attr.__get__(self, type(self))
        raise AttributeError(f"'UserOverride' object has no attribute '{name}'")


    def __setattr__(self, name, value):
        """
        Enables dynamic attribute setting.

        Attributes that are not internal (i.e., not '_data', '_plotconfig', etc.)
        are stored in the internal parameter dictionary or delegated to plotconfig.
        """
        if name in {"_data", "_messages", "_plotconfig"} or name.startswith("__"):
            super().__setattr__(name, value)
        elif name == "plotconfig":
            self._plotconfig = self.plotconfig_validator(value)
        elif name in plotconfig_keys:
            self.plotconfig(**{name: value})
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
        lines.append(f"{'plotconfig'.rjust(max_len)}: {self.plotconfig}")
        print("\n".join(lines))
        return str(self)

    def __str__(self):
        r"""
        Returns a short summary of the override container.

        For example:
            <UserOverride: 3 parameters defined>
        """
        return f"<{type(self).__name__}: {len(self._data)} parameters defined>"

    # --- plotconfig validator ---
    @staticmethod
    def plotconfig_validator(plotconfig=None,
            tscale=None, tunit=None,
            lscale=None, lunit=None,
            Cscale=None, Cunit=None
        ):
        """
        Validate and complete a plot configuration dictionary for postprocessing.

        This static method takes an optional `plotconfig` dictionary and individual keyword overrides
        to build a fully defined configuration used to scale and label time, space, and concentration
        in plots.

        Key behavior:
        -------------
        - Only missing or `None` fields are filled using default values.
        - If `tscale` or `lscale` are provided as `(value, unit)` tuples, they override any value in `plotconfig`.
        - If only scalar values are available in `plotconfig` or as keyword arguments, corresponding unit keys
          (`tunit`, `lunit`) must also be present to enable conversion via `_toSI()`.
        - The function expects and returns SI-based scaling factors internally.

        Parameters:
        -----------
        plotconfig : dict or None
            An optional dictionary that may contain keys:
            "tscale", "tunit", "lscale", "lunit", "Cscale", "Cunit".

        tscale : float or (value, unit) tuple, optional
            Time scale. If given as a tuple, it takes precedence over `plotconfig`.

        tunit : str, optional
            Time unit. Used only if `tscale` is a scalar and needs conversion to SI.

        lscale : float or (value, unit) tuple, optional
            Length scale. If given as a tuple, it takes precedence over `plotconfig`.

        lunit : str, optional
            Length unit. Used only if `lscale` is a scalar and needs conversion to SI.

        Cscale : float, optional
            Concentration scaling factor. Defaults to 1.0 if not provided or found in `plotconfig`.

        Cunit : str, optional
            Concentration unit label. Defaults to "a.u.".

        Returns:
        --------
        dict
            A validated plot configuration dictionary with the following keys:
            - "tscale": float (SI)
            - "tunit" : str
            - "lscale": float (SI)
            - "lunit" : str
            - "Cscale": float
            - "Cunit" : str

        Note:
        -----
        (value, unit) tuples have precedence over existing or default values.
        Default values are only applied when fields are missing or None.
        """

        default = {
            "tscale": (1, "days"),
            "tunit": "days",
            "lscale": (1, "µm"),
            "lunit": "µm",
            "Cscale": 1.0,
            "Cunit": "a.u."
        }
        from patankar.layer import _toSI
        cfg = {} if plotconfig is None else plotconfig.copy()

        def _tofloat(x, keyname="value"):
            """
            Convert a value to float with the following rules:
            - int → float
            - float → float
            - ndarray → first element → float
            - otherwise: raise TypeError
            """
            import numpy as np
            if isinstance(x, (int, float)):
                return float(x)
            elif isinstance(x, np.ndarray):
                try:
                    return float(x.flat[0])  # works on any shape
                except Exception:
                    raise TypeError(f"{keyname} (ndarray) must contain at least one element.")
            else:
                raise TypeError(f"{keyname} must be a float, int, or ndarray, not {type(x).__name__}")

        # TIME SCALE
        if isinstance(tscale, tuple):
            cfg["tscale"] = _toSI(tscale).item()
            cfg["tunit"] = tscale[1]
        elif tscale is not None:
            cfg["tscale"] = tscale
            if tunit is not None:
                cfg["tunit"] = tunit
        elif tscale is None and tunit is not None:
            if not isinstance(tunit,str):
                TypeError(f"tunit must be str not a {type(tunit).__name__}")
            try:
                cfg["tscale"] = _toSI((1,tunit)).item()
            except AttributeError:
                print(f"Warning: the unit [ {tunit} ] is not convertible to SI.\nYou need to adjust also tscale.")
            cfg["tunit"] = tunit
        elif "tscale" not in cfg or cfg.get("tscale") is None:
            cfg["tscale"] = _toSI(default["tscale"]).item()
            cfg["tunit"] = default["tunit"]

        # LENGTH SCALE
        if isinstance(lscale, tuple):
            cfg["lscale"] = _toSI(lscale).item()
            cfg["lunit"] = lscale[1]
        elif lscale is not None:
            cfg["lscale"] = lscale
            if lunit is not None:
                cfg["lunit"] = lunit
        elif lscale is None and lunit is not None:
            if not isinstance(lunit,str):
                TypeError(f"tunit must be str not a {type(tunit).__name__}")
            try:
                cfg["lscale"] = _toSI((1,lunit)).item()
            except AttributeError:
                print(f"Warning: the unit [ {lunit} ] is not convertible to SI.\nYou need to adjust also lscale.")
            cfg["lunit"] = lunit
        elif "lscale" not in cfg or cfg.get("lscale") is None:
            cfg["lscale"] = _toSI(default["lscale"]).item()
            cfg["lunit"] = default["lunit"]

        # CONCENTRATION SCALE (no conversion, just scalars)
        cfg["Cscale"] = Cscale if Cscale is not None else cfg.get("Cscale", default["Cscale"])
        cfg["Cunit"]  = Cunit  if Cunit  is not None else cfg.get("Cunit",  default["Cunit"])

        # Final validation
        cfg["tscale"] = _tofloat(cfg["tscale"], "tscale")
        cfg["lscale"] = _tofloat(cfg["lscale"], "lscale")
        cfg["Cscale"] = _tofloat(cfg["Cscale"], "Cscale")

        # Return validated/updated plotconfig
        return cfg



# Create a global instance to be used throughout SFPPy.
useroverride = UserOverride()

# Here a list of useful overrides for patankar.migration
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

if __name__ == '__main__':
    useroverride.plotconfig(tscale=(1,'s'))
    print(useroverride.plotconfig)            # should print full validated dict
    print(useroverride._plotconfig)           # should match
    print(useroverride.plotconfig.tscale)     # should match _plotconfig["tscale"]
