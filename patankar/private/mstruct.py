#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Module: `struct.py`
Matlab-like structure class with extensions for parameter evaluation, file paths, and automatic management of dependencies in parameter definitions. This module provides the following key classes:

- **`struct`**: A flexible base class that mimics Matlab structures, offering dynamic field creation, indexing, concatenation, and field-level evaluation.
- **`param`**: Derived from `struct`, this class enables dynamic evaluation of fields based on interdependent definitions.
- **`paramauto`**: A further extension of `param` with automatic sorting and resolution of parameter dependencies during operations.
- **`pstr`**: A string subclass specialized for handling file paths and POSIX compatibility.

---

## Purpose
This module aims to streamline the creation and manipulation of structures for scientific computation, data management, and dynamic scripting, particularly in complex workflows.

---

## Key Features
- **Flexible Dynamic Structure**: Provides `struct` with field creation, deletion, and manipulation.
- **Parameter Evaluation**: Supports interdependent parameter evaluation with `param`.
- **Path and String Management**: Handles file paths and POSIX compliance with `pstr`.
- **Automatic Dependency Resolution**: Manages parameter dependencies automatically with `paramauto`.

---

## Evaluation Features (Updated for Pizza 1.0)
- **Dynamic Expressions**: Evaluate expressions within `${...}` placeholders or as standalone scalar expressions.
- **Matrix and Array Support**: Perform advanced operations such as matrix multiplication (`@`), transposition (`.T`), and slicing within `${...}`.
- **Safe Evaluation**: Eliminates the use of `eval`, using `safe_fstring()` and `SafeEvaluator` for secure computation.
- **Comprehensive Function Set**:
  - **Trigonometric Functions**: `sin`, `cos`, `tan`, etc.
  - **Exponential and Logarithmic**: `exp`, `log`, `sqrt`, etc.
  - **Random Functions**: `gauss`, `uniform`, `randint`, etc.
- **Error Handling**: Robust detection of undefined variables, invalid operations, and unsupported expressions.
- **Type Preservation**: Retains original data types (e.g., `float`, `numpy.ndarray`) for accuracy and further computation.
- **Custom Formatting**: Formats arrays and matrices for display with clear distinction between row/column vectors and higher-dimensional arrays.

The implementation in Pizza 1.0 ensures both flexibility and security, making it ideal for scenarios requiring dynamic parameter management and safe expression evaluation.

---

## Examples

### Basic Struct Usage
```python
from struct import struct

s = struct(a=1, b=2, c='${a} + ${b}')
s.a = 10
s["b"] = 5
delattr(s, "c")  # Delete a field
```

---

### Parameter Evaluation with `param`
```python
from struct import param

# Define parameters with dependencies
p = param(a=1, b='${a}*2', c='${b}+5')
evaluated = p.eval()  # Evaluate all fields dynamically
print(evaluated.c)  # Output: 7
```

---

### Path Management with `pstr`
```python
from struct import pstr

# Create and manipulate POSIX-compliant paths
path = pstr("/this/is/a/path/")
combined = path / "file.txt"
print(combined)  # Output: "/this/is/a/path/file.txt"
```

---

### Automatic Dependency Handling with `paramauto`
```python
from struct import paramauto

# Automatically resolve dependencies in parameters
pa = paramauto(a=1, b='${a}+1', c='${b}*2')
pa.disp()
# Output:
# -----------
#       a: 1
#       b: ${a}+1
#        = 2
#       c: ${b}*2
#        = 4
# -----------
```

---

### Evaluation Usage
```python
from pizza.private.mstruct import param
import numpy as np

p = param()
p.a = [1.0, 0.2, 0.03, 0.004]
p.b = np.array([p.a])
p.f = p.b.T @ p.b  # Matrix multiplication
p.g = "${a[1]}"    # Expression referencing `a`
p.h = "${b.T @ b}" # Matrix operation
print(p.eval())
```

---

Created on Sun Jan 23 14:19:03 2022
**Author**: Olivier Vitrac, AgroParisTech
"""

# revision history
# 2022-02-12 fix disp method for empty structures
# 2022-02-12 add type, format
# 2022-02-19 integration of the derived class param()
# 2022-02-20 code optimization, iterable class- major update
# 2022-02-26 clarify in the help the precedence s=s1+s2
# 2022-02-28 display nested structures
# 2022-03-01 implement value as list
# 2022-03-02 display correctly class names (not instances)
# 2022-03-04 add str()
# 2022-03-05 add __copy__ and __deepcopy__ methods
# 2022-03-05 AttributeError replaces KeyError in getattr() exceptions (required for for pdoc3)
# 2022-03-16 Prevent replacement/evaluation if the string is escaped \${parameter}
# 2022-03-19 add struct2dict(), dict2struct()
# 2022-03-20 add zip(), items()
# 2022-03-27 add __setitem__(), fromkeysvalues(), struct2param()
# 2022-03-28 add read() and write()
# 2022-03-29 fix protection for $var, $variable - add keysorted(), tostruct()
# 2022-03-30 specific display p"this/is/my/path" for pstr
# 2022-03-31 add dispmax()
# 2022-04-05 add check(), such that a.check(b) is similar to b+a
# 2022-04-09 manage None and [] values in check()
# 2022-05-14 s[:4], s[(3,5,2)] indexing a structure with a slice, list, tuple generates a substructure
# 2022-05-14 isempty (property) is TRUE for an empty structure
# 2022-05-15 __getitem__ and __set__item are now vectorized, add clear()
# 2022-05-16 add sortdefinitions(), isexpression, isdefined(), isstrdefined()
# 2022-05-16 __add__ and __iadd__ when called explicitly (not with + and +=) accept sordefinitions=True
# 2022-05-16 improved help, add tostatic() - v 0.45 (major version)
# 2022-05-17 new class paramauto() to simplify the management of multiple definitions
# 2022-05-17 catches most common errors in expressions and display explicit error messages - v 0.46
# 2023-01-18 add indexing as a dictionary s["a"] is the same as s.a - 0.461
# 2023-01-19 add % as comment instead of # to enable replacement
# 2023-01-27 param.eval() add % to freeze an interpretation (needed when a list is spanned as a string)
# 2023-01-27 struct.format() will replace {var} by ${var} if var is not defined
# 2023-08-11 display "" as <empty string> if evaluated
# 2024-09-06 add _returnerror as paramm class attribute (default=true) - dscript.lambdaScriptdata overrides it
# 2024-09-12 file management for all OS
# 2024-09-12 repr() improvements
# 2024-10-09 enable @property as attribute if _propertyasattribute is True
# 2024-10-11 add _callable__ and update()
# 2024-10-22 raises an error in escape() if s is not a string
# 2024-10-25 add dellatr()
# 2024-10-26 force silentmode to + and += operators
# 2024-12-08 fix help
# 2025-01-17 enable evaluation with ! and first recursion for lists (v1.002)
# 2025-01-18 fixes and explicit imports, better management of NumpPy arrays
# 2025-01-19 consolidation of slice handling, implicit evaluation and error handling (v1.003)
# 2025-01-20 full implementation and documentation of NumPy shorthands (v.1.004), use debug=True to show evaluation errors
# 2025-01-21 full implementation of 3D and 4D NumPy arrays (v.1.005): better parser (replace_matrix_shorthand) and display (format_array), see final example
# 2025-01-22 implement Matlab syntaxes for row and column vector, and matrices
# 2025-01-23 consolidation of nD matrices, expand ranges start:stop, start:step:stop à la Matlab
# 2025-01-30 "!" forces nemerical evaluation in lists '![1,2,"test","${a[1]}"]
# 2025-01-30 substitutions and calulations applied in lists (not expressions) [1, 2, 'test', '${a[1]}'] yields [1, 2, 'test', 1]
# 2025-01-31 better parsing rules for converting spaces into , (Matlab-like) (version 1.0051)
# 2025-01-31 p=param(a=..,b=..), p("a","b") returns the values of a and b as a struct, p.getval("a") returns the value of a
# 2025-02-01 better separation of param interpolation and global evaluation, consolidation of local evaluation: "the sum p is ${sum(p)}", [${n},${o}]" (version 1.0052)
# 2025-02-03 add _precision, fix generator for NumPy
# 2025-02-06 maintnance version (fixes and extensions) in agreement with the revised version of param_demo()
# 2025-02-07 add struct.numrepl(), major release (version 1.0054), improve evaluation error messages
# 2025-02-19 add importfrom(), add operator <<, add validkeys()
# 2025-02-19 fix strong slowdown when precedence rules are applied with +: defer sorting() for paramauto()

__project__ = "Pizza3"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "GPLv3"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.0062"


# %% Dependencies
# import types     # to check types (not required anymore since only builtin types are used)
import ast         # for safe evaluation (ast.literal_eval is used to evaluate strings starting with !)
import operator    # operators
import re          # regular expression
from pathlib import Path # for path managment (note that pstr uses its own logic)
from pathlib import PurePosixPath as PurePath
from copy import copy as duplicate # to duplicate objects
from copy import deepcopy as duplicatedeep # used by __deepcopy__()
# Import math functions
import builtins, math
import random
import numpy as np

__all__ = ['AttrErrorDict', 'SafeEvaluator', 'evaluate_with_placeholders', 'is_empty', 'is_literal_string', 'param', 'paramauto', 'pstr', 'struct']


# %% Private classes, functions, variables

_list_types = (list,tuple,np.ndarray) # list types recognized as such
_numeric_types = (int,float,str,list,tuple,np.ndarray, np.generic) # numeric types recognized as such


# List of functions recognized by SafeEvaluator()
# some functions might be unavailable with Python versions older than 3.10
_number_theoretic_funcs = [ # --- Number-Theoretic Functions
    "comb",      # Number of ways to choose k items from n items without repetition and without order
    "factorial", # n factorial
    "gcd",       # Greatest common divisor of the integer arguments
    "isqrt",     # Integer square root of a nonnegative integer n
    "lcm",       # Least common multiple of the integer arguments
    "perm"       # Number of ways to choose k items from n items without repetition and with order
]
_floating_arithmetic_funcs = [ # --- Floating Point Arithmetic
    "ceil",      # Ceiling of x, the smallest integer greater than or equal to x
    "fabs",      # Absolute value of x
    "floor",     # Floor of x, the largest integer less than or equal to x
#    "fma",       # Fused multiply-add operation: (x * y) + z (Python 3.13 required)
    "fmod",      # Remainder of division x / y
    "modf",      # Fractional and integer parts of x
    "remainder", # Remainder of x with respect to y
    "trunc"      # Integer part of x
]
_floating_manipulation_funcs = [ # --- Floating Point Manipulation Functions
    "copysign",  # Magnitude (absolute value) of x with the sign of y
    "frexp",     # Mantissa and exponent of x
    "isclose",   # Check if the values a and b are close to each other
    "isfinite",  # Check if x is neither an infinity nor a NaN
    "isinf",     # Check if x is a positive or negative infinity
    "isnan",     # Check if x is a NaN (not a number)
    "ldexp",     # x * (2**i), inverse of function frexp()
    "nextafter", # Floating-point value steps steps after x towards y
    "ulp"        # Value of the least significant bit of x
]
_power_exp_log_funcs = [ # --- Power, Exponential, and Logarithmic Functions
    "cbrt",      # Cube root of x
    "exp",       # e raised to the power x
    "exp2",      # 2 raised to the power x
    "expm1",     # e raised to the power x, minus 1
    "log",       # Logarithm of x to the given base (e by default)
    "log1p",     # Natural logarithm of 1+x (base e)
    "log2",      # Base-2 logarithm of x
    "log10",     # Base-10 logarithm of x
    "pow",       # x raised to the power y
    "sqrt"       # Square root of x
]
_summation_product_funcs = [ # --- Summation and Product Functions
    "dist",      # Euclidean distance between two points p and q given as an iterable of coordinates
    "fsum",      # Sum of values in the input iterable
    "hypot",     # Euclidean norm of an iterable of coordinates
    "prod",      # Product of elements in the input iterable with a start value
    "sumprod"    # Sum of products from two iterables p and q
]
_angular_conversion_funcs = [ # --- Angular Conversion Functions
    "degrees",   # Convert angle x from radians to degrees
    "radians"    # Convert angle x from degrees to radians
]
_trigonometric_funcs = [ # --- Trigonometric Functions
    "acos",      # Arc cosine of x
    "asin",      # Arc sine of x
    "atan",      # Arc tangent of x
    "atan2",     # Arc tangent of y/x
    "cos",       # Cosine of x
    "sin",       # Sine of x
    "tan"        # Tangent of x
]
_hyperbolic_funcs = [ # --- Hyperbolic Functions
    "acosh",     # Inverse hyperbolic cosine of x
    "asinh",     # Inverse hyperbolic sine of x
    "atanh",     # Inverse hyperbolic tangent of x
    "cosh",      # Hyperbolic cosine of x
    "sinh",      # Hyperbolic sine of x
    "tanh"       # Hyperbolic tangent of x
]
_special_funcs = [ # --- Special Functions
    "erf",       # Error function at x
    "erfc",      # Complementary error function at x
    "gamma",     # Gamma function at x
    "lgamma"     # Natural logarithm of the absolute value of the Gamma function at x
]
_constants = [ # --- Constants
    "pi",        # π = 3.141592…
    "e",         # e = 2.718281…
    "tau",       # τ = 2π = 6.283185…
    "inf",       # Positive infinity
    "nan"        # “Not a number” (NaN)
]
# Combine all functions and constants into one list:
_all_math_names = (
    _number_theoretic_funcs +
    _floating_arithmetic_funcs +
    _floating_manipulation_funcs +
    _power_exp_log_funcs +
    _summation_product_funcs +
    _angular_conversion_funcs +
    _trigonometric_funcs +
    _hyperbolic_funcs +
    _special_funcs +
    _constants
)


# Safe f"" to evaluate ${var}, ${expression} and some expressions ${v1}+${v2}
class SafeEvaluator(ast.NodeVisitor):
    """A safe evaluator class for expressions involving math, NumPy, random, and basic operators."""

    def __init__(self, context={}):
        self.context = {**context}
        # Update context with math functions/constants from _all_math_names that exist in math module.
        self.context.update({
            name: getattr(math, name)
            for name in _all_math_names if hasattr(math, name)
        })
        # Add built-in functions relevant for math that are not part of the math module.
        for name in ("abs", "round", "min", "max", "sum", "divmod"):
            self.context[name] = getattr(builtins, name)

        self.context.update({
            "gauss": random.gauss,
            "uniform": random.uniform,
            "randint": random.randint,
            "choice": random.choice
        })
        # Add NumPy as np
        self.context["np"] = np  # Allow 'np.sin', 'np.cos', etc.
        # Define allowed operators
        self.operators = {
            ast.Add: operator.add,
            ast.Sub: operator.sub,
            ast.Mult: operator.mul,
            ast.Div: operator.truediv,
            ast.FloorDiv: operator.floordiv,
            ast.Mod: operator.mod,
            ast.Pow: operator.pow,
            ast.USub: operator.neg,  # Unary subtraction
        }

    def visit_Name(self, node):
        if node.id in self.context:
            return self.context[node.id]
        raise ValueError(f"Variable or function '{node.id}' is not defined")

    def visit_Constant(self, node):
        return node.value

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        op_type = type(node.op)
        if isinstance(left, np.ndarray) and isinstance(right, np.ndarray) and isinstance(node.op, ast.MatMult):
            return np.matmul(left, right)
        if op_type in self.operators:
            return self.operators[op_type](left, right)
        raise ValueError(f"Unsupported operator: {op_type}")

    def visit_UnaryOp(self, node):
        operand = self.visit(node.operand)
        op_type = type(node.op)
        if op_type in self.operators:
            return self.operators[op_type](operand)
        raise ValueError(f"Unsupported unary operator: {op_type}")

    def visit_Call(self, node):
        func = self.visit(node.func)
        if callable(func):
            args = [self.visit(arg) for arg in node.args]
            kwargs = {kw.arg: self.visit(kw.value) for kw in node.keywords}
            return func(*args, **kwargs)
        raise ValueError(f"Function '{ast.dump(node.func)}' is not callable")

    def visit_Attribute(self, node):
        value = self.visit(node.value)
        attr = node.attr
        if hasattr(value, attr):
            # If the attribute is "T", return the transpose of the array
            if attr == "T" and isinstance(value, np.ndarray):
                return value.T
            # Check if the attribute is the '@' matrix multiplication operator
            if attr == "@" and isinstance(value, np.ndarray):
                return value @ value  # or handle accordingly with another operand
            return getattr(value, attr)
        raise ValueError(f"Object '{value}' has no attribute '{attr}'")

    def visit_Subscript(self, node):
        value = self.visit(node.value)
        slice_obj = self.visit(node.slice)
        try:
            return value[slice_obj]
        except Exception as e:
            raise ValueError(f"Invalid index {slice_obj} for object of type {type(value).__name__}: {e}")

    def visit_Index(self, node):
        return self.visit(node.value)

    def visit_Slice(self, node):
        lower = self.visit(node.lower) if node.lower else None
        upper = self.visit(node.upper) if node.upper else None
        step = self.visit(node.step) if node.step else None
        return slice(lower, upper, step)

    def visit_ExtSlice(self, node):
        dims = tuple(self.visit(dim) for dim in node.dims)
        return dims

    def visit_Tuple(self, node):
        return tuple(self.visit(elt) for elt in node.elts)

    def visit_List(self, node):
        return [self.visit(elt) for elt in node.elts]

    def visit_Dict(self, node):
        """
        Evaluate a dictionary expression by safely evaluating each key and value.
        This allows expressions like: {"a": ${v1}+${v2}, "b": ${var}}.
        """
        return {self.visit(key): self.visit(value) for key, value in zip(node.keys, node.values)}

    def generic_visit(self, node):
        raise ValueError(f"Unsupported expression: {ast.dump(node)}")

    def evaluate(self, expression):
        tree = ast.parse(expression, mode='eval')
        return self.visit(tree.body)

# Use SafeEvaluator only to unescaped variables ${var} and expressions ${var1+var2}
# This methodology separates string interpolation from expression evaluation
def evaluate_with_placeholders(text, evaluator, evaluator_nocontext=SafeEvaluator(),raiseerror=False):
    """
    Evaluates only unescaped placeholders of the form ${...} in the input text.
    Escaped placeholders (\\${...}) are left as literal text (after removing the escape).

    Note1 : ${ ... } can be ${var} or an expressions such as ${var1+var2}
    Note2 : a full evaluation is attempted only after the full evaluation using the same
    evaluator without context

    Example:

        context = {'a': 10, 'b': 5}
        evaluator = SafeEvaluator(context)
        evaluator_nocontext = SafeEvaluator()

        text = "Evaluated variable: ${a} and literal: \\${a} and sum: ${a + b}, leave intact a+b, ${a}+${b}"
        processed_text = evaluate_with_placeholders(text, evaluator)
        print(processed_text)

    """
    # Pattern explanation:
    #   (?<!\\)   : Negative lookbehind to ensure the '${' is not preceded by a backslash.
    #   (\$\{([^}]+)\}) : Captures the full placeholder in group 1 and the inner expression in group 2.
    pattern = r'(?<!\\)(\$\{([^}]+)\})'

    def replace_placeholder(match):
        # Extract the inner expression from the placeholder.
        expr = match.group(2)
        # Evaluate the expression using your SafeEvaluator instance.
        try:
            value = evaluator.evaluate(expr)
        except Exception as e:
            raise ValueError(f"Error evaluating expression '{expr}': {e}")
        return str(value)

    # String Interpolation: Replace unescaped placeholders with their evaluated results.
    result = re.sub(pattern, replace_placeholder, text)
    # Finally, unescape the escaped placeholders: replace "\${" with "${".
    result = result.replace(r'\${', '${')
    # Full evaluation if possible
    if raiseerror:
        return evaluator_nocontext.evaluate(result)
    else:
        try:
            return evaluator_nocontext.evaluate(result)
        except Exception:
            return result


# returns True for literal string starting with "$"
def is_literal_string(s):
    """
    Returns True if the first non-blank character in the string is '$'
    and it is not immediately followed by '{' or '['.

    Parameters:
        s (str): The string to check.

    Returns:
        bool: True if the condition is met, otherwise False.
    """
    stripped = s.lstrip()  # Remove leading whitespace
    if not stripped:
        return False
    if stripped[0] != '$':
        return False
    # If there is a character following '$', ensure it's not '{' or '['.
    if len(stripped) > 1 and stripped[1] in ('{', '['):
        return False
    return True

# Class to handle expressions containing operators correctly without being misinterpreted as attribute accesses.
class AttrErrorDict(dict):
    """Custom dictionary that raises AttributeError (as required for the logic of struct)
       instead of KeyError for missing keys and strips quotes from strings."""
    def __getitem__(self, key):
        try:
            value = super().__getitem__(key)
            if isinstance(value, str):
                # Strip surrounding single or double quotes if present
                if (value.startswith("'") and value.endswith("'")) or (value.startswith('"') and value.endswith('"')):
                    return value[1:-1]
                return value
            return value
        except KeyError:
            raise AttributeError(f"Attribute '{key}' not found")

# A generalized isempty defintion
def is_empty(value):
    """Return True if value is considered empty (None, "", [] or ())."""
    return value is None or value == "" or value == [] or value == ()

# %% core struct class
class struct():
    """
    Class: `struct`
    ================

    A lightweight class that mimics Matlab-like structures, with additional features
    such as dynamic field creation, indexing, concatenation, and compatibility with
    evaluated parameters (`param`).

    ---

    ### Features
    - Dynamic creation of fields.
    - Indexing and iteration support for fields.
    - Concatenation and subtraction of structures.
    - Conversion to and from dictionaries.
    - Compatible with `param` and `paramauto` for evaluation and dependency handling.

    ---

    ### Examples

    #### Basic Usage
    ```python
    s = struct(a=1, b=2, c='${a} + ${b} # evaluate me if you can')
    print(s.a)  # 1
    s.d = 11    # Append a new field
    delattr(s, 'd')  # Delete the field
    ```

    #### Using `param` for Evaluation
    ```python
    p = param(a=1, b=2, c='${a} + ${b} # evaluate me if you can')
    p.eval()
    # Output:
    # --------
    #      a: 1
    #      b: 2
    #      c: ${a} + ${b} # evaluate me if you can (= 3)
    # --------
    ```

    ---

    ### Concatenation and Subtraction
    Fields from the right-most structure overwrite existing values.
    ```python
    a = struct(a=1, b=2)
    b = struct(c=3, d="d", e="e")
    c = a + b
    e = c - a
    ```

    ---

    ### Practical Shorthands

    #### Constructing a Structure from Keys
    ```python
    s = struct.fromkeys(["a", "b", "c", "d"])
    # Output:
    # --------
    #      a: None
    #      b: None
    #      c: None
    #      d: None
    # --------
    ```

    #### Building a Structure from Variables in a String
    ```python
    s = struct.scan("${a} + ${b} * ${c} / ${d} --- ${ee}")
    s.a = 1
    s.b = "test"
    s.c = [1, "a", 2]
    s.generator()
    # Output:
    # X = struct(
    #      a=1,
    #      b="test",
    #      c=[1, 'a', 2],
    #      d=None,
    #      ee=None
    # )
    ```

    #### Indexing and Iteration
    Structures can be indexed or sliced like lists.
    ```python
    c = a + b
    c[0]      # Access the first field
    c[-1]     # Access the last field
    c[:2]     # Slice the structure
    for field in c:
        print(field)
    ```

    ---

    ### Dynamic Dependency Management
    `struct` provides control over dependencies, sorting, and evaluation.

    ```python
    s = struct(d=3, e="${c} + {d}", c='${a} + ${b}', a=1, b=2)
    s.sortdefinitions()
    # Output:
    # --------
    #      d: 3
    #      a: 1
    #      b: 2
    #      c: ${a} + ${b}
    #      e: ${c} + ${d}
    # --------
    ```

    For dynamic evaluation, use `param`:
    ```python
    p = param(sortdefinitions=True, d=3, e="${c} + ${d}", c='${a} + ${b}', a=1, b=2)
    # Output:
    # --------
    #      d: 3
    #      a: 1
    #      b: 2
    #      c: ${a} + ${b}  (= 3)
    #      e: ${c} + ${d}  (= 6)
    # --------
    ```

    ---

    ### Overloaded Methods and Operators
    #### Supported Operators
    - `+`: Concatenation of two structures (`__add__`).
    - `-`: Subtraction of fields (`__sub__`).
    - `<<`: Import values from another structure (`__lshift__`)
    - `len()`: Number of fields (`__len__`).
    - `in`: Check for field existence (`__contains__`).

    #### Method Overview
    | Method                | Description                                             |
    |-----------------------|---------------------------------------------------------|
    | `check(default)`      | Populate fields with defaults if missing.               |
    | `clear()`             | Remove all fields.                                      |
    | `dict2struct(dico)`   | Create a structure from a dictionary.                   |
    | `disp()`              | Display the structure.                                  |
    | `eval()`              | Evaluate expressions within fields.                     |
    | `fromkeys(keys)`      | Create a structure from a list of keys.                 |
    | `generator()`         | Generate Python code representing the structure.        |
    | `importfrom()`        | Import undefined values from another struct or dict.    |
    | `items()`             | Return key-value pairs.                                 |
    | `keys()`              | Return all keys in the structure.                       |
    | `read(file)`          | Load structure fields from a file.                      |
    | `scan(string)`        | Extract variables from a string and populate fields.    |
    | `sortdefinitions()`   | Sort fields to resolve dependencies.                    |
    | `struct2dict()`       | Convert the structure to a dictionary.                  |
    | `validkeys()`         | Return valid keys                                       |
    | `values()`            | Return all field values.                                |
    | `write(file)`         | Save the structure to a file.                           |

    ---

    ### Dynamic Properties
    | Property    | Description                            |
    |-------------|----------------------------------------|
    | `isempty`   | `True` if the structure is empty.      |
    | `isdefined` | `True` if all fields are defined.      |

    ---
    """

    # attributes to be overdefined
    _type = "struct"        # object type
    _fulltype = "structure" # full name
    _ftype = "field"        # field name
    _evalfeature = False    # true if eval() is available
    _maxdisplay = 40        # maximum number of characters to display (should be even)
    _propertyasattribute = False
    _precision = 4
    _needs_sorting = False

    # attributes for the iterator method
    # Please keep it static, duplicate the object before changing _iter_
    _iter_ = 0

    # excluded attributes (keep the , in the Tupple if it is singleton)
    _excludedattr = {'_iter_','__class__','_protection','_evaluation','_returnerror','_debug','_precision','_needs_sorting'} # used by keys() and len()


    # Methods
    def __init__(self,debug=False,**kwargs):
        """ constructor, use debug=True to report eval errors"""
        # Optionally extend _excludedattr here
        self._excludedattr = self._excludedattr | {'_excludedattr', '_type', '_fulltype','_ftype'} # addition 2024-10-11
        self._debug = debug
        self.set(**kwargs)

    def zip(self):
        """ zip keys and values """
        return zip(self.keys(),self.values())

    @staticmethod
    def dict2struct(dico,makeparam=False):
        """ create a structure from a dictionary """
        if isinstance(dico,dict):
            s = param() if makeparam else struct()
            s.set(**dico)
            return s
        raise TypeError("the argument must be a dictionary")

    def struct2dict(self):
        """ create a dictionary from the current structure """
        return dict(self.zip())

    def struct2param(self,protection=False,evaluation=True):
        """ convert an object struct() to param() """
        p = param(**self.struct2dict())
        for i in range(len(self)):
            if isinstance(self[i],pstr): p[i] = pstr(p[i])
        p._protection = protection
        p._evaluation = evaluation
        return p

    def set(self,**kwargs):
        """ initialization """
        self.__dict__.update(kwargs)

    def setattr(self,key,value):
        """ set field and value """
        if isinstance(value,list) and len(value)==0 and key in self:
            delattr(self, key)
        else:
            self.__dict__[key] = value

    def getattr(self,key):
        """Get attribute override to access both instance attributes and properties if allowed."""
        if key in self.__dict__:
            return self.__dict__[key]
        elif getattr(self, '_propertyasattribute', False) and \
             key not in self._excludedattr and \
             key in self.__class__.__dict__ and isinstance(self.__class__.__dict__[key], property):
            # If _propertyasattribute is True and it's a property, get its value
            return self.__class__.__dict__[key].fget(self)
        else:
            raise AttributeError(f'the {self._ftype} "{key}" does not exist')

    def hasattr(self, key):
        """Return true if the field exists, considering properties as regular attributes if allowed."""
        return key in self.__dict__ or (
            getattr(self, '_propertyasattribute', False) and
            key not in self._excludedattr and
            key in self.__class__.__dict__ and isinstance(self.__class__.__dict__[key], property)
        )

    def __getstate__(self):
        """ getstate for cooperative inheritance / duplication """
        return self.__dict__.copy()

    def __setstate__(self,state):
        """ setstate for cooperative inheritance / duplication """
        self.__dict__.update(state)

    def __getattr__(self,key):
        """ get attribute override """
        return pstr.eval(self.getattr(key))

    def __setattr__(self,key,value):
        """ set attribute override """
        self.setattr(key,value)

    def __contains__(self,item):
        """ in override """
        return self.hasattr(item)

    def keys(self):
        """ return the fields """
        # keys() is used by struct() and its iterator
        return [key for key in self.__dict__.keys() if key not in self._excludedattr]

    def keyssorted(self,reverse=True):
        """ sort keys by length() """
        klist = self.keys()
        l = [len(k) for k in klist]
        return [k for _,k in sorted(zip(l,klist),reverse=reverse)]

    def values(self):
        """ return the values """
        # values() is used by struct() and its iterator
        return [pstr.eval(value) for key,value in self.__dict__.items() if key not in self._excludedattr]

    @staticmethod
    def fromkeysvalues(keys,values,makeparam=False):
        """ struct.keysvalues(keys,values) creates a structure from keys and values
            use makeparam = True to create a param instead of struct
        """
        if keys is None: raise AttributeError("the keys must not empty")
        if not isinstance(keys,_list_types): keys = [keys]
        if not isinstance(values,_list_types): values = [values]
        nk,nv = len(keys), len(values)
        s = param() if makeparam else struct()
        if nk>0 and nv>0:
            iv = 0
            for ik in range(nk):
                s.setattr(keys[ik], values[iv])
                iv = min(nv-1,iv+1)
            for ik in range(nk,nv):
                s.setattr(f"key{ik}", values[ik])
        return s

    def items(self):
        """ return all elements as iterable key, value """
        return self.zip()

    def __getitem__(self,idx):
        """
            s[i] returns the ith element of the structure
            s[:4] returns a structure with the four first fields
            s[[1,3]] returns the second and fourth elements
        """
        if isinstance(idx,int):
            if idx<len(self):
                return self.getattr(self.keys()[idx])
            raise IndexError(f"the {self._ftype} index should be comprised between 0 and {len(self)-1}")
        elif isinstance(idx,slice):
            out = struct.fromkeysvalues(self.keys()[idx], self.values()[idx])
            if isinstance(self,paramauto):
                return paramauto(**out)
            elif isinstance(self,param):
                return param(**out)
            else:
                return out
        elif isinstance(idx,(list,tuple)):
            k,v= self.keys(), self.values()
            nk = len(k)
            s = param() if isinstance(self,param) else struct()
            for i in idx:
                if isinstance(i,int):
                    if -nk <= i < nk:  # Allow standard Python negative indexing
                        i = i % nk  # Convert negative index to positive equivalent
                        s.setattr(k[i],v[i])
                    else:
                        raise IndexError(f"idx must contain integers in range [-{nk}, {nk-1}], not {i}")
                elif isinstance(i,str):
                    if i in self:
                        s.setattr(i, self.getattr(i))
                    else:
                        raise KeyError((f'idx "{idx}" is not a valid key'))
                else:
                    TypeError("idx must contain only integers or strings")
            return s
        elif isinstance(idx,str):
            return self.getattr(idx)
        else:
            raise TypeError("The index must be an integer or a slice and not a %s" % type(idx).__name__)

    def __setitem__(self,idx,value):
        """ set the ith element of the structure  """
        if isinstance(idx,int):
            if idx<len(self):
                self.setattr(self.keys()[idx], value)
            else:
                raise IndexError(f"the {self._ftype} index should be comprised between 0 and {len(self)-1}")
        elif isinstance(idx,slice):
            k = self.keys()[idx]
            if len(value)<=1:
                for i in range(len(k)): self.setattr(k[i], value)
            elif len(k) == len(value):
                for i in range(len(k)): self.setattr(k[i], value[i])
            else:
                raise IndexError("the number of values (%d) does not match the number of elements in the slice (%d)" \
                       % (len(value),len(idx)))
        elif isinstance(idx,(list,tuple)):
            if len(value)<=1:
                for i in range(len(idx)): self[idx[i]]=value
            elif len(idx) == len(value):
                for i in range(len(idx)): self[idx[i]]=value[i]
            else:
                raise IndexError("the number of values (%d) does not match the number of indices (%d)" \
                                 % (len(value),len(idx)))

    def __len__(self):
        """ return the number of fields """
        # len() is used by struct() and its iterator
        return len(self.keys())

    def __iter__(self):
        """ struct iterator """
        # note that in the original object _iter_ is a static property not in dup
        dup = duplicate(self)
        dup._iter_ = 0
        return dup

    def __next__(self):
        """ increment iterator """
        self._iter_ += 1
        if self._iter_<=len(self):
            return self[self._iter_-1]
        self._iter_ = 0
        raise StopIteration(f"Maximum {self._ftype} iteration reached {len(self)}")

    def __add__(self, s, sortdefinitions=False, raiseerror=True, silentmode=True):
        """
        Add two structure objects, with precedence as follows:

          paramauto > param > struct

        In c = a + b, if b has a higher precedence than a then c will be of b's class,
        otherwise it will be of a's class.

        The new instance is created by copying the fields from the left-hand operand (a)
        and then updating with the fields from the right-hand operand (b).

        If self or s is of class paramauto, the current state of _needs_sorting is propagated
        but not forced to be true.

        """
        if not isinstance(s, struct):
            raise TypeError(f"the second operand must be {self._type}")

        # Define a helper to assign a precedence value.
        def get_precedence(obj):
            if isinstance(obj, paramauto):
                return 2
            elif isinstance(obj, param):
                return 1
            elif isinstance(obj, struct):
                return 0
            else:
                return 0  # fallback for unknown derivations

        # current classes
        leftprecedence = get_precedence(self)
        rightprecedence = get_precedence(s)
        # Determine which class to use for the duplicate.
        # If s (b) has a higher precedence than self (a), use s's class; otherwise, use self's.
        hi_class = self.__class__ if leftprecedence >= rightprecedence else s.__class__
        # Create a new instance of the chosen class by copying self's fields.
        dup = hi_class(**self)
        # Update with the fields from s.
        dup.update(**s)
        if sortdefinitions: # defer sorting by preserving the state of _needs_sorting
            if leftprecedence < rightprecedence == 2: # left is promoted
                dup._needs_sorting = s._needs_sorting
            elif rightprecedence < leftprecedence == 2: # right is promoted
                dup._needs_sorting = self._needs_sorting
            elif leftprecedence == rightprecedence == 2: # left and right are equivalent
                dup._needs_sorting = self._needs_sorting or s._needs_sorting
            # dup.sortdefinitions(raiseerror=raiseerror, silentmode=silentmode)
        return dup

    def __iadd__(self,s,sortdefinitions=False,raiseerror=False, silentmode=True):
        """ iadd a structure
            set sortdefintions=True to sort definitions (to maintain executability)
        """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        self.update(**s)
        if sortdefinitions:
            self._needs_sorting = True
            # self.sortdefinitions(raiseerror=raiseerror,silentmode=silentmode)
        return self

    def __sub__(self,s):
        """ sub a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        dup = duplicate(self)
        listofkeys = dup.keys()
        for k in s.keys():
            if k in listofkeys:
                delattr(dup,k)
        return dup

    def __isub__(self,s):
        """ isub a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        listofkeys = self.keys()
        for k in s.keys():
            if k in listofkeys:
                delattr(self,k)
        return self

    def dispmax(self,content):
        """ optimize display """
        strcontent = str(content)
        if len(strcontent)>self._maxdisplay:
            nchar = round(self._maxdisplay/2)
            return strcontent[:nchar]+" [...] "+strcontent[-nchar:]
        else:
            return content

    def __repr__(self):
        """ display method """
        if self.__dict__=={}:
            print(f"empty {self._fulltype} ({self._type} object) with no {self._type}s")
            return f"empty {self._fulltype}"
        else:
            numfmt = f".{self._precision}g"
            tmp = self.eval() if self._evalfeature else []
            keylengths = [len(key) for key in self.__dict__]
            width = max(10,max(keylengths)+2)
            fmt = "%%%ss:" % width
            fmteval = fmt[:-1]+"="
            fmtcls =  fmt[:-1]+":"
            line = ( fmt % ('-'*(width-2)) ) + ( '-'*(min(40,width*5)) )
            print(line)
            for key,value in self.__dict__.items():
                if key not in self._excludedattr:
                    if isinstance(value,_numeric_types):
                        # old code (removed on 2025-01-18)
                        # if isinstance(value,pstr):
                        #     print(fmt % key,'p"'+self.dispmax(value)+'"')
                        # if isinstance(value,str) and value=="":
                        #     print(fmt % key,'""')
                        # else:
                        #     print(fmt % key,self.dispmax(value))
                        if isinstance(value,np.ndarray):
                            print(fmt % key, struct.format_array(value,numfmt=numfmt))
                        else:
                            print(fmt % key,self.dispmax(value))
                    elif isinstance(value,struct):
                        print(fmt % key,self.dispmax(value.__str__()))
                    elif isinstance(value,(type,dict)):
                        print(fmt % key,self.dispmax(str(value)))
                    else:
                        print(fmt % key,type(value))
                        print(fmtcls % "",self.dispmax(str(value)))
                    if self._evalfeature:
                        if isinstance(self,paramauto):
                            try:
                                if isinstance(value,pstr):
                                    print(fmteval % "",'p"'+self.dispmax(tmp.getattr(key))+'"')
                                elif isinstance(value,str):
                                    if value == "":
                                        print(fmteval % "",self.dispmax("<empty string>"))
                                    else:
                                        print(fmteval % "",self.dispmax(tmp.getattr(key)))
                            except Exception as err:
                                print(fmteval % "",err.message, err.args)
                        else:
                            if isinstance(value,pstr):
                                print(fmteval % "",'p"'+self.dispmax(tmp.getattr(key))+'"')
                            elif isinstance(value,str):
                                if value == "":
                                    print(fmteval % "",self.dispmax("<empty string>"))
                                else:
                                    calcvalue =tmp.getattr(key)
                                    if isinstance(calcvalue, str) and "error" in calcvalue.lower():
                                        print(fmteval % "",calcvalue)
                                    else:
                                        if isinstance(calcvalue,np.ndarray):
                                            print(fmteval % "", struct.format_array(calcvalue,numfmt=numfmt))
                                        else:
                                            print(fmteval % "",self.dispmax(calcvalue))
                            elif isinstance(value,list):
                                calcvalue =tmp.getattr(key)
                                print(fmteval % "",self.dispmax(str(calcvalue)))
            print(line)
            return f"{self._fulltype} ({self._type} object) with {len(self)} {self._ftype}s"

    def disp(self):
        """ display method """
        self.__repr__()

    def __str__(self):
        return f"{self._fulltype} ({self._type} object) with {len(self)} {self._ftype}s"

    @property
    def isempty(self):
        """ isempty is set to True for an empty structure """
        return len(self)==0

    def clear(self):
        """ clear() delete all fields while preserving the original class """
        for k in self.keys(): delattr(self,k)

    def format(self, s, escape=False, raiseerror=True):
        """
            Format a string with fields using {field} as placeholders.
            Handles expressions like ${variable1}.

            Args:
                s (str): The input string to format.
                escape (bool): If True, prevents replacing '${' with '{'.
                raiseerror (bool): If True, raises errors for missing fields.

            Note:
                NumPy vectors and matrices are converted into their text representation (default behavior)
                If expressions such ${var[1,2]} are used an error is expected, the original content will be used instead

            Returns:
                str: The formatted string.
        """
        tmp = self.np2str()
        if raiseerror:
            try:
                if escape:
                    try: # we try to evaluate with all np objects converted in to strings (default)
                        return s.format_map(AttrErrorDict(tmp.__dict__))
                    except: # if an error occurs, we use the orginal content
                        return s.format_map(AttrErrorDict(self.__dict__))
                else:
                    try: # we try to evaluate with all np objects converted in to strings (default)
                        return s.replace("${", "{").format_map(AttrErrorDict(tmp.__dict__))
                    except: # if an error occurs, we use the orginal content
                        return s.replace("${", "{").format_map(AttrErrorDict(self.__dict__))
            except AttributeError as attr_err:
                # Handle AttributeError for expressions with operators
                s_ = s.replace("{", "${")
                if self._debug:
                    print(f"WARNING: the {self._ftype} {attr_err} is undefined in '{s_}'")
                return s_  # Revert to using '${' for unresolved expressions
            except IndexError as idx_err:
                s_ = s.replace("{", "${")
                if self._debug:
                    print(f"Index Error {idx_err} in '{s_}'")
                raise IndexError from idx_err
            except Exception as other_err:
                s_ = s.replace("{", "${")
                raise RuntimeError from other_err
        else:
            if escape:
                try: # we try to evaluate with all np objects converted in to strings (default)
                    return s.format_map(AttrErrorDict(tmp.__dict__))
                except: # if an error occurs, we use the orginal content
                    return s.format_map(AttrErrorDict(self.__dict__))
            else:
                try: # we try to evaluate with all np objects converted in to strings (default)
                    return s.replace("${", "{").format_map(AttrErrorDict(tmp.__dict__))
                except:  # if an error occurs, we use the orginal content
                    return s.replace("${", "{").format_map(AttrErrorDict(self.__dict__))

    def format_legacy(self,s,escape=False,raiseerror=True):
        """
            format a string with field (use {field} as placeholders)
                s.replace(string), s.replace(string,escape=True)
                where:
                    s is a struct object
                    string is a string with possibly ${variable1}
                    escape is a flag to prevent ${} replaced by {}
        """
        if raiseerror:
            try:
                if escape:
                    return s.format(**self.__dict__)
                else:
                    return s.replace("${","{").format(**self.__dict__)
            except KeyError as kerr:
                s_ = s.replace("{","${")
                print(f"WARNING: the {self._ftype} {kerr} is undefined in '{s_}'")
                return s_ # instead of s (we put back $) - OV 2023/01/27
            except Exception as othererr:
                s_ = s.replace("{","${")
                raise RuntimeError from othererr
        else:
            if escape:
                return s.format(**self.__dict__)
            else:
                return s.replace("${","{").format(**self.__dict__)

    def fromkeys(self,keys):
        """ returns a structure from keys """
        return self+struct(**dict.fromkeys(keys,None))

    @staticmethod
    def scan(s):
        """ scan(string) scan a string for variables """
        if not isinstance(s,str): raise TypeError("scan() requires a string")
        tmp = struct()
        #return tmp.fromkeys(set(re.findall(r"\$\{(.*?)\}",s)))
        found = re.findall(r"\$\{(.*?)\}",s);
        uniq = []
        for x in found:
            if x not in uniq: uniq.append(x)
        return tmp.fromkeys(uniq)

    @staticmethod
    def isstrexpression(s):
        """ isstrexpression(string) returns true if s contains an expression  """
        if not isinstance(s,str): raise TypeError("s must a string")
        return re.search(r"\$\{.*?\}",s) is not None

    @property
    def isexpression(self):
        """ same structure with True if it is an expression """
        s = param() if isinstance(self,param) else struct()
        for k,v in self.items():
            if isinstance(v,str):
                s.setattr(k,struct.isstrexpression(v))
            else:
                s.setattr(k,False)
        return s

    @staticmethod
    def isstrdefined(s,ref):
        """ isstrdefined(string,ref) returns true if it is defined in ref  """
        if not isinstance(s,str): raise TypeError("s must a string")
        if not isinstance(ref,struct): raise TypeError("ref must be a structure")
        if struct.isstrexpression(s):
            k = struct.scan(s).keys()
            allfound,i,nk = True,0,len(k)
            while (i<nk) and allfound:
                allfound = k[i] in ref
                i += 1
            return allfound
        else:
            return False


    def isdefined(self,ref=None):
        """ isdefined(ref) returns true if it is defined in ref """
        s = param() if isinstance(self,param) else struct()
        k,v,isexpr = self.keys(), self.values(), self.isexpression.values()
        nk = len(k)
        if ref is None:
            for i in range(nk):
                if isexpr[i]:
                    s.setattr(k[i],struct.isstrdefined(v[i],self[:i]))
                else:
                    s.setattr(k[i],True)
        else:
            if not isinstance(ref,struct): raise TypeError("ref must be a structure")
            for i in range(nk):
                if isexpr[i]:
                    s.setattr(k[i],struct.isstrdefined(v[i],ref))
                else:
                    s.setattr(k[i],True)
        return s


    def sortdefinitions(self,raiseerror=True,silentmode=False):
        """ sortdefintions sorts all definitions
            so that they can be executed as param().
            If any inconsistency is found, an error message is generated.

            Flags = default values
                raiseerror=True show erros of True
                silentmode=False no warning if True
        """
        find = lambda xlist: [i for i, x in enumerate(xlist) if x]
        findnot = lambda xlist: [i for i, x in enumerate(xlist) if not x]
        k,v,isexpr =  self.keys(), self.values(), self.isexpression.values()
        istatic = findnot(isexpr)
        idynamic = find(isexpr)
        static = struct.fromkeysvalues(
            [ k[i] for i in istatic ],
            [ v[i] for i in istatic ],
            makeparam = False)
        dynamic = struct.fromkeysvalues(
            [ k[i] for i in idynamic ],
            [ v[i] for i in idynamic ],
            makeparam=False)
        current = static # make static the current structure
        nmissing, anychange, errorfound = len(dynamic), False, False
        while nmissing:
            itst, found = 0, False
            while itst<nmissing and not found:
                teststruct = current + dynamic[[itst]] # add the test field
                found = all(list(teststruct.isdefined()))
                ifound = itst
                itst += 1
            if found:
                current = teststruct # we accept the new field
                dynamic[ifound] = []
                nmissing -= 1
                anychange = True
            else:
                if raiseerror:
                    raise KeyError('unable to interpret %d/%d expressions in "%ss"' % \
                                   (nmissing,len(self),self._ftype))
                else:
                    if (not errorfound) and (not silentmode):
                        print('WARNING: unable to interpret %d/%d expressions in "%ss"' % \
                              (nmissing,len(self),self._ftype))
                    current = teststruct # we accept the new field (even if it cannot be interpreted)
                    dynamic[ifound] = []
                    nmissing -= 1
                    errorfound = True
        if anychange:
            self.clear() # reset all fields and assign them in the proper order
            k,v = current.keys(), current.values()
            for i in range(len(k)):
                self.setattr(k[i],v[i])


    def generator(self, printout=False):
        """
        Generate Python code of the equivalent structure.

        This method converts the current structure (an instance of `param`, `paramauto`, or `struct`)
        into Python code that, when executed, recreates an equivalent structure. The generated code is
        formatted with one field per line.

        By default (when `printout` is False), the generated code is returned as a raw string that starts
        directly with, for example, `param(` (or `paramauto(` or `struct(`), with no "X = " prefix or leading
        newline. When `printout` is True, the generated code is printed to standard output and includes a prefix
        "X = " to indicate the variable name.

        Parameters:
            printout (bool): If True, the generated code is printed to standard output with the "X = " prefix.
                             If False (default), the code is returned as a raw string starting with, e.g.,
                             `param(`.

        Returns:
            str: The generated Python code representing the structure (regardless of whether it was printed).
        """
        nk = len(self)
        tmp = self.np2str()
        # Compute the field format based on the maximum key length (with a minimum width of 10)
        fmt = "%%%ss =" % max(10, max([len(k) for k in self.keys()]) + 2)
        # Determine the appropriate class string for the current instance.
        if isinstance(self, param):
            classstr = "param"
        elif 'paramauto' in globals() and isinstance(self, paramauto):
            classstr = "paramauto"
        else:
            classstr = "struct"

        lines = []
        if nk == 0:
            # For an empty structure.
            if printout:
                lines.append(f"X = {classstr}()")
            else:
                lines.append(f"{classstr}()")
        else:
            # Header: include "X = " only if printing.
            if printout:
                header = f"X = {classstr}("
            else:
                header = f"{classstr}("
            lines.append(header)
            # Iterate over keys to generate each field line.
            for i, k in enumerate(self.keys()):
                v = getattr(self, k)
                if isinstance(v, np.ndarray):
                    vtmp = getattr(tmp, k)
                    field = fmt % k + " " + vtmp
                elif isinstance(v, (int, float)) or v is None:
                    field = fmt % k + " " + str(v)
                elif isinstance(v, str):
                    field = fmt % k + " " + f'"{v}"'
                elif isinstance(v, (list, tuple, dict)):
                    field = fmt % k + " " + str(v)
                else:
                    field = fmt % k + " " + "/* unsupported type */"
                # Append a comma after each field except the last one.
                if i < nk - 1:
                    field += ","
                lines.append(field)
            # Create a closing line that aligns the closing parenthesis.
            closing_line = fmt[:-1] % ")"
            lines.append(closing_line)
        result = "\n".join(lines)
        if printout:
            print(result)
            return None
        return result


    # copy and deep copy methpds for the class
    def __copy__(self):
        """ copy method """
        cls = self.__class__
        copie = cls.__new__(cls)
        copie.__dict__.update(self.__dict__)
        return copie

    def __deepcopy__(self, memo):
        """ deep copy method """
        cls = self.__class__
        copie = cls.__new__(cls)
        memo[id(self)] = copie
        for k, v in self.__dict__.items():
            setattr(copie, k, duplicatedeep(v, memo))
        return copie


    # write a file
    def write(self, file, overwrite=True, mkdir=False):
        """
            write the equivalent structure (not recursive for nested struct)
                write(filename, overwrite=True, mkdir=False)

            Parameters:
            - file: The file path to write to.
            - overwrite: Whether to overwrite the file if it exists (default: True).
            - mkdir: Whether to create the directory if it doesn't exist (default: False).
        """
        # Create a Path object for the file to handle cross-platform paths
        file_path = Path(file).resolve()

        # Check if the directory exists or if mkdir is set to True, create it
        if mkdir:
            file_path.parent.mkdir(parents=True, exist_ok=True)
        elif not file_path.parent.exists():
            raise FileNotFoundError(f"The directory {file_path.parent} does not exist.")
        # If overwrite is False and the file already exists, raise an exception
        if not overwrite and file_path.exists():
            raise FileExistsError(f"The file {file_path} already exists, and overwrite is set to False.")
        # Convert to static if needed
        if isinstance(p,(param,paramauto)):
            tmp = self.tostatic()
        else:
            tmp = self
        # Open and write to the file using the resolved path
        with file_path.open(mode="w", encoding='utf-8') as f:
            print(f"# {self._fulltype} with {len(self)} {self._ftype}s\n", file=f)
            for k, v in tmp.items():
                if v is None:
                    print(k, "=None", file=f, sep="")
                elif isinstance(v, (int, float)):
                    print(k, "=", v, file=f, sep="")
                elif isinstance(v, str):
                    print(k, '="', v, '"', file=f, sep="")
                else:
                    print(k, "=", str(v), file=f, sep="")


    # read a file
    @staticmethod
    def read(file):
        """
            read the equivalent structure
                read(filename)

            Parameters:
            - file: The file path to read from.
        """
        # Create a Path object for the file to handle cross-platform paths
        file_path = Path(file).resolve()
        # Check if the parent directory exists, otherwise raise an error
        if not file_path.parent.exists():
            raise FileNotFoundError(f"The directory {file_path.parent} does not exist.")
        # If the file does not exist, raise an exception
        if not file_path.exists():
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        # Open and read the file
        with file_path.open(mode="r", encoding="utf-8") as f:
            s = struct()  # Assuming struct is defined elsewhere
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip()
                expr = line.split(sep="=")
                if len(line) > 0 and line[0] != "#" and len(expr) > 0:
                    lhs = expr[0]
                    rhs = "".join(expr[1:]).strip()
                    if len(rhs) == 0 or rhs == "None":
                        v = None
                    else:
                        v = eval(rhs)
                    s.setattr(lhs, v)
        return s

    # argcheck
    def check(self,default):
        """
        populate fields from a default structure
            check(defaultstruct)
            missing field, None and [] values are replaced by default ones

            Note: a.check(b) is equivalent to b+a except for [] and None values
        """
        if not isinstance(default,struct):
            raise TypeError("the first argument must be a structure")
        for f in default.keys():
            ref = default.getattr(f)
            if f not in self:
                self.setattr(f, ref)
            else:
                current = self.getattr(f)
                if ((current is None)  or (current==[])) and \
                    ((ref is not None) and (ref!=[])):
                        self.setattr(f, ref)


    # update values based on key:value
    def update(self, **kwargs):
        """
        Update multiple fields at once, while protecting certain attributes.

        Parameters:
        -----------
        **kwargs : dict
            The fields to update and their new values.

        Protected attributes defined in _excludedattr are not updated.

        Usage:
        ------
        s.update(a=10, b=[1, 2, 3], new_field="new_value")
        """
        protected_attributes = getattr(self, '_excludedattr', ())
        for key, value in kwargs.items():
            if key in protected_attributes:
                print(f"Warning: Cannot update protected attribute '{key}'")
            else:
                self.setattr(key, value)


    # override () for subindexing structure with key names
    def __call__(self, *keys):
        """
        Extract a sub-structure based on the specified keys,
        keeping the same class type.

        Parameters:
        -----------
        *keys : str
            The keys for the fields to include in the sub-structure.

        Returns:
        --------
        struct
            A new instance of the same class as the original, containing
            only the specified keys.

        Usage:
        ------
        sub_struct = s('key1', 'key2', ...)
        """
        # Create a new instance of the same class
        sub_struct = self.__class__()

        # Get the full type and field type for error messages
        fulltype = getattr(self, '_fulltype', 'structure')
        ftype = getattr(self, '_ftype', 'field')

        # Add only the specified keys to the new sub-structure
        for key in keys:
            if key in self:
                sub_struct.setattr(key, self.getattr(key))
            else:
                raise KeyError(f"{fulltype} does not contain the {ftype} '{key}'.")

        return sub_struct


    def __delattr__(self, key):
        """ Delete an instance attribute if it exists and is not a class or excluded attribute. """
        if key in self._excludedattr:
            raise AttributeError(f"Cannot delete excluded attribute '{key}'")
        elif key in self.__class__.__dict__:  # Check if it's a class attribute
            raise AttributeError(f"Cannot delete class attribute '{key}'")
        elif key in self.__dict__:  # Delete only if in instance's __dict__
            del self.__dict__[key]
        else:
            raise AttributeError(f"{self._type} has no attribute '{key}'")


    # A la Matlab display method of vectors, matrices and ND-arrays
    @staticmethod
    def format_array(value,numfmt=".4g"):
        """
        Format NumPy array for display with distinctions for scalars, row/column vectors, and ND arrays.
        Recursively formats multi-dimensional arrays without introducing unwanted commas.

        Args:
            value (np.ndarray): The NumPy array to format.
            numfmt: numeric format to be used for the string conversion (default=".4g")

        Returns:
            str: A formatted string representation of the array.
        """
        dtype_str = {
            np.float64: "double",
            np.float32: "single",
            np.int32: "int32",
            np.int64: "int64",
            np.complex64: "complex single",
            np.complex128: "complex double",
        }.get(value.dtype.type, str(value.dtype))  # Default to dtype name if not in the map

        max_display = 10  # Maximum number of elements to display

        def format_recursive(arr):
            """
            Recursively formats the array based on its dimensions.

            Args:
                arr (np.ndarray): The array or sub-array to format.

            Returns:
                str: Formatted string of the array.
            """
            if arr.ndim == 0:
                return f"{arr.item()}"

            if arr.ndim == 1:
                if len(arr) <= max_display:
                    return "[" + " ".join(f"{v:{numfmt}}" for v in arr) + "]"
                else:
                    return f"[{len(arr)} elements]"

            if arr.ndim == 2:
                if arr.shape[1] == 1:
                    # Column vector
                    if arr.shape[0] <= max_display:
                        return "[" + " ".join(f"{v[0]:{numfmt}}" for v in arr) + "]T"
                    else:
                        return f"[{arr.shape[0]}×1 vector]"
                elif arr.shape[0] == 1:
                    # Row vector
                    if arr.shape[1] <= max_display:
                        return "[" + " ".join(f"{v:{numfmt}}" for v in arr[0]) + "]"
                    else:
                        return f"[1×{arr.shape[1]} vector]"
                else:
                    # General matrix
                    return f"[{arr.shape[0]}×{arr.shape[1]} matrix]"

            # For higher dimensions
            shape_str = "×".join(map(str, arr.shape))
            if arr.size <= max_display:
                # Show full content
                if arr.ndim > 2:
                    # Represent multi-dimensional arrays with nested brackets
                    return "[" + " ".join(format_recursive(subarr) for subarr in arr) + f"] ({shape_str} {dtype_str})"
            return f"[{shape_str} array ({dtype_str})]"

        if value.size == 0:
            return "[]"

        if value.ndim == 0 or value.size == 1:
            return f"{value.item()} ({dtype_str})"

        if value.ndim == 1 or value.ndim == 2:
            # Use existing logic for vectors and matrices
            if value.ndim == 1:
                if len(value) <= max_display:
                    formatted = "[" + " ".join(f"{v:{numfmt}}" for v in value) + f"] ({dtype_str})"
                else:
                    formatted = f"[{len(value)}×1 {dtype_str}]"
            elif value.ndim == 2:
                rows, cols = value.shape
                if cols == 1:  # Column vector
                    if rows <= max_display:
                        formatted = "[" + " ".join(f"{v[0]:{numfmt}}" for v in value) + f"]T ({dtype_str})"
                    else:
                        formatted = f"[{rows}×1 {dtype_str}]"
                elif rows == 1:  # Row vector
                    if cols <= max_display:
                        formatted = "[" + " ".join(f"{v:{numfmt}}" for v in value[0]) + f"] ({dtype_str})"
                    else:
                        formatted = f"[1×{cols} {dtype_str}]"
                else:  # General matrix
                    formatted = f"[{rows}×{cols} {dtype_str}]"
            return formatted

        # For higher-dimensional arrays
        if value.size <= max_display:
            formatted = format_recursive(value)
        else:
            shape_str = "×".join(map(str, value.shape))
            formatted = f"[{shape_str} array ({dtype_str})]"

        return formatted


    # convert all NumPy entries to "nestable" expressions
    def np2str(self):
        """ Convert all NumPy entries of s into their string representations, handling both lists and dictionaries. """
        out = struct()
        def format_numpy_result(value):
            """
            Converts a NumPy array or scalar into a string representation:
            - Scalars and single-element arrays (any number of dimensions) are returned as scalars without brackets.
            - Arrays with more than one element are formatted with proper nesting and commas to make them valid `np.array()` inputs.
            - If the value is a list or dict, the conversion is applied recursively.
            - Non-ndarray inputs that are not list/dict are returned without modification.

            Args:
                value (np.ndarray, scalar, list, dict, or other): The value to format.

            Returns:
                str, list, dict, or original type: A properly formatted string for NumPy arrays/scalars,
                a recursively converted list/dict, or the original value.
            """
            if isinstance(value, dict):
                # Recursively process each key in the dictionary.
                new_dict = {}
                for k, v in value.items():
                    new_dict[k] = format_numpy_result(v)
                return new_dict
            elif isinstance(value, list):
                # Recursively process each element in the list.
                return [format_numpy_result(x) for x in value]
            elif isinstance(value, tuple):
                return tuple(format_numpy_result(x) for x in value)
            elif isinstance(value, struct):
                return value.npstr()
            elif np.isscalar(value):
                # For scalars: if numeric, use str() to avoid extra quotes.
                if isinstance(value, (int, float, complex, str)) or value is None:
                    return value
                else:
                    return repr(value)
            elif isinstance(value, np.ndarray):
                # Check if the array has exactly one element.
                if value.size == 1:
                    # Extract the scalar value.
                    return repr(value.item())
                # Convert the array to a nested list.
                nested_list = value.tolist()
                # Recursively format the nested list into a valid string.
                def list_to_string(lst):
                    if isinstance(lst, list):
                        return "[" + ",".join(list_to_string(item) for item in lst) + "]"
                    else:
                        return repr(lst)
                return list_to_string(nested_list)
            else:
                # Return the input unmodified if not a NumPy array, list, dict, or scalar.
                return str(value) # str() preferred over repr() for concision
        # Process all entries in self.
        for key, value in self.items():
            out.setattr(key, format_numpy_result(value))
        return out

    # minimal replacement of placeholders by numbers or their string representations
    def numrepl(self, text):
        r"""
        Replace all placeholders of the form ${key} in the given text by the corresponding
        numeric value from the instance fields, under the following conditions:

        1. 'key' must be a valid field in self (i.e., if key in self).
        2. The value corresponding to 'key' is either:
             - an int,
             - a float, or
             - a string that represents a valid number (e.g., "1" or "1.0").

        Only when these conditions are met, the placeholder is substituted.
        The conversion preserves the original type: if the stored value is int, then the
        substitution will be done as an integer (e.g., 1 and not 1.0). Otherwise, if it is
        a float then it will be substituted as a float.

        Any placeholder for which the above conditions are not met remains unchanged.

        Placeholders are recognized by the pattern "${<key>}" where <key> is captured as all
        text until the next "}" (optionally allowing whitespace inside the braces).
        """
        # Pattern: match "${", then optional whitespace, capture all characters until "}",
        # then optional whitespace, then "}".
        placeholder_pattern = re.compile(r"\$\{\s*([^}]+?)\s*\}")

        def replace_match(match):
            key = match.group(1)
            # Check if the key exists in self.
            if key in self:
                value = self[key]
                # If the value is already numeric, substitute directly.
                if isinstance(value, (int, float)):
                    return str(value)
                # If the value is a string, try to interpret it as a numeric value.
                elif isinstance(value, str):
                    s = value.strip()
                    # Check if s is a valid integer representation.
                    if re.fullmatch(r"[+-]?\d+", s):
                        try:
                            num = int(s)
                            return str(num)
                        except ValueError:
                            # Should not occur because the regex already matched.
                            return match.group(0)
                    # Check if s is a valid float representation (including scientific notation).
                    elif re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?", s):
                        try:
                            num = float(s)
                            return str(num)
                        except ValueError:
                            return match.group(0)
            # If key not in self or value is not numeric (or numeric string), leave placeholder intact.
            return match.group(0)

        # Replace all placeholders in the text using the replacer function.
        return placeholder_pattern.sub(replace_match, text)

    # import method
    def importfrom(self, s, nonempty=True, replacedefaultvar=True):
        """
        Import values from 's' into self according to the following rules:

        - Only fields that already exist in self are considered.
        - If s is a dictionary, it is converted to a struct via struct(**s).
        - If the current value of a field in self is empty (None, "", [] or ()),
          then that field is updated from s.
        - If nonempty is True (default), then only non-empty values from s are imported.
        - If replacedefaultvar is True (default), then if a field in self exactly equals
          "${key}" (with key being the field name), it is replaced by the corresponding
          value from s if it is empty.
        - Raises a TypeError if s is not a dict or struct (i.e. if it doesn’t support keys()).
        """
        # If s is a dictionary, convert it to a struct instance.
        if isinstance(s, dict):
            s = struct(**s)
        elif not hasattr(s, "keys"):
            raise TypeError(f"s must be a struct or a dictionary not a type {type(s).__name__}")

        for key in self.keys():
            if key in s:
                s_value = getattr(s, key)
                current_value = getattr(self, key)
                if is_empty(current_value) or (replacedefaultvar and current_value == "${" + key + "}"):
                    if nonempty:
                        if not is_empty(s_value):
                            setattr(self, key, s_value)
                    else:
                        setattr(self, key, s_value)

    # importfrom with copy
    def __lshift__(self, other):
        """
        Allows the syntax:

            s = s1 << s2

        where a new instance is created as a copy of s1 (preserving its type, whether
        struct, param, or paramauto) and then updated with the values from s2 using
        importfrom.
        """
        # Create a new instance preserving the type of self.
        new_instance = type(self)(**{k: getattr(self, k) for k in self.keys()})
        # Import values from other (s2) into the new instance.
        new_instance.importfrom(other)
        return new_instance

    # returns only valid keys
    def validkeys(self, list_of_keys):
        """
        Validate and return the subset of keys from the provided list that are valid in the instance.

        Parameters:
        -----------
        list_of_keys : list
            A list of keys (as strings) to check against the instance’s attributes.

        Returns:
        --------
        list
            A list of keys from list_of_keys that are valid (i.e., exist as attributes in the instance).

        Raises:
        -------
        TypeError
            If list_of_keys is not a list or if any element in list_of_keys is not a string.

        Example:
        --------
        >>> s = struct()
        >>> s.foo = 42
        >>> s.bar = "hello"
        >>> valid = s.validkeys(["foo", "bar", "baz"])
        >>> print(valid)   # Output: ['foo', 'bar'] assuming 'baz' is not defined in s
        """
        # Check that list_of_keys is a list
        if not isinstance(list_of_keys, list):
            raise TypeError("list_of_keys must be a list")

        # Check that every entry in the list is a string
        for key in list_of_keys:
            if not isinstance(key, str):
                raise TypeError("Each key in list_of_keys must be a string")

        # Assuming valid keys are those present in the instance's __dict__
        return [key for key in list_of_keys if key in self]

# %% param class with scripting and evaluation capabilities
class param(struct):
    """
    Class: `param`
    ==============

    A class derived from `struct` that introduces dynamic evaluation of field values.
    The `param` class acts as a container for evaluated parameters, allowing expressions
    to depend on other fields. It supports advanced evaluation, sorting of dependencies,
    and text formatting.

    ---

    ### Features
    - Inherits all functionalities of `struct`.
    - Supports dynamic evaluation of field expressions.
    - Automatically resolves dependencies between fields.
    - Includes utility methods for text formatting and evaluation.

    ### Shorthands for `p=param(...)`
    - `s = p.eval()` returns the full evaluated structure
    - `p.getval("field")` returns the evaluation for the field "field"
    - `s = p()` returns the full evaluated structure as `p.eval()`
    - `s = p("field1","field2"...)` returns the evaluated substructure for fields "field1", "field2"

    ---

    ### Examples

    #### Basic Usage with Evaluation
    ```python
    s = param(a=1, b=2, c='${a} + ${b} # evaluate me if you can', d="$this is a string", e="1000 # this is my number")
    s.eval()
    # Output:
    # --------
    #      a: 1
    #      b: 2
    #      c: ${a} + ${b} # evaluate me if you can (= 3)
    #      d: $this is a string (= this is a string)
    #      e: 1000 # this is my number (= 1000)
    # --------

    s.a = 10
    s.eval()
    # Output:
    # --------
    #      a: 10
    #      b: 2
    #      c: ${a} + ${b} # evaluate me if you can (= 12)
    #      d: $this is a string (= this is a string)
    #      e: 1000 # this is my number (= 1000)
    # --------
    ```

    #### Handling Text Parameters
    ```python
    s = param()
    s.mypath = "$/this/folder"
    s.myfile = "$file"
    s.myext = "$ext"
    s.fullfile = "$${mypath}/${myfile}.${myext}"
    s.eval()
    # Output:
    # --------
    #    mypath: $/this/folder (= /this/folder)
    #    myfile: $file (= file)
    #     myext: $ext (= ext)
    #  fullfile: $${mypath}/${myfile}.${myext} (= /this/folder/file.ext)
    # --------
    ```

    ---

    ### Text Evaluation and Formatting

    #### Evaluate Strings
    ```python
    s = param(a=1, b=2)
    result = s.eval("this is a string with ${a} and ${b}")
    print(result)  # "this is a string with 1 and 2"
    ```

    #### Prevent Evaluation
    ```python
    definitions = param(a=1, b="${a}*10+${a}", c="\\${a}+10", d='\\${myparam}')
    text = definitions.formateval("this is my text ${a}, ${b}, \\${myvar}=${c}+${d}")
    print(text)  # "this is my text 1, 11, \\${myvar}=\\${a}+10+${myparam}"
    ```

    ---

    ### Advanced Usage

    #### Rearranging and Sorting Definitions
    ```python
    s = param(
        a=1,
        f="${e}/3",
        e="${a}*${c}",
        c="${a}+${b}",
        b=2,
        d="${c}*2"
    )
    s.sortdefinitions()
    s.eval()
    # Output:
    # --------
    #      a: 1
    #      b: 2
    #      c: ${a} + ${b} (= 3)
    #      d: ${c} * 2 (= 6)
    #      e: ${a} * ${c} (= 3)
    #      f: ${e} / 3 (= 1.0)
    # --------
    ```

    #### Internal Evaluation and Recursion with !
    ```python
    p=param()
    p.a = [0,1,2]
    p.b = '![1,2,"test","${a[1]}"]'
    p
    # Output:
    #  -------------:----------------------------------------
    #          a: [0, 1, 2]
    #          b: ![1,2,"test","${a[1]}"]
    #           = [1, 2, 'test', '1']
    #  -------------:----------------------------------------
    # Out: parameter list (param object) with 2 definitions
    ```

    #### Error Handling
    ```python
    p = param(b="${a}+1", c="${a}+${d}", a=1)
    p.disp()
    # Output:
    # --------
    #      b: ${a} + 1 (= 2)
    #      c: ${a} + ${d} (= < undef definition "${d}" >)
    #      a: 1
    # --------
    ```

    Sorting unresolved definitions raises errors unless explicitly suppressed:
    ```python
    p.sortdefinitions(raiseerror=False)
    # WARNING: unable to interpret 1/3 expressions in "definitions"
    ```

    ---

    ### Utility Methods
    | Method                 | Description                                             |
    |------------------------|---------------------------------------------------------|
    | `eval()`               | Evaluate all field expressions.                         |
    | `formateval(string)`   | Format and evaluate a string with field placeholders.   |
    | `protect(string)`      | Escape variable placeholders in a string.               |
    | `sortdefinitions()`    | Sort definitions to resolve dependencies.               |
    | `escape(string)`       | Protect escaped variables in a string.                  |
    | `safe_fstring(string)` | evaluate safely complex mathemical expressions.         |

    ---

    ### Overloaded Methods and Operators
    #### Supported Operators
    - `+`: Concatenation of two parameter lists, sorting definitions.
    - `-`: Subtraction of fields.
    - `len()`: Number of fields.
    - `in`: Check for field existence.

    ---

    ### Notes
    - The `paramauto` class simplifies handling of partial definitions and inherits from `param`.
    - Use `paramauto` when definitions need to be stacked irrespective of execution order.
    """

    # override
    _type = "param"
    _fulltype = "parameter list"
    _ftype = "definition"
    _evalfeature = True    # This class can be evaluated with .eval()
    _returnerror = True    # This class returns an error in the evaluation string (added on 2024-09-06)


    # magic constructor
    def __init__(self,_protection=False,_evaluation=True,
                 sortdefinitions=False,debug=False,**kwargs):
        """ constructor """
        super().__init__(debug=debug,**kwargs)
        self._protection = _protection
        self._evaluation = _evaluation
        self._needs_sorting = False # defers sorting
        if sortdefinitions: self.sortdefinitions()

    # escape definitions if needed
    @staticmethod
    def escape(s):
        """
            escape \\${} as ${{}} --> keep variable names
            convert ${} as {} --> prepare Python replacement

            Examples:
                escape("\\${a}")
                returns ('${{a}}', True)

                escape("  \\${abc} ${a} \\${bc}")
                returns ('  ${{abc}} {a} ${{bc}}', True)

                escape("${a}")
                Out[94]: ('{a}', False)

                escape("${tata}")
                returns ('{tata}', False)

        """
        if not isinstance(s,str):
            raise TypeError(f'the argument must be string not {type(s)}')
        se, start, found = "", 0, True
        while found:
            pos0 = s.find(r"\${",start)
            found = pos0>=0
            if found:
                pos1 = s.find("}",pos0)
                found = pos1>=0
                if found:
                    se += s[start:pos0].replace("${","{")+"${{"+s[pos0+3:pos1]+"}}"
                    start=pos1+1
        result = se+s[start:].replace("${","{")
        if isinstance(s,pstr): result = pstr(result)
        return result,start>0

    # protect variables in a string
    def protect(self,s=""):
        """ protect $variable as ${variable} """
        if isinstance(s,str):
            t = s.replace(r"\$","££") # && is a placeholder
            escape = t!=s
            for k in self.keyssorted():
                t = t.replace("$"+k,"${"+k+"}")
            if escape: t = t.replace("££",r"\$")
            if isinstance(s,pstr): t = pstr(t)
            return t, escape
        raise TypeError(f'the argument must be string not {type(s)}')


    # lines starting with # (hash) are interpreted as comments
    # ${variable} or {variable} are substituted by variable.value
    # any line starting with $ is assumed to be a string (no interpretation)
    # ^ is accepted in formula(replaced by **))
    def eval(self,s="",protection=False):
        """
            Eval method for structure such as MS.alias

                s = p.eval() or s = p.eval(string)

                where :
                    p is a param object
                    s is a structure with evaluated fields
                    string is only used to determine whether definitions have been forgotten

        """
        # handle deferred sorting
        if self._needs_sorting:
            self.sortdefinitions(raiseerror=False, silentmode=True)
        # the argument s is only used by formateval() for error management
        tmp = struct(debug=self._debug)
        # evaluator without context
        evaluator_nocontext = SafeEvaluator() # for global evaluation without context

        # main string evaluator
        def evalstr(value,key=""):
            # replace ${variable} (Bash, Lammps syntax) by {variable} (Python syntax)
            # use \${variable} to prevent replacement (espace with \)
            # Protect variables if required
            ispstr = isinstance(value,pstr)
            valuesafe = pstr.eval(value,ispstr=ispstr) # value.strip()
            if valuesafe=="${"+key+"}": # circular reference (it cannot be evaluated)
                return valuesafe
            if protection or self._protection:
                valuesafe, escape0 = self.protect(valuesafe)
            else:
                escape0 = False
            # replace ${var} by {var} once basic substitutions have been applied
            valuesafe_priorescape = tmp.numrepl(valuesafe) # minimal substitution
            valuesafe, escape = param.escape(valuesafe_priorescape)
            escape = escape or escape0
            # replace "^" (Matlab, Lammps exponent) by "**" (Python syntax)
            valuesafe = pstr.eval(valuesafe.replace("^","**"),ispstr=ispstr)
            # Remove all content after #
            # if the first character is '#', it is not comment (e.g. MarkDown titles)
            poscomment = valuesafe.find("#")
            if poscomment>0: valuesafe = valuesafe[0:poscomment].strip()
            # Matrix shorthand replacement
            # $[[1,2,${a}]]+$[[10,20,30]] --> np.array([[1,2,${a}]])+np.array([[10,20,30]])
            valuesafe = param.replace_matrix_shorthand(valuesafe)
            # Literal string starts with $ (no interpretation), ! (evaluation)
            if not self._evaluation:
                return pstr.eval(tmp.format(valuesafe,escape),ispstr=ispstr)
            elif valuesafe.startswith("!"): # <---------- FORECED LITERAL EVALUATION (error messages are returned)
                try:
                    #vtmp = ast.literal_eval(valuesafe[1:])
                    evaluator = SafeEvaluator(tmp)
                    vtmp = evaluate_with_placeholders(valuesafe[1:],evaluator,evaluator_nocontext)
                    if isinstance(vtmp,list):
                        for i,item in enumerate(vtmp):
                            if isinstance(item,str) and not is_literal_string(item):
                                try:
                                    vtmp[i] = tmp.format(item, raiseerror=False) # in case substitions/interpolations are needed
                                    try:
                                        vtmp[i] = evaluator_nocontext.evaluate(vtmp[i]) # full evaluation without context
                                    except Exception as othererr:
                                        if self._debug:
                                            print(f"DEBUG {key}: Error evaluating: {vtmp[i]}\n< {othererr} >")
                                except Exception as ve:
                                    vtmp[i] = f"Error in <{item}>: {ve.__class__.__name__} - {str(ve)}"
                    return vtmp
                except (SyntaxError, ValueError) as e:
                    return f"Error: {e.__class__.__name__} - {str(e)}"
            elif valuesafe.startswith("$") and not escape:
                return tmp.format(valuesafe[1:].lstrip()) # discard $
            elif valuesafe.startswith("%"):
                return tmp.format(valuesafe[1:].lstrip()) # discard %
            else: # string empty or which can be evaluated
                if valuesafe=="":
                    return valuesafe # empty content
                else:
                    if isinstance(value,pstr): # keep path
                        return pstr.topath(tmp.format(valuesafe,escape=escape))
                    elif escape:  # partial evaluation
                        return tmp.format(valuesafe,escape=True)
                    else: # full evaluation (if it fails the last string content is returned) <---------- FULL EVALUTION will be tried
                        try:
                            resstr = tmp.format(valuesafe,raiseerror=False)
                        except (KeyError,NameError) as nameerr:
                            try: # nested indexing (guess)
                                resstr = param.safe_fstring(
                                param.replace_matrix_shorthand(valuesafe_priorescape),tmp)
                                evaluator = SafeEvaluator(tmp)
                                reseval = evaluate_with_placeholders(resstr,
                                          evaluator,evaluator_nocontext,raiseerror=True)
                            except Exception as othererr:
                                if self._returnerror: # added on 2024-09-06
                                    strnameerr = str(nameerr).replace("'","")
                                    if self._debug:
                                        print(f'Key Error for "{key}" < {othererr} >')
                                    return '< undef %s "${%s}" >' % (self._ftype,strnameerr)
                                else:
                                    return value #we keep the original value
                            else:
                                return reseval
                        except SyntaxError as commonerr:
                            return "Syntax Error < %s >" % commonerr
                        except TypeError  as commonerr:
                            try: # nested indexing (guess)
                                resstr = param.safe_fstring(
                                param.replace_matrix_shorthand(valuesafe_priorescape),tmp)
                                evaluator = SafeEvaluator(tmp)
                                reseval = evaluate_with_placeholders(resstr,
                                          evaluator,evaluator_nocontext,raiseerror=True)
                            except Exception as othererr:
                                if self._debug:
                                    print(f'Type Error for "{key}" < {othererr} >')
                                return "Type Error < %s >" % commonerr
                            else:
                                return reseval
                        except (IndexError,AttributeError):
                            try:
                                resstr = param.safe_fstring(
                                    param.replace_matrix_shorthand(valuesafe_priorescape),tmp)
                            except Exception as fstrerr:
                                return "Index Error < %s >" % fstrerr
                            else:
                                try:
                                    # reseval = eval(resstr)
                                    # reseval = ast.literal_eval(resstr)
                                    # Use SafeEvaluator to evaluate the final expression
                                    evaluator = SafeEvaluator(tmp)
                                    reseval = evaluator.evaluate(resstr)
                                except Exception as othererr:
                                    #tmp.setattr(key,"Mathematical Error around/in ${}: < %s >" % othererr)
                                    if self._debug:
                                        print(f"DEBUG {key}: Error evaluating: {resstr}\n< {othererr} >")
                                    return resstr
                                else:
                                    return reseval
                        except ValueError as valerr: # forced evaluation within ${}
                            try:
                                evaluator = SafeEvaluator(tmp)
                                reseval = evaluate_with_placeholders(valuesafe_priorescape,evaluator,evaluator_nocontext,raiseerror=True)
                            except SyntaxError as synerror:
                                if self._debug:
                                    print(f"DEBUG {key}: Error evaluating: {valuesafe_priorescape}\n< {synerror} >")
                                return evaluate_with_placeholders(valuesafe_priorescape,evaluator,evaluator_nocontext,raiseerror=False)
                            except Exception as othererr:
                                if self._debug:
                                    print(f"DEBUG {key}: Error evaluating: {valuesafe_priorescape}\n< {othererr} >")
                                return "Error in ${}: < %s >" % valerr
                            else:
                                return reseval

                        except Exception as othererr:
                            return "Error in ${}: < %s >" % othererr
                        else:
                            try:
                                # reseval = eval(resstr)
                                evaluator = SafeEvaluator(tmp)
                                reseval = evaluate_with_placeholders(resstr,evaluator,evaluator_nocontext)
                            except Exception as othererr:
                                #tmp.setattr(key,"Eval Error < %s >" % othererr)
                                if self._debug:
                                    print(f"DEBUG {key}: Error evaluating: {resstr}\n< {othererr} >")
                                return resstr.replace("\n",",") # \n replaced by ,
                            else:
                                return reseval

        # evalstr() refactored for error management
        def safe_evalstr(x,key=""):
            xeval = evalstr(x,key)
            if isinstance(xeval,str):
                try:
                    evaluator = SafeEvaluator(tmp)
                    return evaluate_with_placeholders(xeval,evaluator,evaluator_nocontext)
                except Exception as e:
                    if self._debug:
                        print(f"DEBUG {key}: Error evaluating '{x}': {e}")
                    return xeval  # default fallback value
            else:
                return xeval

        # Evaluate all DEFINITIONS
        for key,value in self.items():
            # strings are assumed to be expressions on one single line
            if isinstance(value,str):
                tmp.setattr(key,evalstr(value,key))
            elif isinstance(value,_numeric_types): # already a number
                if isinstance(value,list):
                    valuelist = [safe_evalstr(x,key) if isinstance(x,str) else x for x in value]
                    tmp.setattr(key,valuelist)
                else:
                    tmp.setattr(key, value) # store the value with the key
            elif isinstance(value, dict):
                # For dictionaries, evaluate each entry using its own key (sub_key)
                new_dict = {}
                for sub_key, sub_value in value.items():
                    if isinstance(sub_value, str):
                        new_dict[sub_key] = safe_evalstr(sub_value,key)
                    elif isinstance(sub_value, list):
                        # If an entry is a list, apply safe_evalstr to each string element within it
                        new_dict[sub_key] = [safe_evalstr(x, sub_key) if isinstance(x, str) else x
                                             for x in sub_value]
                    else:
                        new_dict[sub_key] = sub_value
                tmp.setattr(key, new_dict)
            else: # unsupported types
                if s.find("{"+key+"}")>=0:
                    print(f'*** WARNING ***\n\tIn the {self._ftype}:"\n{s}\n"')
                else:
                    print(f'unable to interpret the "{key}" of type {type(value)}')
        return tmp

    # formateval obeys to following rules
    # lines starting with # (hash) are interpreted as comments
    def formateval(self,s,protection=False,fullevaluation=True):
        """
            format method with evaluation feature

                txt = p.formateval("this my text with ${variable1}, ${variable2} ")

                where:
                    p is a param object

                Example:
                    definitions = param(a=1,b="${a}",c="\\${a}")
                    text = definitions.formateval("this my text ${a}, ${b}, ${c}")
                    print(text)

        """
        tmp = self.eval(s,protection=protection)
        evaluator = SafeEvaluator(tmp) # used when fullevaluation=True
        evaluator_nocontext = SafeEvaluator() # for global evaluation without context
        # Do all replacements in s (keep comments)
        if len(tmp)==0:
            return s
        else:
            ispstr = isinstance(s,pstr)
            ssafe, escape = param.escape(s)
            slines = ssafe.split("\n")
            slines_priorescape = s.split("\n")
            for i in range(len(slines)):
                poscomment = slines[i].find("#")
                if poscomment>=0:
                    while (poscomment>0) and (slines[i][poscomment-1]==" "):
                        poscomment -= 1
                    comment = slines[i][poscomment:len(slines[i])]
                    slines[i]  = slines[i][0:poscomment]
                else:
                    comment = ""
                # Protect variables if required
                if protection or self._protection:
                    slines[i], escape2 = self.protect(slines[i])
                # conversion
                if ispstr:
                    slines[i] = pstr.eval(tmp.format(slines[i],escape=escape),ispstr=ispstr)
                else:
                    if fullevaluation:
                        try:
                            resstr =tmp.format(slines[i],escape=escape)
                        except:
                            resstr = param.safe_fstring(slines[i],tmp,varprefix="")
                        try:
                            #reseval = evaluator.evaluate(resstr)
                            reseval = evaluate_with_placeholders(slines_priorescape[i],evaluator,evaluator_nocontext,raiseerror=True)
                            slines[i] = str(reseval)+" "+comment if comment else str(reseval)
                        except:
                            slines[i] = resstr + comment
                    else:
                        slines[i] = tmp.format(slines[i],escape=escape)+comment
                # convert starting % into # to authorize replacement in comments
                if len(slines[i])>0:
                    if slines[i][0] == "%": slines[i]="#"+slines[i][1:]
            return "\n".join(slines)


    # return the value instead of formula
    def getval(self,key):
        """ returns the evaluated value """
        s = self.eval()
        return getattr(s,key)

    # override () for subindexing structure with key names
    def __call__(self, *keys):
        """
        Extract an evaluated sub-structure based on the specified keys,
        keeping the same class type.

        Parameters:
        -----------
        *keys : str
            The keys for the fields to include in the sub-structure.

        Returns:
        --------
        struct
            An evaluated instance of class struct, containing
            only the specified keys with evaluated values.

        Usage:
        ------
        sub_struct = p('key1', 'key2', ...)
        """
        s = self.eval()
        if keys:
            return s(*keys)
        else:
            return s

    # returns the equivalent structure evaluated
    def tostruct(self,protection=False):
        """
            generate the evaluated structure
                tostruct(protection=False)
        """
        return self.eval(protection=protection)

    # returns the equivalent paramauto instance
    def toparamauto(self):
        """
            convert a param instance into a paramauto instance
                toparamauto()
        """
        return paramauto(**self)

    # returns the equivalent structure evaluated
    def tostatic(self):
        """ convert dynamic a param() object to a static struct() object.
            note: no interpretation
            note: use tostruct() to interpret them and convert it to struct
            note: tostatic().struct2param() makes it reversible
        """
        return struct.fromkeysvalues(self.keys(),self.values(),makeparam=False)

    # Matlab vector/list conversion
    @staticmethod
    def expand_ranges(text,numfmt=".4g"):
        """
        Expands MATLAB-style ranges in a string.

        Args:
            text: The input string containing ranges.
            numfmt: numeric format to be used for the string conversion (default=".4g")

        Returns:
            The string with ranges expanded, or the original string if no valid ranges are found
            or if expansion leads to more than 100 elements. Returns an error message if the input
            format is invalid.
        """
        def expand_range(match):
            try:
                parts = match.group(1).split(':')
                if len(parts) == 2:
                    start, stop = map(float, parts)
                    step = 1.0
                elif len(parts) == 3:
                    start, step, stop = map(float, parts)
                else:
                    return match.group(0)  # Return original if format is invalid
                if step == 0:
                    return "Error: <Step cannot be zero.>"
                if (stop - start) / step > 1e6:
                    return "Error: <Range is too large.>"
                if step > 0:
                    num_elements = int(np.floor((stop - start)/step)+1)
                else:
                     num_elements = int(np.floor((start - stop)/-step)+1)
                if num_elements > 100:
                    return match.group(0)  # Return original if too many elements
                expanded_range = np.arange(start, stop + np.sign(step)*1e-9, step) #adding a small number to include the stop in case of integer steps
                return '[' + ','.join(f'{x:{numfmt}}' for x in expanded_range) + ']'
            except ValueError:
                return "Error: <Invalid range format.>"
        pattern = r'(\b(?:-?\d+(?:\.\d*)?|-?\.\d+)(?::(?:-?\d+(?:\.\d*)?|-?\.\d+)){1,2})\b'
        expanded_text = re.sub(pattern, expand_range, text)
        #check for errors generated by the function
        if "Error:" in expanded_text:
            return expanded_text
        return expanded_text

    # Matlab syntax conversion
    @staticmethod
    def convert_matlab_like_arrays(text):
        """
        Converts Matlab-like array syntax (including hybrid notations) into
        a NumPy-esque list syntax in multiple passes.

        Steps:
          1) Convert 2D Matlab arrays (containing semicolons) into Python-like nested lists.
          2) Convert bracketed row vectors (no semicolons or nested brackets) into double-bracket format.
          3) Replace spaces with commas under specific conditions and remove duplicates.

        Args:
            text (str): Input string that may contain Matlab-like arrays.

        Returns:
            str: Transformed text with arrays converted to a Python/NumPy-like syntax.

        Examples:

            examples = [
                "[1, 2  ${var1}          ; 4, 5     ${var2}]",
                "[1,2,3]",
                "[1 2 ,  3]",
                "[1;2; 3]",
                "[[-0.5, 0.5;-0.5, 0.5],[ -0.5,  0.5; -0.5,  0.5]]",
                "[[1,2;3,4],[5,6; 7,8]]",
                "[1, 2, 3; 4, 5, 6]",  # Hybrid
                "[[Already, in, Python]]",  # Already Python-like?
                "Not an array"
            ]

            for ex in examples:
                converted = param.convert_matlab_like_arrays(ex)
                print(f"Matlab: {ex}\nNumPy : {converted}\n")

        """
        # --------------------------------------------------------------------------
        # Step 1: Detect innermost [ ... ; ... ] blocks and convert them
        # --------------------------------------------------------------------------
        def convert_matrices_with_semicolons(txt):
            """
            Repeatedly find the innermost bracket pair that contains a semicolon
            and convert it to a Python-style nested list, row by row.
            """
            # Pattern to find innermost [ ... ; ... ] without nested brackets
            pattern = r'\[[^\[\]]*?;[^\[\]]*?\]'
            while True:
                match = re.search(pattern, txt)
                if not match:
                    break  # No more [ ... ; ... ] blocks to convert
                inner_block = match.group(0)
                # Remove the outer brackets
                inner_content = inner_block[1:-1].strip()
                # Split into rows by semicolon
                rows = [row.strip() for row in inner_content.split(';')]
                converted_rows = []
                for row in rows:
                    # Replace multiple spaces with a single space
                    row_clean = re.sub(r'\s+', ' ', row)
                    # Split row by commas or spaces
                    row_elems = re.split(r'[,\s]+', row_clean)
                    row_elems = [elem for elem in row_elems if elem]  # Remove empty strings
                    # Join elements with commas and encapsulate in brackets
                    converted_rows.append("[" + ",".join(row_elems) + "]")
                # Join the row lists and encapsulate them in brackets
                replacement = "[" + ",".join(converted_rows) + "]"
                # Replace the original Matlab matrix with the Python list
                txt = txt[:match.start()] + replacement + txt[match.end():]
            return txt
        # --------------------------------------------------------------------------
        # Step 2: Convert row vectors without semicolons or nested brackets
        #         into double-bracket format, e.g. [1,2,3] -> [[1,2,3]]
        # --------------------------------------------------------------------------
        def convert_row_vectors(txt):
            """
            Convert [1,2,3] or [1 2 3] into [[1,2,3]] if the bracket does not contain
            semicolons, nor nested brackets. We do this iteratively, skipping any
            bracket blocks that don't qualify, rather than stopping.
            """
            # We only want bracket blocks that are NOT preceded by '[' or ','
            # do not contain semicolons or nested brackets
            # and are not followed by ']' or ','
            pattern = r"(?<!\[)(?<!,)\([^();]*\)(?!\s*\])(?!\s*\,)"
            startpos = 0
            while True:
                match = re.search(pattern, txt[startpos:])
                if not match:
                    break  # No more bracket blocks to check
                # Compute absolute positions in txt
                mstart = startpos + match.start()
                mend   = startpos + match.end()
                block  = txt[mstart:mend]
                # we need to be sure that [ ] are present around the block even if separated by spaces
                if mstart == 0 or mend == len(txt) - 1:
                    break
                if not (re.match(r"\[\s*$", txt[:mstart]) and re.match(r"^\s*\]", txt[mend+1:])):
                    break
                # Double-check that this bracket does not contain semicolons or nested brackets
                # If it does, we skip it (just advance the search) to avoid messing up matrices.
                # That is, we do not transform it into double brackets.
                if ';' in block or '[' in block[1:-1] or ']' in block[1:-1]:
                    # Move beyond this match and keep searching
                    startpos = mend
                    continue
                # It's a pure row vector (no semicolons, no nested brackets)
                new_block = "[" + block + "]"  # e.g. [1,2,3] -> [[1,2,3]]
                txt = txt[:mstart] + new_block + txt[mend:]
                # Update search position to avoid re-matching inside the newly inserted text
                startpos = mstart + len(new_block)
            return txt
        # --------------------------------------------------------------------------
        # Step 3: Replace spaces with commas under specific conditions and clean up
        # --------------------------------------------------------------------------
        def replace_spaces_safely(txt):
            """
            Replace spaces not preceded by ',', '[', or whitespace and followed by digit or '$' with commas.
            Then remove multiple consecutive commas and trailing commas before closing brackets.
            """
            # 1) Replace spaces with commas if not preceded by [,\s
            #    and followed by digit or $
            txt = re.sub(r'(?<![,\[\s])\s+(?=[\d\$])', ',', txt)
            # 2) Remove multiple consecutive commas
            txt = re.sub(r',+', ',', txt)
            # 3) Remove trailing commas before closing brackets
            txt = re.sub(r',+\]', ']', txt)
            return txt

        def replace_spaces_with_commas(txt):
            """
            Replaces spaces with commas only when they're within array shorthands and not preceded by a comma, opening bracket, or whitespace,
            and are followed by a digit or '$'. Also collapses multiple commas into one and strips leading/trailing commas.

            Parameters:
            ----------
            txt : str
                The text to process.

            Returns:
            -------
            str
                The processed text with appropriate commas.
            """
            # Replace spaces not preceded by ',', '[', or whitespace and followed by digit or '$' with commas
            txt = re.sub(r'(?<![,\[\s])\s+(?=[\d\$])', ',', txt)
            # Remove multiple consecutive commas
            txt = re.sub(r',+', ',', txt)
            # Strip leading and trailing commas
            txt = txt.strip(',')
            # Replace residual multiple consecutive spaces with a single space
            return re.sub(r'\s+', ' ', txt)

        # --------------------------------------------------------------------------
        # Apply Step 1: Convert matrices with semicolons
        # --------------------------------------------------------------------------
        text_cv = convert_matrices_with_semicolons(text)
        # --------------------------------------------------------------------------
        # Apply Step 2: Convert row vectors (no semicolons/nested brackets)
        # --------------------------------------------------------------------------
        text_cv = convert_row_vectors(text_cv)
        # --------------------------------------------------------------------------
        # Apply Step 3: Replace spaces with commas and clean up
        # --------------------------------------------------------------------------
        if text_cv != text:
            return replace_spaces_with_commas(text_cv) # old method: replace_spaces_safely(text_cv)
        else:
            return text



    @classmethod
    def replace_matrix_shorthand(cls,valuesafe):
        """
        Transforms custom shorthand notations for NumPy arrays within a string into valid NumPy array constructors.
        Supports up to 4-dimensional arrays and handles variable references.

        **Shorthand Patterns:**
        - **1D**: `$[1 2 3]` → `np.atleast_2d(np.array([1,2,3]))`
        - **2D**: `$[[1 2],[3 4]]` → `np.array([[1,2],[3,4]])`
        - **3D**: `$[[[1 2],[3 4]],[[5 6],[7 8]]]` → `np.array([[[1,2],[3,4]],[[5,6],[7,8]]])`
        - **4D**: `$[[[[1 2]]]]` → `np.array([[[[1,2]]]])`
        - **Variable References**: `@{var}` → `np.atleast_2d(np.array(${var}))`

        **Parameters:**
        ----------
        valuesafe : str
            The input string containing shorthand notations for NumPy arrays and variable references.

        **Returns:**
        -------
        str
            The transformed string with shorthands replaced by valid NumPy array constructors.

        **Raises:**
        -------
        ValueError
            If there are unmatched brackets in any shorthand.

        **Examples:**
        --------
        >>> # 1D shorthand
        >>> s = "$[1 2 3]"
        >>> param.replace_matrix_shorthand(s)
        'np.atleast_2d(np.array([1,2,3]))'

        >>> # 2D shorthand with mixed spacing
        >>> s = "$[[1, 2], [3 4]]"
        >>> param.replace_matrix_shorthand(s)
        'np.array([[1,2],[3,4]])'

        >>> # 3D array with partial spacing
        >>> s = "$[[[1  2], [3 4]], [[5 6], [7 8]]]"
        >>> param.replace_matrix_shorthand(s)
        'np.array([[[1,2],[3,4]],[[5,6],[7,8]]])'

        >>> # 4D array
        >>> s = "$[[[[1 2]]]]"
        >>> param.replace_matrix_shorthand(s)
        'np.array([[[[1,2]]]])'

        >>> # Combined with variable references
        >>> s = "@{a} + $[[${b}, 2],[ 3  4]]"
        >>> param.replace_matrix_shorthand(s)
        'np.atleast_2d(np.array(${a})) + np.array([[${b},2],[3,4]])'

        >>> # Complex ND array with scaling
        >>> s = '$[[[-0.5, -0.5],[-0.5, -0.5]],[[ 0.5,  0.5],[ 0.5,  0.5]]]*0.001'
        >>> param.replace_matrix_shorthand(s)
        'np.array([[[-0.5,-0.5],[-0.5,-0.5]],[[0.5,0.5],[0.5,0.5]]])*0.001'
        """

        numfmt = f".{cls._precision}g"

        def replace_spaces_with_commas(txt):
            """
            Replaces spaces with commas only when they're not preceded by a comma, opening bracket, or whitespace,
            and are followed by a digit or '$'. Also collapses multiple commas into one and strips leading/trailing commas.

            Parameters:
            ----------
            txt : str
                The text to process.

            Returns:
            -------
            str
                The processed text with appropriate commas.
            """
            # Replace spaces not preceded by ',', '[', or whitespace and followed by digit or '$' with commas
            txt = re.sub(r'(?<![,\[\s])\s+(?=[\d\$])', ',', txt)
            # Remove multiple consecutive commas
            txt = re.sub(r',+', ',', txt)
            return txt.strip(',')

        def build_pass_list(string):
            """
            Determines which dimensions (1D..4D) appear in the string by searching for:
             - 4D: $[[[[
             - 3D: $[[[
             - 2D: $[[
             - 1D: $[

            Returns a sorted list in descending order, e.g., [4, 3, 2, 1].

            Parameters:
            ----------
            string : str
                The input string to scan.

            Returns:
            -------
            list
                A list of integers representing the dimensions found, sorted descending.
            """
            dims_found = set()
            if re.search(r'\$\[\[\[\[', string):
                dims_found.add(4)
            if re.search(r'\$\[\[\[', string):
                dims_found.add(3)
            if re.search(r'\$\[\[', string):
                dims_found.add(2)
            if re.search(r'\$\[', string):
                dims_found.add(1)
            return sorted(dims_found, reverse=True)

        # Step 0: convert eventual Matlab syntax for row and column vectors into NumPy syntax
        valuesafe = param.expand_ranges(valuesafe,numfmt)  # expands start:stop and start:step:stop syntax
        valuesafe = param.convert_matlab_like_arrays(valuesafe) # vectors and matrices conversion
        # Step 1: Handle @{var} -> np.atleast_2d(np.array(${var}))
        valuesafe = re.sub(r'@\{([^\{\}]+)\}', r'np.atleast_2d(np.array(${\1}))', valuesafe)
        # Step 2: Build pass list from largest dimension to smallest
        pass_list = build_pass_list(valuesafe)
        # Step 3: Define patterns and replacements for each dimension
        dimension_patterns = {
            4: (r'\$\[\[\[\[(.*?)\]\]\]\]', 'np.array([[[[{content}]]]])'),   # 4D
            3: (r'\$\[\[\[(.*?)\]\]\]', 'np.array([[[{content}]]])'),         # 3D
            2: (r'\$\[\[(.*?)\]\]', 'np.array([[{content}]])'),               # 2D
            1: (r'\$\[(.*?)\]', 'np.atleast_2d(np.array([{content}]))')       # 1D
        }
        # Step 4: Iterate over each dimension and perform replacements
        for dim in pass_list:
            pattern, replacement_fmt = dimension_patterns[dim]
            # Find all non-overlapping matches for the current dimension
            matches = list(re.finditer(pattern, valuesafe))
            for match in matches:
                full_match = match.group(0)       # Entire matched shorthand
                inner_content = match.group(1)    # Content inside the brackets
                # Replace spaces with commas as per rules
                processed_content = replace_spaces_with_commas(inner_content.strip())
                # Create the replacement string
                replacement = replacement_fmt.format(content=processed_content)
                # Replace the shorthand in the string
                valuesafe = valuesafe.replace(full_match, replacement)
        # Step 5: Verify that all shorthands have been replaced by checking for remaining '$['
        if re.search(r'\$\[', valuesafe):
            raise ValueError("Unmatched or improperly formatted brackets detected in the input string.")
        return valuesafe



    # Safe fstring
    @staticmethod
    def safe_fstring(template, context,varprefix="$"):
        """Safely evaluate expressions in ${} using SafeEvaluator."""
        evaluator = SafeEvaluator(context)
        # Process template string in combination with safe_fstring()
        # it is required to have an output compatible with eval()

        def process_template(valuesafe):
            """
            Processes the input string by:
            1. Stripping leading and trailing whitespace.
            2. Removing comments (any text after '#' unless '#' is the first character).
            3. Replacing '^' with '**'.
            4. Replacing '{' with '${' if '{' is not preceded by '$'. <-- not applied anymore (brings confusion)

            Args:
                valuesafe (str): The input string to process.

            Returns:
                str: The processed string.
            """
            # Step 1: Strip leading and trailing whitespace
            valuesafe = valuesafe.strip()
            # Step 2: Remove comments
            # This regex removes '#' and everything after it if '#' is not the first character
            # (?<!^) is a negative lookbehind that ensures '#' is not at the start of the string
            valuesafe = re.sub(r'(?<!^)\#.*', '', valuesafe)
            # Step 3: Replace '^' with '**'
            valuesafe = re.sub(r'\^', '**', valuesafe)
            # Step 4: Replace '{' with '${' if '{' is not preceded by '$'
            # (?<!\$)\{ matches '{' not preceded by '$'
            # valuesafe = re.sub(r'(?<!\$)\{', '${', valuesafe)
            # Optional: Strip again to remove any trailing whitespace left after removing comments
            valuesafe = valuesafe.strip()
            return valuesafe

        # Adjusted display for NumPy arrays
        def serialize_result(result):
            """
            Serialize the result into a string that can be evaluated in Python.
            Handles NumPy arrays by converting them to lists with commas.
            Handles other iterable types appropriately.
            """
            if isinstance(result, np.ndarray):
                return str(result.tolist())
            elif isinstance(result, (list, tuple, dict)):
                return str(result)
            else:
                return str(result)
        # Regular expression to find ${expr} patterns
        escaped_varprefix = re.escape(varprefix)
        pattern = re.compile(escaped_varprefix+r'\{([^{}]+)\}')
        def replacer(match):
            expr = match.group(1)
            try:
                result = evaluator.evaluate(expr)
                serialized = serialize_result(result)
                return serialized
            except Exception as e:
                return f"<Error: {e}>"
        return pattern.sub(replacer, process_template(template))

# %% str class for file and paths
# this class guarantees that paths are POSIX at any time


class pstr(str):
    """
    Class: `pstr`
    =============

    A specialized string class for handling paths and filenames, derived from `struct`.
    The `pstr` class ensures compatibility with POSIX-style paths and provides enhanced
    operations for path manipulation.

    ---

    ### Features
    - Maintains POSIX-style paths.
    - Automatically handles trailing slashes.
    - Supports path concatenation using `/`.
    - Converts seamlessly back to `str` for compatibility with string methods.
    - Includes additional utility methods for path evaluation and formatting.

    ---

    ### Examples

    #### Basic Usage
    ```python
    a = pstr("this/is/mypath//")
    b = pstr("mylocalfolder/myfile.ext")
    c = a / b
    print(c)  # this/is/mypath/mylocalfolder/myfile.ext
    ```

    #### Keeping Trailing Slashes
    ```python
    a = pstr("this/is/mypath//")
    print(a)  # this/is/mypath/
    ```

    ---

    ### Path Operations

    #### Path Concatenation
    Use the `/` operator to concatenate paths:
    ```python
    a = pstr("folder/subfolder")
    b = pstr("file.txt")
    c = a / b
    print(c)  # folder/subfolder/file.txt
    ```

    #### Path Evaluation
    Evaluate or convert paths while preserving the `pstr` type:
    ```python
    result = pstr.eval("some/path/afterreplacement", ispstr=True)
    print(result)  # some/path/afterreplacement
    ```

    ---

    ### Advanced Usage

    #### Using String Methods
    Methods like `replace()` convert `pstr` back to `str`. To retain the `pstr` type:
    ```python
    new_path = pstr.eval(a.replace("mypath", "newpath"), ispstr=True)
    print(new_path)  # this/is/newpath/
    ```

    #### Handling POSIX Paths
    The `pstr.topath()` method ensures the path remains POSIX-compliant:
    ```python
    path = pstr("C:\\Windows\\Path")
    posix_path = path.topath()
    print(posix_path)  # C:/Windows/Path
    ```

    ---

    ### Overloaded Operators

    #### Supported Operators
    - `/`: Concatenates two paths (`__truediv__`).
    - `+`: Concatenates strings as paths, resulting in a `pstr` object (`__add__`).
    - `+=`: Adds to an existing `pstr` object (`__iadd__`).

    ---

    ### Utility Methods

    | Method          | Description                                  |
    |------------------|----------------------------------------------|
    | `eval(value)`    | Evaluates the path or string for compatibility with `pstr`. |
    | `topath()`       | Returns the POSIX-compliant path.           |

    ---

    ### Notes
    - Use `pstr` for consistent and safe handling of file paths across different platforms.
    - Converts back to `str` when using non-`pstr` specific methods to ensure compatibility.
    """

    def __repr__(self):
        result = self.topath()
        if result[-1] != "/" and self[-1] == "/":
            result += "/"
        return result

    def topath(self):
        """ return a validated path """
        value = pstr(PurePath(self))
        if value[-1] != "/" and self [-1]=="/":
            value += "/"
        return value


    @staticmethod
    def eval(value,ispstr=False):
        """ evaluate the path of it os a path """
        if isinstance(value,pstr):
            return value.topath()
        elif isinstance(value,PurePath) or ispstr:
            return pstr(value).topath()
        else:
            return value

    def __truediv__(self,value):
        """ overload / """
        operand = pstr.eval(value)
        result = pstr(PurePath(self) / operand)
        if result[-1] != "/" and operand[-1] == "/":
            result += "/"
        return result

    def __add__(self,value):
        return pstr(str(self)+value)

    def __iadd__(self,value):
        return pstr(str(self)+value)


# %% class paramauto() which enforces sortdefinitions = True, raiseerror=False
class paramauto(param):
    """
    Class: `paramauto`
    ==================

    A subclass of `param` with enhanced handling for automatic sorting and evaluation
    of definitions. The `paramauto` class ensures that all fields are sorted to resolve
    dependencies, allowing seamless stacking of partially defined objects.

    ---

    ### Features
    - Inherits all functionalities of `param`.
    - Automatically sorts definitions for dependency resolution.
    - Simplifies handling of partial definitions in dynamic structures.
    - Supports safe concatenation of definitions.

    ---

    ### Examples

    #### Automatic Dependency Sorting
    Definitions are automatically sorted to resolve dependencies:
    ```python
    p = paramauto(a=1, b="${a}+1", c="${a}+${b}")
    p.disp()
    # Output:
    # --------
    #      a: 1
    #      b: ${a} + 1 (= 2)
    #      c: ${a} + ${b} (= 3)
    # --------
    ```

    #### Handling Missing Definitions
    Unresolved dependencies raise warnings but do not block execution:
    ```python
    p = paramauto(a=1, b="${a}+1", c="${a}+${d}")
    p.disp()
    # Output:
    # --------
    #      a: 1
    #      b: ${a} + 1 (= 2)
    #      c: ${a} + ${d} (= < undef definition "${d}" >)
    # --------
    ```

    ---

    ### Concatenation and Inheritance
    Concatenating `paramauto` objects resolves definitions:
    ```python
    p1 = paramauto(a=1, b="${a}+2")
    p2 = paramauto(c="${b}*3")
    p3 = p1 + p2
    p3.disp()
    # Output:
    # --------
    #      a: 1
    #      b: ${a} + 2 (= 3)
    #      c: ${b} * 3 (= 9)
    # --------
    ```

    ---

    ### Utility Methods

    | Method                | Description                                            |
    |-----------------------|--------------------------------------------------------|
    | `sortdefinitions()`   | Automatically sorts fields to resolve dependencies.    |
    | `eval()`              | Evaluate all fields, resolving dependencies.           |
    | `disp()`              | Display all fields with their resolved values.         |

    ---

    ### Overloaded Operators

    #### Supported Operators
    - `+`: Concatenates two `paramauto` objects, resolving dependencies.
    - `+=`: Updates the current object with another, resolving dependencies.
    - `len()`: Number of fields.
    - `in`: Check for field existence.

    ---

    ### Advanced Usage

    #### Partial Definitions
    The `paramauto` class simplifies handling of partially defined fields:
    ```python
    p = paramauto(a="${d}", b="${a}+1")
    p.disp()
    # Warning: Unable to resolve dependencies.
    # --------
    #      a: ${d} (= < undef definition "${d}" >)
    #      b: ${a} + 1 (= < undef definition "${d}" >)
    # --------

    p.d = 10
    p.disp()
    # Dependencies are resolved:
    # --------
    #      d: 10
    #      a: ${d} (= 10)
    #      b: ${a} + 1 (= 11)
    # --------
    ```

    ---

    ### Notes
    - The `paramauto` class is computationally more intensive than `param` due to automatic sorting.
    - It is ideal for managing dynamic systems with complex interdependencies.

    ### Examples
                    p = paramauto()
                    p.b = "${aa}"
                    p.disp()
                yields
                    WARNING: unable to interpret 1/1 expressions in "definitions"
                      -----------:----------------------------------------
                                b: ${aa}
                                 = < undef definition "${aa}" >
                      -----------:----------------------------------------
                      p.aa = 2
                      p.disp()
                yields
                    -----------:----------------------------------------
                             aa: 2
                              b: ${aa}
                               = 2
                    -----------:----------------------------------------
                    q = paramauto(c="${aa}+${b}")+p
                    q.disp()
                yields
                    -----------:----------------------------------------
                             aa: 2
                              b: ${aa}
                               = 2
                              c: ${aa}+${b}
                               = 4
                    -----------:----------------------------------------
                    q.aa = 30
                    q.disp()
                yields
                    -----------:----------------------------------------
                             aa: 30
                              b: ${aa}
                               = 30
                              c: ${aa}+${b}
                               = 60
                    -----------:----------------------------------------
                    q.aa = "${d}"
                    q.disp()
                yields multiple errors (recursion)
                WARNING: unable to interpret 3/3 expressions in "definitions"
                  -----------:----------------------------------------
                           aa: ${d}
                             = < undef definition "${d}" >
                            b: ${aa}
                             = Eval Error < invalid [...] (<string>, line 1) >
                            c: ${aa}+${b}
                             = Eval Error < invalid [...] (<string>, line 1) >
                  -----------:----------------------------------------
                    q.d = 100
                    q.disp()
                yields
                  -----------:----------------------------------------
                            d: 100
                           aa: ${d}
                             = 100
                            b: ${aa}
                             = 100
                            c: ${aa}+${b}
                             = 200
                  -----------:----------------------------------------


            Example:

                p = paramauto(b="${a}+1",c="${a}+${d}",a=1)
                p.disp()
            generates:
                WARNING: unable to interpret 1/3 expressions in "definitions"
                  -----------:----------------------------------------
                            a: 1
                            b: ${a}+1
                             = 2
                            c: ${a}+${d}
                             = < undef definition "${d}" >
                  -----------:----------------------------------------
            setting p.d
                p.d = 2
                p.disp()
            produces
                  -----------:----------------------------------------
                            a: 1
                            d: 2
                            b: ${a}+1
                             = 2
                            c: ${a}+${d}
                             = 3
                  -----------:----------------------------------------

    """

    def __add__(self,p):
        return super().__add__(p,sortdefinitions=True,raiseerror=False)
        self._needs_sorting = True

    def __iadd__(self,p):
        return super().__iadd__(p,sortdefinitions=True,raiseerror=False)
        self._needs_sorting = True

    def __repr__(self):
        self.sortdefinitions(raiseerror=False)
        #super(param,self).__repr__()
        super().__repr__()
        return str(self)

    def setattr(self,key,value):
        """ set field and value """
        if isinstance(value,list) and len(value)==0 and key in self:
            delattr(self, key)
        else:
            self.__dict__[key] = value
            self.__dict__["_needs_sorting"] = True

    def sortdefinitions(self, raiseerror=True, silentmode=True):
        if self._needs_sorting:
            super().sortdefinitions(raiseerror=raiseerror, silentmode=silentmode)
            self._needs_sorting = False

# %% DEBUG
# ===================================================
# main()
# ===================================================
# for debugging purposes (code called as a script)
# the code is called from here
# ===================================================
if __name__ == '__main__':
# =============================================================================
#     # very advanced
#     import os
#     from fitness.private.loadods import alias
#     local = "C:/Users/olivi/OneDrive/Data/Olivier/INRA/Etudiants & visiteurs/Steward Ouadi/python/test/output/"
#     odsfile = "fileid_conferences_FoodRisk.ods"
#     fullfodsfile = os.path.join(local,odsfile)
#     p = alias(fullfodsfile)
#     p.disp()
# =============================================================================
# new feature
    a = struct(a=1,b=2)
    a["b"]
# path example
    s0 = struct(a=pstr("/tmp/"),b=pstr("test////"),c=pstr("${a}/${b}"),d=pstr("${a}/${c}"),e=pstr("$c/$a"))
    s = struct.struct2param(s0,protection=True)
    s.disp()
    s.a/s.b
    str(pstr.topath(f"{s.a}/{s.b}"))
    s.eval()
    # escape example
    definitions = param(a=1,b="${a}*10+${a}",c=r"\${a}+10",d=r'\${myparam}')
    text = definitions.formateval(r"this my text ${a}, ${b}, \${myvar}=${c}+${d}")
    print(text)

    definitions = param(a=1,b="$a*10+$a",c=r"\$a+10",d=r'\$myparam')
    text = definitions.formateval(r"this my text $a, $b, \$myvar=$c+$d",protection=True)
    print(text)
    # assignment
    s = struct(a=1,b=2)
    s[1] = 3
    s.disp()
    # conversion
    s = {"a":1, "b":2}
    t=struct.dict2struct(s)
    t.disp()
    sback = t.struct2dict()
    sback.__repr__()
    # file definition
    p=struct.fromkeysvalues(["a","b","c","d"],[1,2,3]).struct2param()
    ptxt = p.protect("$c=$a+$b")
    definitions.write("../../tmp/test.txt")
    # populate/inherit fields
    default = struct(a=1,b="2",c=[1,2,3])
    tst = struct(a=10)
    tst.check(default)
    tst.disp()
    # multiple assigment
    a = struct(a=1,b=2,c=3,d=4)
    b = struct(a=10,b=20,c=30,d=40)
    a[:2] = b[1:3]
    a[:2] = b[(1,3)]
    # reorganize definitions to enable param.eval()
    s = param(
        a = 1,
        f = "${e}/3",
        e = "${a}*${c}",
        c = "${a}+${b}",
        b = 2,
        d = "${c}*2"
        )
    #s[0:2] = [1,2]
    s.isexpression
    struct.isstrdefined("${a}+${b}",s)
    s.isdefined()
    s.sortdefinitions()
    s.disp()
    p = param(b="${a}+1",c="${a}+${d}",a=1)
    p.disp()

# features 2025
    p=param()
    p.a = [0,1,2]
    p.b = '![1,2,"test","${a[1]}"]'
    p

# Mathematical expressions
    # Example: param.safe_fstring()
    # Sample context with a NumPy array
    context = param(
        f = np.array([
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16]
        ])
    )
    # Example expressions
    expressions = [
        "${a[1]}",                # Should return 0.2 (assuming 'a' is defined in context)
        "${b[0,1]} + ${a[0]}",    # Should return 1.2 (assuming 'b' and 'a' are defined)
        "${f[0:2,1]}"               # Should return the second column of 'f'
    ]
    # Assuming 'a' and 'b' are defined in the context
    context.update(
        a =[1.0, 0.2, 0.03, 0.004],
        b = np.array([[1, 0.2, 0.03, 0.004]])
    )
    for expr in expressions:
        result = param.safe_fstring(expr, context)
        print(f"Expression: {expr} => Result: {result}")

    # OUTPUT
    #   -------------:----------------------------------------
    # Expression: ${a[1]} => Result: 0.2
    # Expression: ${b[0,1]} + ${a[0]} => Result: 0.2 + 1.0
    # Expression: ${f[0:2,1]} => Result: [2, 6]
    #   -------------:----------------------------------------

    # Example with matrix operations
    p=param()
    p.a = [1.0, .2, .03, .004]
    p.b = np.array([p.a])
    p.c = p.a*2
    p.d = p.b*2
    p.e = p.b.T
    p.f = p.b.T@p.b # Matrix multiplication for (3x1) @ (1x3)
    p.g = "${a[1]}"
    p.h = "${b[0,1]} + ${a[0]}"
    p.i = "${f[0,1]}"
    p.j = "${f[:,1]}"
    p.k = "@{j}+1"  # note that "@{j}+1" and "${j}+1" do not have the same meaning
    p.l = "${b.T}"
    p.m = "${b.T @ b}"    # evaluate fully the matrix operation
    p.n = "${b.T} @ ${b}" # concatenate two string-results separated by @
    p.o ="the result is: ${b[0,1]} + ${a[0]}"
    p.p = "the value of a[0] is ${a[0]}"
    p.q = "1+1"
    print(repr(p))

    # OUTPUT
    #   -------------:----------------------------------------
    #               a: [1.0, 0.2, 0.03, 0.004]
    #               b: [1 0.2 0.03 0.004] (double)
    #               c: [1.0, 0.2, 0.03, 0.0 [...] 0, 0.2, 0.03, 0.004]
    #               d: [2 0.4 0.06 0.008] (double)
    #               e: [1 0.2 0.03 0.004]T (double)
    #               f: [4×4 double]
    #               g: ${a[1]}
    #                = 0.2
    #               h: ${b[0,1]} + ${a[0]}
    #                = 1.2
    #               i: ${f[0,1]}
    #                = 0.2
    #               j: ${f[:,1]}
    #                = [0.2, 0.040000000000 [...] 0001, 0.006, 0.0008]
    #               k: ${j}+1
    #                = [0.2, 0.040000000000 [...] 01, 0.006, 0.0008]+1
    #               l: ${b.T}
    #                = [[1.   ], [0.2  ], [0.03 ], [0.004]]
    #               m: ${b.T @ b}
    #                = [[1.0, 0.2, 0.03, 0. [...] , 0.00012, 1.6e-05]]
    #               n: ${b.T} @ ${b}
    #                = [[1.   ], [0.2  ], [ [...]  0.2   0.03  0.004]]
    #               o: the result is: ${b[0,1]} + ${a[0]}
    #                = the result is: 0.2 + 1.0
    #               p: the value of a[0] is ${a[0]}
    #                = the value of a[0] is 1.0
    #               q: 1+1
    #                = 2
    #   -------------:----------------------------------------
    # parameter list (param object) with 17 definitions


# Example with new NumPy shorthands
    p = param(debug=True);
    p.a = 1.0
    p.b = "10.0"
    p.c = "$[${a},2,3]*${b}" # Create a Numpy vector from an operation
    p.n = "$[0,0,1]"         # another one
    p.o1 = "@{n}"          # create a copy
    p.o2 = "$[${a},2,3]" # create a Numpy vector
    p.o3 = "@{o1} @ @{o2}.T" # multiplication between two vectots
    p.d = "@{n}.T @ $[[${a},2,3]]" # another one
    p.f = "($[${a},2,3]*${b}) @ np.array([[0,0,1]]).T" # another one using explicitly NumPy
    p.nT = "@{n}.T" # transpose of a vector/matrix
    p.m = "${n.T}*2" # this operation is illegal and will be kept as a string
    p.o = "@{n}.T*2" # this one is the correct one
    p.p = "$[[1,2],[3,4]]" # Create a 2D Numpy array
    p.q = "${p[1,1]}"   # index a 2D NumPy array
    p.r = "${p[:,1]}"   # this is a valid syntax to get the slice as a list
    p.s = "@{p}[:,1]+1" # use this syntax if you need apply an operation to the slice
    # more advanced
    p.V1 = "$[1.0,0.2,0.03]"
    p.V2 = "@{V1}+1"
    p.V3 = "@{V1}.T @ @{V2}"
    p.V4 = "np.diag(@{V3})"
    p.V5 = "np.linalg.eig(@{V3})"
    p.out = "the first eigenvalue is: ${V5.eigenvalues[0]}"
    print(repr(p))


    # OUTPUT with DEBUG MESSAGES
    # DEBUG m: Error in Evaluating: [[0]
    #  [0]
    #  [1]]*2
    # < Invalid index 1 for object of type int: 'int' object is not subscriptable >
    # DEBUG out: Error in Evaluating: the first eigenvalue is: 2.0
    # < invalid syntax (<unknown>, line 1) >
    #   -------------:----------------------------------------
    #               a: 1.0
    #               b: 10.0
    #                = 10.0
    #               c: $[${a},2,3]*${b}
    #                = [10 20 30] (double)
    #               n: $[0,0,1]
    #                = [0 0 1] (int64)
    #              o1: @{n}
    #                = [0 0 1] (int64)
    #              o2: $[${a},2,3]
    #                = [1 2 3] (double)
    #              o3: @{o1} @ @{o2}.T
    #                = 3.0 (double)
    #               d: @{n}.T @ $[[${a},2,3]]
    #                = [3×3 double]
    #               f: ($[${a},2,3]*${b}) @ [...] p.array([[0,0,1]]).T
    #                = 30.0 (double)
    #              nT: @{n}.T
    #                = [0 0 1]T (int64)
    #               m: ${n.T}*2
    #                = [[0], [0], [1]]*2
    #               o: @{n}.T*2
    #                = [0 0 2]T (int64)
    #               p: $[[1,2],[3,4]]
    #                = [2×2 int64]
    #               q: ${p[1,1]}
    #                = 4
    #               r: ${p[:,1]}
    #                = [2, 4]
    #               s: @{p}[:,1]+1
    #                = [3 5] (int64)
    #              V1: $[1.0,0.2,0.03]
    #                = [1 0.2 0.03] (double)
    #              V2: @{V1}+1
    #                = [2 1.2 1.03] (double)
    #              V3: @{V1}.T @ @{V2}
    #                = [3×3 double]
    #              V4: np.diag(@{V3})
    #                = [2 0.24 0.0309] (double)
    #              V5: np.linalg.eig(@{V3})
    #                = EigResult(eigenvalue [...] 332, -0.06057363]]))
    #             out: the first eigenvalue is: ${V3[0,0]}
    #                = the first eigenvalue is: 2.0
    #   -------------:----------------------------------------
    # parameter list (param object) with 22 definitions

# Advanced NumPy example

    p = param(debug=True)
    p.p = "$[[1, 2], [3, 4]]"      # Create a 2D NumPy array
    p.q = "${p[1, 1]}"             # Indexing: retrieves 4
    p.r = "@{p}[:,1] + 1"          # Add 1 to the second column
    p.s = "@{p}[:, 1].reshape(-1, 1) @ @{r}" # perform p(:,1)'*s in Matlab sense
    p.t = "np.linalg.eig(@{s})"
    p.w = "${t.eigenvalues[0]} + ${t.eigenvalues[1]}" # sum of eigen values
    p.x = "$[[0,${t.eigenvalues[0]}+${t.eigenvalues[1]}]]" # horizontal concat à la Matlab
    print(repr(p))

    # Output
  #   -------------:----------------------------------------
  #               p: $[[1, 2], [3, 4]]
  #                = [2×2 int64]
  #               q: ${p[1, 1]}
  #                = 4
  #               r: @{p}[:,1] + 1
  #                = [3 5] (int64)
  #               s: @{p}[:, 1].reshape(-1, 1) @ @{r}
  #                = [2×2 int64]
  #               t: np.linalg.eig(@{s})
  #                = EigResult(eigenvalue [...] 576, -0.89442719]]))
  #               w: ${t.eigenvalues[0]}  [...]  ${t.eigenvalues[1]}
  #                = 26.0
  #               x: $[[0,${t.eigenvalues [...] {t.eigenvalues[1]}]]
  #                = [0 26] (double)
  #   -------------:----------------------------------------
  # parameter list (param object) with 7 definitions


#%% Math example with DSCRIPT (pending) - v 1.005
    from pizza.dscript import dscript
    from pizza.private.mstruct import param

    D = dscript(name="math example")

    # The definitions are given with hybrid Matlab/NumPy notations
    D.DEFINITIONS.l = [1e-3, 2e-3, 3e-3]    # l is defined as a list
    D.DEFINITIONS.a = "$[1 2 3]"            # a is defined with  Matlab notations
    D.DEFINITIONS.b = "$[1:3]"              # b is defined with  Matlab notations
    D.DEFINITIONS.c = "$[0.1:0.1:0.9]"      # c is defined with  Matlab notations
    D.DEFINITIONS.scale = "@{l}*2*@{a}"     # l is rescaled
    D.DEFINITIONS.x0 = "$[[[-0.5, -0.5],[-0.5, -0.5]],[[ 0.5,  0.5],[ 0.5,  0.5]]]*${scale[0,0]}*${a[0,0]}"
    D.DEFINITIONS.y0 = "$[[[-0.5, -0.5],[0.5, 0.5]],[[ -0.5,  -0.5],[ 0.5,  0.5]]]*${scale[0,1]}*${a[0,1]}"
    D.DEFINITIONS.z0 = "$[[-0.5 0.5 ;-0.5 0.5],[ -0.5,  0.5;  -0.5,  0.5]]*${l[2]}*${a[0,2]}"
    D.DEFINITIONS.X0 = "@{x0}.flatten()"
    D.DEFINITIONS.Y0 = "@{y0}.flatten()"
    D.DEFINITIONS.Z0 = "@{z0}.flatten()"
    T = D.DEFINITIONS.eval()
    D.DEFINITIONS.x0
    T.x0
    print(repr(D.DEFINITIONS))

    # Output
    #  -------------:----------------------------------------
    #              l: [0.001, 0.002, 0.003]
    #              a: $[1 2 3]
    #              = [1 2 3] (int64)
    #             b: $[1:3]
    #              = [1 2 3] (int64)
    #             c: $[0.1:0.1:0.9]
    #              = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9] (double)
    #         scale: @{l}*2*@{a}
    #              = [0.002 0.008 0.018] (double)
    #            x0: $[[[-0.5, -0.5],[-0. [...] cale[0,0]}*${a[0,0]}
    #              = [[2×2 matrix] [2×2 matrix]] (2×2×2 double)
    #            y0: $[[[-0.5, -0.5],[0.5 [...] cale[0,1]}*${a[0,1]}
    #              = [[2×2 matrix] [2×2 matrix]] (2×2×2 double)
    #            z0: $[[-0.5 0.5 ;-0.5 0. [...] ]]*${l[2]}*${a[0,2]}
    #              = [[2×2 matrix] [2×2 matrix]] (2×2×2 double)
    #            X0: @{x0}.flatten()
    #              = [-0.001 -0.001 -0.001 -0.001 0.001 0.001 0.001 0.001] (double)
    #            Y0: @{y0}.flatten()
    #              = [-0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008] (double)
    #            Z0: @{z0}.flatten()
    #              = [-0.0045 0.0045 -0.0045 0.0045 -0.0045 0.0045 -0.0045 0.0045] (double)
    #  -------------:----------------------------------------
