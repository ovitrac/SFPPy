# Running `SFPPy` without installation🍏⏩🍎

## Principles
`SFFPy` can be run without any specific installation (e.g., in [Spyder](https://www.spyder-ide.org/)):

- from the 📂`SFPPy/ `
- from anywhere if the 📂`SFPPy/ `is in Python `sys.path`

> The  📂`SFPPy/ ` is folder containing the folder 📁`patankar/` and 📄`example1.py`,📄`example2.py`…
>
>  📂`SFPPy/ ` is the folder is created when you unzip the last downloaded release ([download page](https://github.com/ovitrac/SFPPy/releases)).  Depending on the version, this folder is named `SFPPy_v1.30/` or `SFPPy_v1.30min/`.

## Update `sys.path` temporarily to include  📂`SFPPy/ `

In any Python 📄 script which is not located in  📂`SFPPy/ `, add at the beginning;


```python
# %% Temporarily sys.path modification
import sys

# Check if SFPPy path is already in sys.path
SFPPy_path = '/replace/it/by/your/path/SFPPy_v1.30' # <---- update it
if SFPPy_path in sys.path:
    print(f"SFPPy installation '{SFPPy_path}' is already in sys.path.")
else:
    print(f"Add SFPPy path '{SFPPy_path}' to sys.path.")
    sys.path.append(SFPPy_path)
```

    Add SFPPy path '/replace/it/by/your/path/SFPPy_v1.30'to sys.path.

## Let’s check if you can import substances (just an example)

You script can continue afterwards as usual:


```python
# %% Place your code below
from patankar.loadpubchem import migrant
list_of_migrants = ["toluene","limonene","BHT","Irganox 1076"]
m = [migrant(s) for s in list_of_migrants]
print(m)
```

> ⚠️ All imports of `SFPPy` modules must include `patankar`. 
>
>  ❌ Do not add 📁`patankar/` to `sys.path` since the same import rules are used internally. 📁`patankar/`is the core library of `SFPPy` but not the application itself, which include examples, its documentation…

## List of modules

|                               |                                                              |
| ----------------------------: | ------------------------------------------------------------ |
|          `patankar.migration` | 🏗️: Core solver implementing the Patankar finite-volume method for mass transfer modeling. |
|           `patankar.geometry` | 📐: Defines and processes 3D packaging geometries, including surface-area-to-volume calculations<br><br>**Use `from patankar.geometry import help_geometry` and `help_geometry()` to list all predefined geometries.** |
|               `patankar.food` | 🍎: Models food layers and their interactions with packaging materials. <br><br>💡**Use `from patankar.food import help_food` and `help_food()` to list all predefined food, simulants and contact conditions.** |
|              `patankar.layer` | 📜: Manages materials and multilayer packaging structures.<br><br>💡 **Use `from patankar.layer import help_layer` and `help_layer()` to list all predefined polymers and materials.** |
|           `patankar.property` | 📊: Computes essential physical and chemical properties (*e.g.,* diffusivities, partitioning coefficients). |
|        `patankar.loadpubchem` | 🔬: Interfaces with PubChem to retrieve molecular properties and cache them locally. |
| `patankar.private.somemodule` | 🔒Contains internal modules and utilities not intended for direct usage by end-users. For example: `patankar.private.pubchempy` |

> Each module contains many classes and subclasses. Import only needed classes (do not use `from.patankar.module import *` or `import patankar.module`). The methods are documented and attached to objects.
>
> Use helpers: `help_geometry()`, `help_food()`, `help_layer()`to access 

***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>