# Private Modules for SFPPy 🍏⏩🍎

This directory contains internal modules that are used internally by the SFPPy library. These modules provide additional utilities and extended functionalities but are not intended for direct interaction by the user.

## 📁 Files and Subdirectories
- `mstruct.py` 📦: Tools for handling structured data.
- `pint/` 📏: Standard **Pint** library used for SI unit conversions in SFPPy.
- `pubchempy.py` 🔬: Interface for retrieving chemical data from **PubChem**.
- `chemspipy/` ⚠️: Previously used for **ChemSpider** integration but now deprecated (requires API tokens).
- `__pycache__/` ⚙️: Compiled Python bytecode files.

## 🔹 Notes
- Modules in this folder are accessed indirectly through public-facing modules in `patankar/`.
- **Do not modify** files here unless you are extending or debugging SFPPy internals.
- **Avoid importing modules directly from `private/`**; instead, use the official `patankar.private.module` structure.


