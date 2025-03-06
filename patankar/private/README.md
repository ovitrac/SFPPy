# Private Modules for SFPPy 🍏⏩🍎

This directory contains internal modules that are used internally by the SFPPy library. These modules provide additional utilities and extended functionalities but are not intended for direct interaction by the user.

## 📁 Files and Subdirectories
- `mstruct.py` 📦: Tools for handling structured data.
- `pint/` 📏: Standard **Pint** library used for SI unit conversions in SFPPy.
- `pubchempy.py` 🔬: Interface for retrieving chemical data from **PubChem**.
- `chemspipy/` ⚠️: Previously used for **ChemSpider** integration but now deprecated (requires API tokens).
- `toxtree/` ☠️: Local installation folder of the private copy of Toxtree.

## 🔹 Notes
- Modules in this folder are accessed indirectly through public-facing modules in `patankar/`.
- **Do not modify** files here unless you are extending or debugging SFPPy internals.
- **Avoid importing modules directly from `private/`**; instead, use the official `patankar.private.module` structure.



***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>

