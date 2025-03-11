# SFPPy Core Library  🍏⏩🍎

This directory contains the core Python modules for the SFPPy framework.
These modules provide foundational functionalities for building, running, and analyzing migration simulations.

## 📁 Files and Subdirectories (Overview, Not Exhaustive)
- `__init__.py`: Initializes the SFPPy Python package.
- `migration.py` 🏗️: Core solver implementing the Patankar finite-volume method for mass transfer modeling.
- `geometry.py` 📐: Defines and processes 3D packaging geometries, including surface-area-to-volume calculations.
- `food.py` 🍎: Models food layers and their interactions with packaging materials.
- `layer.py` 📜: Manages materials and multilayer packaging structures.
- `property.py` 📊: Computes essential physical and chemical properties (e.g., diffusivities, partitioning coefficients).
- `loadpubchem.py` 🔬: Interfaces with PubChem to retrieve molecular properties and cache them locally.
- `private/`: Contains internal modules and utilities not intended for direct usage by end-users.

## 🔹 Notes
- **Do not include this directory in your `PYTHONPATH`**. Instead, always import modules as `patankar.module`.
- The `private/` subdirectory holds internal components and low-level implementations that support the public API but should not be accessed directly.

🚀 **SFPPy is designed to simplify compliance testing and migration modeling—keeping it modular, scalable, and efficient!**