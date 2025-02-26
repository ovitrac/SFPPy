# SFPPy - Python Framework for Food Contact Compliance and Risk Assessment 🍏

## 🛠️ Overview

SFPPy is a Python-based framework for **compliance testing of food contact materials** and **recycled plastic safety assessment** under:

- 🇺🇸 **US FDA regulations**
- 🇪🇺 **European Union regulations** (EFSA, EU 10/2011, etc.)
- 🇨🇳 **Chinese GB standards**
- 🌍 **Other international guidelines**

This project **translates well-established chemical migration models** from MATLAB and other languages to **pure Python**, ensuring minimal dependencies.

## 📁 Main Modules (Located in `patankar/`)

- **`migration.py`** 🏗️ - Core solver using a Patankar finite-volume method for mass transfer modeling.
- **`geometry.py`** 📐 - Defines 3D packaging geometries and calculates volume/surface area.
- **`food.py`** 🍎 - Models food layers and their interactions with packaging.
- **`layer.py`** 📜 - Defines materials and layers for multilayer packaging.
- **`property.py`** 📊 - Computes physical and chemical properties (e.g., diffusion, partitioning).
- **`loadpubchem.py`** 🔬 - Retrieves molecular properties from PubChem (cached locally).

### Why Patankar?

<details>
  <summary>📜 Click to expand</summary>

> 💡 The `patankar` folder is named in honor of **Suhas V. Patankar**, who developed and popularized the **finite volume method**, which this project adapts for **mass transfer problems with an arbitrary number of Rankine discontinuities**.
>
> 🔧 The modules include a knowledge management system via extensible classes, allowing easy expansion to cover additional cases and implement new prediction methods.

</details>

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/ovitrac/SFPPy.git
cd SFPPy

# Install dependencies
pip install -r requirements.txt
```

## 💡 Usage Snippets

### Snippet 1: Simple Migration Simulation

<details>
  <summary>📜 Click to expand</summary>

```python
from patankar.food import ethanol  # food database
from patankar.layer import layer  # material database
from patankar.migration import senspatankar  # solver

# Define medium and layers
simulant = ethanol()
A = layer(layername="layer 1 (contact)", D=1e-15, l=50e-6, C0=0)  # SI units
B = layer(layername="layer 2", D=(1e-9, "cm**2/s"), l=(100, "um"))
multilayer = A + B  # layer A is contact (food is on the left)

# Run solver
solution = senspatankar(multilayer, simulant)
solution.plotCF()  # concentration kinetic in the simulant (F) for default times
solution.plotCx()  # concentration profile in the multilayer packaging
```

📝 **Notations**: $D$ is the diffusivity, $l$ is the thickness layer, and $C_0$ is the initial concentration.

</details>

### Snippet 2: Retrieving Molecular Properties

<details>
  <summary>🔍 Click to expand</summary>

```python
from patankar.loadpubchem import migrant  # connect to pubchem for missing substances

m = migrant(name="bisphenol A")
print(m.M, m.logP)  # Molecular weight & logP value
```

💡 The examples show how to inject `m` into layers (e.g., `multilayer` in snippet 1) to get customized simulations for specific substances and polymers.

</details>

### Snippet 3: Defining a Custom Packaging Shape

<details>
  <summary>📦 Click to expand</summary>

```python
from patankar.geometry import Packaging3D

pkg = Packaging3D('bottle', body_radius=(5, 'cm'), body_height=(0.2, 'm'),
                  neck_radius=(19, "mm"), neck_height=(40, "mm"))
vol, area = pkg.get_volume_and_area()
print("Volume (m³):", vol)
print("Surface Area (m²):", area)
```

💡 The examples show how to use either `pkg` or its properties to achieve mass transfer simulation for a specific geometry.

⚠️ **Note**: To efficiently simulate the migration of substances from packaging materials, SFPPy **unfolds complex 3D packaging geometries** into an equivalent **1D representation**. This transformation assumes that **substance desorption is predominantly governed by diffusion within the walls** of the packaging.

🔍 The `geometry.py` module provides tools to compute **surface-area-to-volume ratios**, extract wall thicknesses, and generate equivalent **1D models** for mass transfer simulations.

</details>

## 📖 Case Studies

The project includes three case studies: `example1.py`, `example2.py`, and `example3.py`, illustrating how real-world scenarios—featuring multiple materials, geometries, substances, food types, and usage conditions—can be numerically resolved.

### Example 1: **Mass Transfer from Monolayer Materials**

- 🥪 Simulates the migration of **Irganox 1076** and **Irgafos 168** from a **100 µm LDPE film** into a **fatty sandwich** over **10 days at 7°C**.
- 📈 Evaluates **migration kinetics** and their implications for food safety.

### Example 2: **Mass Transfer in Recycled PP Bottles**

- 🍼 Investigates **toluene migration** from a **300 µm thick recycled PP bottle** into a **fatty liquid food**.
- 🛡️ Assesses the **effect of a PET functional barrier** (FB) of varying thickness on reducing migration.

### Example 3: **Advanced Migration Simulation with Variants**

- 📦 Simulates migration in a **trilayer (ABA) multilayer system**, with **PET (A) and recycled PP (B)**.
- 🔥 Evaluates migration behavior across **storage with set-off, hot-filling, and long-term storage conditions**.
- ⚙️ Explores **variants** where the migrant and layer thickness are modified to assess performance.

⚠️ **Disclaimer**: These examples do not discuss sources of uncertainty. Please refer to our publications for details on the limitations of the presented approaches and assumptions.

## 📜 License

This project is licensed under the **MIT License**.

## 🤝 Contributors

**INRAE** - [Olivier Vitrac](mailto:olivier.vitrac@agroparistech.fr)  
*This project is part of the SFPPy initiative, aiming to bring the SafeFoodPackaging Portal version 3 (SFPP3) to the general public.*

$2025-02-12$
