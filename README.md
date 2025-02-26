# SFPPy - Python Framework for Food Contact Compliance and Risk Assessment 🍏⏩🍎

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

> 💡 The `patankar` folder is named in honor of **Suhas V. Patankar**, who developed and popularized the **[finite volume method](https://catatanstudi.wordpress.com/wp-content/uploads/2010/02/numerical-heat-transfer-and-fluid-flow.pdf)**, which this project adapts for **mass transfer problems with an arbitrary number of Rankine discontinuities**.
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

SFPPy is fully object-oriented and supports multiple syntax styles, ranging from a functional approach to a more abstract, operator-driven paradigm—all in a **Pythonic** manner. The snippets below demonstrate both approaches.

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

<small>💡 The examples show how to inject `m` into layers (e.g., `multilayer` in snippet 1) to get customized simulations for specific substances and polymers.</small>

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

<small>💡 The examples show how to use either `pkg` or its properties to achieve mass transfer simulation for a specific geometry.</small>

<small>⚠️ **Note**: To efficiently simulate the migration of substances from packaging materials, SFPPy **unfolds complex 3D packaging geometries** into an equivalent **1D representation**. This transformation assumes that **substance desorption is predominantly governed by diffusion within the walls** of the packaging.</small>

<small>🔍 The `geometry.py` module provides tools to compute **surface-area-to-volume ratios**, extract wall thicknesses, and generate equivalent **1D models** for mass transfer simulations.</small>

</details>

### Snippet 4: Using  **⏩**  as Mass Transfer Operator in Chained Simulations
<details>
 <summary>📦 Click to expand</summary>

📌 **SFPPy** leverages **multiple inheritance** to define food contact conditions by combining **storage conditions**, **food types**, and **physical properties**.  

📌 Additionally, **two operators** play a key role in SFPPy’s intuitive syntax:  

- **➕** for **combining layers** and **merging results**  
- **⏩** for naturally representing **mass transfer**  

With these operators, **mass transfer** can be abstracted into a simple, visual representation:  

1. **🍏⏩🍎**  
   _(Direct transfer from green to red, symbolizing migration.)_  

2. **🍏⏩🟠⏩🍎**  
   _(Includes an intermediate step, depicting progressive migration.)_  

3. **🍏⏩🟡⏩🟠⏩🍎**  
   _(More detailed, illustrating multiple contamination stages over time.)_  

4. **🍏⚡⏩🍎**  
   _(Emphasizes **active food transformation**, with accelerated mass transfer.)_  

🌟 **SFPPy** makes this abstraction possible with simple, expressive code.

```python
from patankar.layer import gPET, PP
from patankar.food import ambient, hotfilled, realfood, fat, liquid, stacked
from patankar.loadpubchem import migrant

# Define migrant and packaging layers (ABA: PET-PP-PET)
m = migrant("limonene")
A = gPET(l=(20, "um"), migrant=m, C0=0)
B = PP(l=(500, "um"), migrant=m, C0=200)  
ABA = A + B + A  # the most left layer is contact (food on the left)

# Define storage and processing conditions:
# 1:storage in stacks >> 2:hot-filled container >> 3:long-term storage of packaged food
class contact1(stacked, ambient): name = "1:setoff"; contacttime = (4, "months")
class contact2(hotfilled, realfood, liquid, fat): name = "2:hotfilling"
class contact3(ambient, realfood, liquid, fat): name = "3:storage"; contacttime = (6, "months")

# Instantiate and simulate with ⏩
medium1, medium2, medium3 = contact1(), contact2(), contact3()
medium1 >> ABA >> medium1 >> medium2 >> medium3  # Automatic chaining

# Merge all kinetics into a single one and plot the migration kinetics
sol123 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation
sol123.plotCF()

```
### **🧩 How It Works**

Each **contact class** inherits attributes from **multiple base classes**, allowing flexible combinations of:

1. **📌 Storage Conditions**:  
   - `ambient`: Defines standard storage at room temperature  
   - `hotfilled`: Represents high-temperature filling processes  
   - `stacked`: Models setoff migration when packaging layers are stacked  

2. **🥘 Food Types & Interactions**:  
   - `realfood`: Represents actual food matrices  
   - `liquid`: Specifies that the food is a liquid  
   - `fat`: Indicates a fatty food, influencing partitioning behavior  



<small>🔬 **By combining these components, SFPPy allows streamlined, physics-based simulations with minimal code.** 🚀</small>

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
- 🍏⏩🍎 Example 3 showcases the mass transfer operator ⏩.



> ⚠️ **Disclaimer**: These examples do not discuss sources of uncertainty. Please refer to our publications for details on the limitations of the presented approaches and assumptions.



## **🌟 Why SFPPy?**

✔ **SFPPy:** is **free** and **opensource**. 
✔ **SFPPy:** accepts any unit as `(value,"unit")` or `([value1,value2...],"unit")`. 
✔ **Operator-based chaining:** `>>` handles **automatic mass transfer and property propagation**
✔ **Minimal code for complex simulations:** `+` joins layers and merges results across storage conditions
✔ **Pythonic abstraction:** Works with **PubChem**, predefined **polymer materials**, and **3D packaging geometries**
✔ **Built-in visualization & export:** Supports **Excel (`.xlsx`), CSV, PDF, PNG** and **Matlab** (if its really needed)

🔬 **SFPPy powers scalable, real-world safe food packaging simulations.**

## 📜 License

This project is licensed under the **MIT License**.

## 🤝 Contributors

**INRAE** - [Olivier Vitrac](mailto:olivier.vitrac@agroparistech.fr)  
*This project is part of the SFPPy initiative, aiming to bring the SafeFoodPackaging Portal version 3 (SFPP3) to the general public.*

$2025-02-12$
