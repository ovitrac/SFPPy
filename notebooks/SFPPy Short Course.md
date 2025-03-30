# 🍏⏩🍎 SFPPy Short Course

### Comprehensive Scientific Concepts and Practical Implementation


---

## 🎯 Objectives and Scope of SFPPy
SFPPy models the **migration** of chemical substances from polymeric materials into food or simulants. It calculates the concentration in food ($C_F$) and compares it to regulatory limits, particularly the **Specific Migration Limit (SML)**, under realistic yet conservative assumptions.

---

## 📚 Fundamental Concepts

### 1️⃣ Mass Transfer Governing Equations
he migration of substances (migrants) through polymer materials is a slow process driven by random molecular motions, often described as a "random walk." Microscopically, these movements occur due to the thermal agitation of polymer chains and atoms, creating transient pockets or voids that migrants can occupy and explore. The frequency of void formation and renewal depends significantly on the polymer state: faster in rubbery polymers and slower in glassy ones. The molecular size and rigidity of the migrants further influence their diffusion rates—larger and more rigid substances diffuse more slowly compared to smaller, flexible molecules.

Macroscopically, migrants distribute evenly in all directions, forming concentration gradients only when substances are allowed to exit the material. Given that packaging materials are typically much thinner relative to the dimensions of the food in contact, concentration gradients are effectively one-dimensional (1D), directed toward the food.

Accordingly, **SFPPy** solves **Fick's second law** for each polymeric layer $j$ and each substance $i$:
$$
\frac{\partial C_j}{\partial t} = D_j \frac{\partial^2 C_j}{\partial x^2}
$$

- $C_j$: migrant concentration in layer $j$ [mg/kg] or [mol/m³]
- $D_j$: diffusion coefficient in layer $j$ [m²/s]
- $t$: time [s]
- $x$: position [m]

### 2️⃣ Boundary Conditions (Robin or third kind)

At the polymer-food interface (polymer side):

$$
- D_1 \left.\frac{\partial C_1}{\partial x}\right|_{x=0^+} = h\left(\frac{k_1}{k_0} C_{\text{interface}} - C_F\right)
$$

- $h$: mass transfer coefficient [m/s] (boundary layer in food)
- $k_i$: Henry-like (sorption) coefficient for layer $i$
- $C_{\text{interface}}$: concentration at the polymer-food interface, polymer side
- $C_F$: concentration in food [mg/kg]

The sign depends on the interface position: positive on left side (food is on the left in SFPPy).

### 3️⃣ Sorption Isotherms and Partitioning
**SFPPy** uses linear sorption isotherms with Henry-like coefficients:
$$
p = k \cdot C
$$

- $p$: partial pressure (chemical potential analogy)
- $k$: Henry-like coefficient

Partition coefficient between layers A and B:

$$
K_{AB} = \frac{k_B}{k_A}
$$

---

## 🔔 SFPPy Automatic Property Predictions

**SFPPy** automatically selects the best predictive models from molecular descriptors ⌬:

| **Diffusivity in polymer ($D$)** 🏎️💨, 🚗💨🛻💨🚛💨                  | **Henry-like coefficients ($k$):** 🧲                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| 📌 `Dpiringer` – Legacy overestimation model (many polymers🖇, temperature effects🌡️) | 📌 `kFH` – Flory-Huggins model for polymers, adhesives, simulants 🥫 with entropic correction |
| 📌 `DFV` – Hole Free-Volume model, accurate for recycled plastics ♻️ and decontamination processes 💦 | 📌 Combined rules considering texture and interactions 🍯      |
| 📌 `Dwelle` – Refined for PET and PS around glass transition $T_g$ |                                                              |

---

## 🌐 International Regulatory Alignment

While focusing on 🇪🇺 **EU regulations (EFSA, EU 10/2011)**, SFPPy methods and principles are equally applicable to:

- 🇺🇸 **US FDA regulations**
- 🇨🇳 **Chinese GB standards**
- 🌍 **Other international guidelines**

---

## 🧪 Managing Chemical Information in SFPPy

**SFPPy** facilitates managing chemical information effectively, enabling users to add, change, and compare substances seamlessly.

### 1️⃣ Adding a Substance
Define substances easily by retrieving their chemical properties directly:

```python
from patankar.migrant import migrant

substance = migrant("irgafos 168")
```

### 2️⃣ Changing the Substance
You can update substances in existing layers effortlessly:

```python
# Define a polymer layer initially containing Irganox 1076
LDPE_layer = polymer.LDPE(l=(100, "µm"), substance=migrant("Irganox 1076"), C0=1000)

# Change to another substance, Irgafos 168
LDPE_layer_with_new_substance = LDPE_layer.update(substance=migrant("irgafos 168"))
```

### 3️⃣ Comparing Substances
SFPPy allows comparison between multiple substances in the same configuration:

```python
from patankar.migration import CFSimulationContainer as store # Store results
from patankar.migration import senspatankar as solver # Mass transfer solver

# Create storage for comparison
results_store = store(name="Comparison", SML=10)

# Add first substance (Irganox 1076)
sim1 = solver(LDPE_layer, FOODlayer.update(substance=migrant("Irganox 1076")))
results_store.add(sim1, "Irganox 1076", "r")  # red

# Add second substance (Irgafos 168)
sim2 = solver(LDPE_layer_with_new_substance, FOODlayer.update(substance=migrant("irgafos 168")))
results_store.add(sim2, "Irgafos 168", "b")  # blue

# Compare the migration kinetics
results_store.plotCF(plotSML=False)
```

> This structured approach simplifies the handling of chemical data, making simulations precise, flexible, and regulatory-compliant.

---

## 🧩 Practical Examples in SFPPy

### 📖 Monolayer Migration (`example1.py`)

Sequential Information and mass transfer can be induced with operator <kbd>>></kbd>:

```python
F = foodlayer(contacttime=(10,"days"), volume=(1,"L"), surfacearea=(6,"dm**2"))
P = layer(D=(1e-10,"cm**2/s"), C0=1000, k=0.1, l=(100,"µm"))
sim = F >> P >> F
sim.plotCF()
```

Chemical information can be introduced with the operator <kbd>%</kbd>:

```python
from patankar.loadpubchem import migrant
m = migrant("BHT")
(m % F >> P >> F).plotCF()
```

Food-Packaging geometry `G` can be inherited with the operator <kbd><<</kbd>

```python
from patankar.geometry import Packaging3D
G = Packaging3D('box_container', height=(8, "cm"), width=(10, "cm"), length=(19, "cm"))
(m % F << G >> P >> F)
```



### 📚 Advanced Multilayer Chained Scenario (`example2` and `example3`)

A complex material can be assembled with operator `+`

```python
# create a tri layer material (ABA) with an initial contamination located in B
from patankar.layer import PP, gPET
A = gPET(l=(0.5,"mm"),C0=0) # glassy PET (very barrier of not plasticized)
B = PP(l=(20,"µm"),C0=100)
material = A + B + A 
```

Simulation for different contact conditions can be chained with <kbd>>></kbd>:

```python
res = material >> contact1 >> contact2 >> contact3
res.plotCF() # concentration kinetics in F
res.plotCx() # concentration profiles in material
```
or equivalently:

```python
(res = material >> contact1 >> contact2 >> contact3).plotCF()
```



---

## 🔬 Sensitivity Analysis and Parameter Fitting (`example4)
```python
D = layerLink("D", indices=0, values=1e-10)
k = layerLink("k", indices=0, values=0.1)
P.Dlink = D
P.klink = k

sim = F.migration(P)
exp_data = sim.pseudoexperiment(npoints=30, std_relative=0.01)
sim.fit(exp_data)
```

Dynamic fitting and validation using pseudo-experimental data.

---

## 📗 Further Documentation

- 📚 [Migration Modeling Book Chapter](https://ovitrac.github.io/SFPPy/MigrationModeling/)
- 💻 [SFPPy GitHub Repository](https://github.com/ovitrac/SFPPy)
- 🚀 [SFPPyLite](https://ovitrac.github.io/SFPPyLite)

---

🍏⏩🍎 **End of SFPPy Short Course Outline**