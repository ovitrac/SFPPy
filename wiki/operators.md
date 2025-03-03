# **SFPPy Syntax: Operators `+` and `>>` in Migration Simulations** 🚀

## **Introduction**

SFPPy enables sophisticated mass transfer simulations in food packaging systems using a compact, Pythonic syntax. This page details the use of two key operators:

- **`+` (Addition Operator)** ➕: For **combining layers or results**
- **`>>` (Right Shift Operator)** ⏩: For **chaining migration simulations**

These operators simplify simulations while maintaining full control over the process.

> Note: the transfer of food contact conditions `food >> material` can be also represented as ``food @ material``.

---

## **1. The ➕ Operator: Layer Combination and Result Merging**

### **1.1. Assembling Multilayer Packaging** 🏗️

In SFPPy, packaging materials are modeled as **layers**, which can be combined using the `+` operator.

```python
from patankar.layer import gPET, PP

A = gPET(l=(20, "um"))  # 20 µm PET
B = PP(l=(500, "um"))   # 500 µm PP

ABA = A + B + A  # Creates a trilayer structure (ABA)
```

### **1.2. Merging Simulation Results** 📊

After running simulations for different food contact conditions, results can be **combined** using `+`.

```python
fullsolution = foodcontact1.lastsimulation + foodcontact2.lastsimulation
```

This aggregates migration profiles from multiple steps into a single dataset.

---

## **2. The ⏩ Operator: Chaining Migration Steps**

### **2.1. General Usage** 🔗

The `>>` operator simplifies **sequential simulations**, ensuring automatic propagation of **geometry, temperature, and mass transfer states**.

```python
mypackaging >> myfood >> mymaterial >> myfood
```
or equivalently

```python
mypackaging @ myfood >> mymaterial >> myfood
```

This executes:

1. **Geometry matching** between `mypackaging` and `myfood`
2. **Temperature propagation** to `mymaterial`
3. **Mass transfer simulation** into `myfood`

### **2.2. Handling Multiple Food Contact Conditions** 🍽️

For sequential food contact scenarios:

```python
mypackaging >> foodcontact2 >> foodcontact1 >> mymaterial >> foodcontact1 >> foodcontact2
```

This ensures:

- **Geometry adaptation** for both food contacts
- **Thermal equilibration**
- **Mass transfer simulation**
- **State propagation** across conditions

Each step stores results, accessible via:

```python
foodcontact1.lastsimulation
foodcontact2.lastsimulation
```

To explicitly store the final result:

```python
result = ... step n >> step n+1
```

### **2.3. Updating Properties Mid-Simulation** 🔄

SFPPy allows **dynamic property updates** while chaining simulations:

```python
foodcontact1 >> mymaterial.update(substance="new migrant", l=new_thickness) >> foodcontact1
```

---

## **3. Example: Advanced Simulation Workflow** 🧪

Consider a **trilayer (ABA) packaging** with **limonene** migration:

```python
from patankar.food import fat, ambient, hotfilled
from patankar.loadpubchem import loadpubchem

# Define the migrant
m = loadpubchem('limonene')

# Define packaging structure
A = gPET(l=(20, "um"))
B = PP(l=(500, "um"), migrant=m, CP0=200)
ABA = A + B + A

# Define food contact conditions
medium1 = ambient(contacttemperature=(20, "degC"), contacttime=(4, "months"))
medium2 = hotfilled(contacttemperature=(90, "degC"))
medium3 = fat(contacttemperature=(30, "degC"), contacttime=(6, "months"))

# Simulate migration
medium1 >> ABA >> medium1 >> medium2 >> medium3

# Combine results
sol123 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation
sol123.plotCF()
```

This example demonstrates how **minimal syntax** manages complex simulation steps efficiently. ⚡

---

## **Conclusion** 🎯

SFPPy’s `+` and `>>` operators provide a powerful, intuitive syntax for modeling **mass transfer simulations** in food packaging systems. By leveraging these operators, users can **reduce code complexity while preserving flexibility**. 🏆