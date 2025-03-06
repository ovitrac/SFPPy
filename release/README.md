## SFPPy <small>v1.3</small>  🍏⏩🍎  

**SFPPy** is a **Python-based framework** for **compliance testing of food contact materials** and **recycled plastic safety assessment**, supporting regulations from:  

- 🇺🇸 **US FDA regulations**  
- 🇪🇺 **European Union (EFSA, EU 10/2011, etc.)**  
- 🇨🇳 **Chinese GB standards**  
- 🌍 **Other international guidelines**  

SFPPy brings **well-established chemical migration models** to **pure Python**, offering **object-oriented** and **scalable** capabilities.  
🔍 Read the detailed internal documentation [here](https://ovitrac.github.io/SFPPy).

---

## 📦 Download & Installation  

- **Full version** (`< 1 MB` with documentation) - for Windows, Linux, macOS:  
  📥 **[Download SFPP](https://github.com/ovitrac/SFPPy/releases/)**  

#### 📂 Example Files  

Included in the **root folder**:  

- `example1.py`  
- `example2.py`  
- `example3.py`  
- `example4.py` 

---

## 🚀 Quick Start & Exploration  

📚 **Read the Online Documentation:**  

- [SFPPy Documentation](https://ovitrac.github.io/Pizza3/) 
- [Migration Modeling Guide](https://ovitrac.github.io/SFPPy/MigrationModeling/) 
- [FitNESS E-learning Platform](https://fitness.agroparistech.fr/) - ([Details](https://pubs.acs.org/doi/10.1021/acs.jchemed.4c00137))  

⚡ **Create complex migration scenarios with concise, efficient Python code:**  

<details>
  <summary>🔍 Click to expand example</summary>


```python
from patankar.layer import gPET, PP
from patankar.food import ambient, hotfilled, realfood, fat, liquid, stacked
from patankar.loadpubchem import migrant

# 🔬 Look up substances on PubChem
m = migrant("limonene")

# 🏗️ Create a multilayer material (ABA: PET-PP-PET)
A = gPET(l=(20, "um"), migrant=m, C0=0)
B = PP(l=(500, "um"), migrant=m, C0=200)  
ABA = A + B + A  # The leftmost layer is in contact with food

# 🍽️ Define contact & storage conditions
class contact1(stacked, ambient): name = "1: Setoff"; contacttime = (4, "months")
class contact2(hotfilled, realfood, liquid, fat): name = "2: Hot Filling"
class contact3(ambient, realfood, liquid, fat): name = "3: Storage"; contacttime = (6, "months")

# 🔄 Simulate a multiple-step migration process
medium1, medium2, medium3 = contact1(), contact2(), contact3()
medium1 >> ABA >> medium1 >> medium2 >> medium3  # Automatic chaining

# 📊 Merge & visualize migration kinetics
sol123 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation
sol123.plotCF()
```

</details>

---

## 💡 Why Use SFPPy?  

✅ **Regulation-Ready**: Supports compliance testing for **FDA, EFSA, and GB standards**  
✅ **Fast & Scalable**: Optimized finite-volume solver for **1D mass transfer modeling**  
✅ **Modular & Object-Oriented**: Easily define **multilayer materials & food contact conditions**  
✅ **PubChem Integration**: Fetch molecular properties directly from **PubChem**  
✅ **Automatic Chaining**: Simulate multi-step processes effortlessly  

🚀 **Try SFPPy now and streamline your food safety assessments!**  

---

### 🔗 Additional Resources  

📖 **Publication**:  
[AIChE Journal - SFPPy Migration Modeling](https://doi.org/10.1002/aic.14056)  

📌 **GitHub Pages**:  
[https://ovitrac.github.io/SFPPy/](https://ovitrac.github.io/SFPPy/)



***



<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>
