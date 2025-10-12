### **🚀 SFPPy - Python Framework for Food Contact Compliance & Risk Assessment 🍏⏩🍎**



# **SFPPy** <small>v1.50</small> Release 🍏⏩🍎  

---

> [!IMPORTANT]
>
> This page provides the **latest stable version** of `SFPPy` for **local installation** on your machine.
>
> If you're looking for an **online experience**, there are two alternatives:
>
> - 🟢 **Google Colab**: No packaging needed —  follow the `git clone` instructions shown on the [main project page](https://github.com/ovitrac/SFPPy) or watch the [video walkthrough](https://ovitrac.github.io/SFPPy/wikipages/colab/colab_demo.html).  
>
> - 🌐 **SFPPyLite**: Use SFPPy **directly in your browser** with no setup at all. Powered by WebAssembly and JupyterLite, it gives you access to most of SFPPy’s core capabilities in a single click.  
>   <a href="https://ovitrac.github.io/SFPPylite/lab/index.html?path=demo.ipynb" target="_blank">
>     <img src="https://img.shields.io/badge/launch-SFPPylite%20in%20your%20browser-blueviolet?logo=jupyter&style=for-the-badge" alt="🧪 Try it online!">
>   </a>
>
> 🔄 Both `SFPPy` and `SFPPyLite` are aligned and current in version **1.44**.

---

## **About SFPPy**  

`SFPPy` is a Python-based framework for 𓌉◯𓇋 compliance testing of food contact materials and 🎢⌬♻️ recycled plastic safety assessment. It supports regulations from:

- 🇺🇸 **US FDA regulations**  
- 🇪🇺 **European Union (EFSA, EU 10/2011, etc.)**  
- 🇨🇳 **Chinese GB standards**  
- 🌍 **Other international guidelines**  

## Release Notes

📌 *`SFPPy` is AI-ready. **Yes**, GPT robots can develop chemical migration scenarios using the Pythonic language, operator and class models using the `SFPPy` frameworks.* The limitations are discussed [here](https://ovitrac.github.io/SFPPy/wikipages/#discussion/GPT_assistance.html). `Example5.py` demonstrates how a scientific agent can be built with `ollama`, **RAG**, and a proper knowledge base (provided as `docs/KP.zip`).

📌 *`SFPPy` is integrated in notebooks that keep records of the entire evaluation process.* Starting from v1.37, a [compliance template](https://ovitrac.github.io/SFPPy/wikipages/#notebooks/comply.html) notebook for compliance is shipped with `SFPPy`.

📌 *`SFPPy` is equipped with a graphical interface.* `SFPPy-GUI` supports many of the capabilities of `SFPPy` and can easily be extended with forms or customized widgets. The manual is available [here](https://ovitrac.github.io/SFPPy/wikipages/#gui/gui_manual.html). 

📌 *`SFPPy` can be used in hybrid environments (online/offline), scripts/notebooks, Windows/Linux.* All examples are provided as `exampleX.py` and `notebooks/exampleX.ipynb`, with X = 1–4. 

📌 Starting from version **1.3**, `SFPPy` starts all simulations from the sole knowledge of the chemical name of the substance, inferring the diffusivities (🏎️💨, 🚗💨🛻💨🚛💨) and partitioning (🧲) from the best models applicable. The collection of models can be augmented as shown in `patankar.properties`.

📌 *`SFPPy` uses PubChem as the primary source of chemical information, including chemical structure. An internet connection enables `SFPPy` to retrieve information automatically.* The information is cached locally to avoid overloading PubChem servers.

📌 *Starting from version **1.37**, `SFPPy` is shipped with all substances of Annex I of the Regulation (EU) 10/2011 and provides toxicological alerts (if found) for all substances used as inputs.* It uses ToxTree for that.

📌 *`SFPPy` is agnostic to the system considered and runs indifferently on 🪟 Windows, 🐧 Linux, 🍎 macOS.* SFPPy depends exclusively on standard libraries: SciPy, NumPy, Matplotlib, Pandas, PIL.

📌 *`SFPPy` is very fast and designed to test chained complex scenarios involving one or several substances.*

📌 *As of version **1.40**, `SFPPy` and `SFPPyLite` offer similar functionality. Users are encouraged to work through **Jupyter notebooks** rather than standalone Python scripts, as notebooks offer a more interactive and reproducible workflow. Global behavior—such as discretization resolution, plotting preferences, and more—can be customized via user overrides using:* 
 `from patankar.useroverride import useroverride`

📌 *As of version **1.41**, `SFPPy` and `SFPPyLite` are shipped with US FCN  (🇺🇸), EU Annex I (🇪🇺), and GB Appendix A  (🇨🇳) databases, all managed independently with their database manager and updater. All individual substances and their toxicological assessment are shipped with SPFPPy.  

📌 Version **1.44** brings several fixes and additions identified/requested by users.

📌  *As of version **1.50**, `SFPPy`  brings **design-for-compliance** to packaging with tiers **M0–M3** , *explicit performance indicators*, *open data*, *AI-assisted reasoning*. It automates **risk assessment** for the 100s of chromatogram peaks in **recycled-material extracts**. **Read our 📄 technical paper [here](https://ovitrac.github.io/SFPPy/wikipages/SafeByDesign/sfppydesign.html)**.

---

### 🗜️ Download SFPPy v1.50 – Stable Release

| Feature / Content                             | ✅ **Full Version**<br><kbd>129 MB</kbd>                      | ⚪ **Light Version**<br><kbd>66 MB</kbd>                      |
| --------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 🧪 Includes `Toxtree` integration              | ✅ Yes                                                        | ❌ No                                                         |
| 📚 Full documentation                          | ✅ Yes                                                        | ✅ Yes                                                        |
| 🧬 Advanced risk assessment (`migrantToxtree`) | ✅ Supported                                                  | ❌ Not supported                                              |
| 🔬 Routine migration modeling (`migrant`)      | ✅ Supported                                                  | ✅ Supported                                                  |
| 📥 Download ZIP                                | [Download SFPPy](https://github.com/ovitrac/SFPPy/releases/download/v1.5000-r01/SFPPy_v1.50.zip) | [Download SFPPymin](https://github.com/ovitrac/SFPPy/releases/download/v1.4400-r01/SFPPy_v1.50min.zip) |
| 📦 View all versions                           | [Releases](https://github.com/ovitrac/SFPPy/releases)        | [Releases](https://github.com/ovitrac/SFPPy/releases)        |

---

💡 **Note**: The **light version** does **not** include a private installation of `Toxtree`, used by `patankar.loadpubchem.migrantToxtree` for advanced risk assessment. However, `patankar.loadpubchem.migrant` remains **fully functional** for routine migration simulations.

---

> [!WARNING]
>
> ❗ *`Toxtree` requires a Java runtime environment ([`JRE`](https://www.java.com/en/download/manual.jsp)) if Java ☕ is not already installed.* If needed, `Toxtree` can be manually installed in `patankar/private/toxtree/` as described in `patankar/private/toxtree/README.md`.



---

## 📂 **Example Files**  

`SFPPy` includes **four example scripts** (and more notebooks) to help you get started:

| **Example file (located in the root folder 🏠)**              | **Description** (see 🔍 [Wiki Pages](https://ovitrac.github.io/SFPPy/wikipages/)) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [`example1.py`](https://ovitrac.github.io/SFPPy/wikipages/#examples/example1.html) | *Monolayer materials*                                        |
| [`example2.py`](https://ovitrac.github.io/SFPPy/wikipages/#examples/example2.html) | *Recycled bottles with functional barriers*                  |
| [`example3.py`](https://ovitrac.github.io/SFPPy/wikipages/#examples/example3.html) | *Chained simulations with variations*                        |
| [`example4.py`](https://ovitrac.github.io/SFPPy/wikipages/#examples/example4.html) | *Fitting experimental curves to extract diffusion and partition coefficients* |
| [`example5.py`](https://ovitrac.github.io/SFPPy/wikipages/#examples/example4.html) | *Retrieval-Augmented Generation (RAG) over Legal Knowledge Base* |

> [!TIP]
>
> ✨ **All scripts in the root folder 🏠 can be executed directly without modifying the Python path.**  
> If running from another location, use the CLI:  
>
> ```bash
> sfppy example1.py  # Replace example1.py with your script
> ```

---

> [!IMPORTANT]
>
> **Notebooks**: Keep code, assumptions, results, and interpretations together — ideal for regulatory workflows. If you run `notebooks/example1.ipynb`, you can save automatically all results and code with:
>
> ```python
> export_notebook(filename="example1",outputfolder="reports");
> ```



## 🛠️ **Running SFPPy from Anywhere 🌐**  

To test releases effortlessly:

```bash
# Navigate to the extracted SFPPy folder
cd /path/to/SFPPy  # 🏠 Unzipped folder

# Install SFPPy in editable mode (dependencies included)
pip install -e .

# Run example1 from anywhere
sfppy example1.py  # Replace example1.py with your script
```

For detailed installation instructions, refer to the [installation guide](https://ovitrac.github.io/SFPPy/wikipages/#installation/installation.html) 🧑‍🔧.

> [!NOTE]
>
> All notebooks located in 📂`notebooks/` are equipped with a bootstrap mechanism enabling them to start without any installation.



---

## 🚀 **Quick Start & Exploration**  

📚 **Read the Online Documentation:**  

- 📖 [SFPPy Documentation](https://ovitrac.github.io/SFPPy/)  
- 📚 [SFPPy Wiki Pages](https://ovitrac.github.io/SFPPy/wikipages/)  
- 📑 [Migration Modeling Guide](https://ovitrac.github.io/SFPPy/MigrationModeling/)  
- 🎓 [FitNESS E-learning Platform](https://fitness.agroparistech.fr/) ([Details](https://pubs.acs.org/doi/10.1021/acs.jchemed.4c00137))  

🚀 *Try SFPPy now and streamline your food safety assessments!*



---

## 🔗 **Additional Resources**

📖 **Publication on Chained Simulation**  
🧬 *AIChE Journal*: [**Migration Modeling with SFPPy**](https://doi.org/10.1002/aic.14056)  
A peer-reviewed reference for the core modeling engine and chaining logic implemented in SFPPy.

---

📌 **GitHub Pages**

[![SFPPy main](https://img.shields.io/badge/SFPPy-main-4CAF50?logo=github&style=for-the-badge)](https://ovitrac.github.io/SFPPy/)  
The official development hub for `SFPPy`, including installation, examples, and documentation.

[![SFPPyLite in your browser](https://img.shields.io/badge/SFPPy_Lite-in_your_browser-FF4D4D?logo=github&style=for-the-badge)](https://ovitrac.github.io/SFPPylite/)  
`SFPPyLite` is a sister project that runs entirely **in-browser**, without installation or servers. Most use cases—setup, simulation, visualization, and export—can be handled interactively in JupyterLite.  
See [wiki pages](https://ovitrac.github.io/SFPPy/wikipages/#SFPPylite/access.html) for more.

---

> [!TIP]
>
> 🧠 **Why is SFPPy AI-ready?**  
> Read more on [Generative Simulation](https://github.com/ovitrac/generativeSimulation).  
> GPTs and modern AI tools understand and generate Python code natively. The `SFPPy` API is designed for **semantic modeling**—using intuitive structures like operator chaining (`>>`), integrated unit handling, and chemical object orientation.

---



<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  <b>Contact</b> <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> <b>for questions |</b>
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> <b>|</b>
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>