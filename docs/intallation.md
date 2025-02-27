# **Installing SFPPy** 🍏⏩🍎

SFPPy supports **pip (default installation method)** and **Conda (for Anaconda/Miniconda users)**.  
Follow the instructions below to set up SFPPy on your system.

### **🔹 Need a Python Distribution?**
We recommend installing **[Anaconda](https://www.anaconda.com/download)** and following the **Installing with Conda** section.  
- **Windows users**: Python 3.11+ is also available in the **[Windows Store](https://apps.microsoft.com/detail/9nrwmjp3717k?ocid=webpdpshare)**.

---

## 📌 **Table of Contents**
[TOC]

---

## **📌 Installing with Pip (Default)**
For most users, installing SFPPy with `pip` is the simplest method.

> **What is Pip?**  
> Pip is a package manager for Python that installs libraries and dependencies directly from the Python Package Index (PyPI).

### **🔹 Step 1: Clone the SFPPy Repository**
```bash
git clone https://github.com/ovitrac/SFPPy.git
cd SFPPy
```

### **🔹 Step 2: Install Dependencies**
```bash
pip install -r requirements.txt
```

💡 **Windows users:** If `pip` is missing, try:
```bash
py -m pip install -r requirements.txt
```

### **🔹 Step 3: Verify Installation**
To ensure SFPPy is installed correctly, run:

#### **Linux/macOS**
```bash
python3 -c "import patankar; print('SFPPy is installed successfully')"
```

#### **Windows**
```bash
py -c "import patankar; print('SFPPy is installed successfully')"
```

---

## **📌 Installing with Conda (Anaconda/Miniconda)**
If you prefer **Conda**, follow these steps to create an isolated SFPPy environment.

> **Need to Install Anaconda or Miniconda?**  
> - [Anaconda Download](https://www.anaconda.com/download) – full distribution with many built-in tools.  
> - [Miniconda Installation Guide](https://docs.anaconda.com/miniconda/install/#quick-command-line-install) – lightweight alternative to Anaconda.

### **🔹 Step 1: Create the SFPPy Environment**
From the `SFPPy` directory, initialize Conda (usually with `conda init`), then run:
```bash
conda env create -f environment.yml
```
This creates an environment named **`sfppy`**.

### **🔹 Step 2: Activate the SFPPy Environment**
```bash
conda activate sfppy
```

### **🔹 Step 3: Verify Installation**
```bash
python -c "import patankar; print('SFPPy is installed successfully')"
```

---

## **📌 Manual Installation of Dependencies**
If you encounter issues, install the required dependencies manually:

#### **For Pip Users**
```bash
pip install numpy matplotlib scipy pandas openpyxl
```

#### **For Conda Users**
```bash
conda install numpy matplotlib scipy pandas openpyxl
```

---

## **📌 Troubleshooting Installation Issues**
If you experience issues, refer to the table below:

| **Issue**                                         | **Possible Cause**                | **Solution**                                                 |
| ------------------------------------------------- | --------------------------------- | ------------------------------------------------------------ |
| `python: command not found`                       | Python is not installed           | Install Python 3.x from [python.org](https://www.python.org/downloads/). |
| `ModuleNotFoundError: No module named 'patankar'` | `PYTHONPATH` is not set correctly | Run `export PYTHONPATH=$(pwd)` before executing the script.  |
| `Environment.yml not found`                       | Missing file in repository        | Ensure you are in the `SFPPy` directory before running `conda env create -f environment.yml`. |

For additional support, visit the **[Main Help Page](./index.html)**.

---

## **📌 Recommended Tools**
For an improved user experience, we recommend:

- **[Jupyter Notebook](https://docs.jupyter.org/en/latest/)** – ideal for interactive coding and visualization.
- **[Spyder](https://docs.spyder-ide.org/current/index.html)** – a Python IDE similar to MATLAB.

📌 **Windows users:** Follow this [Jupyter installation guide](https://www.geeksforgeeks.org/install-jupyter-notebook-in-windows/) to set up Jupyter.

---

## **📌 Summary**
| **Installation Method** | **Best For**                                   | **Installation Command**                          |
|------------------|-------------------------------------|------------------------------------|
| **Pip (Default)** | Most users with Python installed   | `pip install -r requirements.txt`  |
| **Conda (Anaconda/Miniconda)** | Users preferring Conda environments | `conda env create -f environment.yml`  |

---