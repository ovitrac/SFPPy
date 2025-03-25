### SFPPy Utilities  🍏⏩🍎

---

# Utilities for Backup, Maintenance, and Synchronization + deployment via notebooks

---

## 🗃️ Notebook utilities

- **Import** `utils.nbutils`

- **Purpose**: Various utilities to deploy SFPPy notebooks

  

## 🧰 Available Scripts

> **All scripts are intended to be run from the `SFPPy/utils` directory.**




### 🔄 Backup Utilities

- **Run:** `./backupme.sh -y`  
  **Purpose:** Backup all important code and documentation.  
  **Output:**
  - `history/backupme.README.html`, `backupme.README.md` – Reports with links
  - `history/utils_backup_user@host_YYYY_MM_DD__HH-SS.zip` – Compressed backup

---

### 🧬 Code Analysis & Metadata Tools

- **Run:** `./generate_diagrams.sh`  
  **Purpose:** Generate the inheritance tree for all classes in `pizza/`.

- **Run:** `./generate_all.sh`  
  **Purpose:** Refresh all `__all__` metadata in Python codes.  
  **Output:** Updates `__all__` global variables across the project.

---

### 📄 Documentation Generators

- **Run:** `./generate_matlab_docs.py`  
  **Purpose:** Refresh Matlab/Octave documentation.  
  **Output:** `html/index_matlab.html` (single file, no dependencies)

- **Run:** `./pdocme.sh`  
  **Purpose:** Regenerate full HTML documentation.  
  **Output:** `html/index_matlab.html` (with additional linked files)

---

### 📦 Manifest and Release Tools

- **Run:** `./create_default_manifest.sh`  
  **Purpose:** Generate the main manifest (with hashes).  
  **Output:** `../Pizza3.manifest`

- **Run:** `./generate_simple_manifest.py`  
  **Purpose:** Generate a simpler version of the manifest.  
  **Output:** `../Pizza3.simple.manifest`

- **Run:** `./generate_release.sh`  
  **Purpose:** Build a release based on `Pizza3.simple.manifest`.  
  **Output:** `../release/*.zip`

---

## 🔁 Full Documentation Refresh Procedure

Run the following:

```bash
./refresh_alldocs.sh
```

> **Legacy Manual Steps (Obsolete):**
```bash
cd utils
rm -rf ../html/
./generate_matlab_docs.py
./generate_post_docs.py
./generate_diagrams.sh
./pdocme.sh
```

---

## 🚀 Release Creation Procedure

Run the following:

```bash
./refresh_allreleases.sh
```

> **Legacy Manual Steps (Obsolete):**
```bash
cd utils
./generate_simple_manifest.py
# edit and run
./generate_release.sh
```

---

## ⚙️ Setup and Packaging Refresh

To refresh the following packaging files:

- `setup.py`
- `requirements.txt`
- `MANIFEST.in` *(requires `generate_simple_manifest.py` first)*

Use these scripts from within `utils/`:

```bash
./generate_requirements.py
./generate_manifest_in.py
./generate_setup.py
```

---

## 📁 Project Directory Tree (Partial)

```plaintext
Pizza3/
│
├── utils/
│   ├── generate_requirements.py   # SETUP script
│   ├── generate_manifest_in.py    # SETUP script
│   └── generate_setup.py          # SETUP script
│
├── patanaker/
│   ├── __init__.py (if any)
│   ├── private/
│   │   ├── __init__.py (if any)
│   │   └── pint/
│   │       ├── __init__.py (if any)
│   │       └── ...
│   └── ... (other modules)
│
├── example1.py
├── example2.py
├── example3.py
...
├── README.md
├── LICENSE
├── SFPPy.simple.manifest
├── requirements.txt              # generated from utils/
├── MANIFEST.in                   # generated from utils/
└── setup.py                      # generated from utils/
```

---

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>
