#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Automatic Python module class and method documentation generation with Mermaid diagrams
#
# Inputs:
#  - Takes a list of modules from STDIN (one per line, in dot notation)
#  - First argument: output Markdown file
#  - Second argument: path to modules_details.json
#
# Outputs:
#  - A single Markdown file containing:
#    1) A tabulated list of methods for each module (method name, first paragraph of docstring, number of lines, __version__)
#    2) A Mermaid graph diagram for the class inheritance of that module
#    3) Links to class examples in HTML documentation
#    4) Date of generation

# INRAE\\Olivier Vitrac
# Email: olivier.vitrac@agroparistech.fr
# Last Revised:** 2025-01-09

import sys
import os
import importlib
import inspect
from datetime import datetime
import json

# Configuration Preambule
ConfigurationPreambule = """
## **Running SFPPy**

> **SFPPy** enables the **compliance assessment of food contact materials** and facilitates **risk assessments** through **migration modeling**.  
> The term *"[Migration Modeling]("MigrationModeling/")"* refers to the **simulation of mass transfer** (primarily diffusive) from packaging materials into food.

#### **What this guide covers:**
1. **Running SFPPy/example1.py** (or `example2.py`, `example3.py`) immediately after cloning the repository.
2. **Using the object-oriented `patankar` library** to solve custom migration modeling problems.

SFPPy runs on **Windows, Linux, macOS**, and any system supporting **Python 3.x**.  
For the best experience, we recommend using [Jupyter](https://docs.jupyter.org/en/latest/) notebooks, which allow you to execute and document your workflow in a single document—especially useful for compliance testing.

---

### **📌 Prerequisites**
Before running SFPPy, make sure you have:

✅ **A basic understanding of Migration Modeling** → Read this [document]("MigrationModeling/").  
✅ **Python 3.x (≥3.8) installed**.  
✅ **Dependencies installed**:
   - If your Python distribution does **not include scientific libraries** (NumPy, SciPy, Matplotlib, Pandas), install them using:
     ```bash
     pip install -r requirements.txt
     ```
   - Alternatively, follow the **[installation guide](installation.html)** for Conda-based setup.  
✅ **An interactive environment** (recommended):
   - Install [Jupyter](https://docs.jupyter.org/en/latest/) or [Spyder](https://docs.spyder-ide.org/current/index.html) for a better experience.

---

### **📌 Running a Script with SFPPy**

This section explains how to **run `SFPPy/yourscript.py`**, assuming your scripts are stored inside the `SFPPy` directory alongside the example scripts (`example1.py`, `example2.py`, `example3.py`).

Using this approach **avoids system-wide modifications** and ensures that all SFPPy modules are found.

---

#### **🔹 Step 1: Navigate to the SFPPy directory**
```bash
cd SFPPy
```
This ensures that you are inside the **SFPPy project folder** before running a script.

---

#### **🔹 Step 2: Run Your Script**

##### **Option A: Standard Python Execution**
Use the appropriate command based on your operating system:

| **OS**                   | **Shell**      | **Command** ^(*)^                                   |
| ------------------------ | -------------- | --------------------------------------------------- |
| **Linux/macOS**          | `bash`         | `export PYTHONPATH=$(pwd) && python3 yourscript.py` |
| **Windows (cmd.exe)**    | Command Prompt | `set PYTHONPATH=%CD% && py -3 yourscript.py`        |
| **Windows (PowerShell)** | PowerShell     | `$env:PYTHONPATH="$PWD"; & py -3 yourscript.py`     |

<small>^*^ Replace `yourscript.py` with your script (e.g., `example1.py`).</small>

---

#### **🔹 Alternative: Running SFPPy Inside a Conda Environment**

**Option B: Standard Python Execution**. If you installed SFPPy via Conda, follow these steps:

##### **Step 2a: Create the SFPPy Conda Environment**
Run the following command **from the `SFPPy` directory**:
```bash
conda env create -f environment.yml
```
This creates a Conda environment named **`sfppy`** with all required dependencies.

##### **Step 2b: Activate the Environment**
```bash
conda activate sfppy
```

Once the environment is activated, run your script:
```bash
python yourscript.py # or `conda activate sfppy && python yourscript.py`
```
"""

def class_name(cls):
    """Return a string representing the class"""
    return cls.__name__

def classes_tree(module, base_module=None):
    """Extract classes and their inheritances starting from a base module."""
    if base_module is None:
        base_module = module.__name__
    module_classes = set()
    inheritances = []
    def inspect_class(cls):
        if class_name(cls) not in module_classes:
            if cls.__module__.startswith(base_module):
                module_classes.add(class_name(cls))
                for base in cls.__bases__:
                    inheritances.append((class_name(base), class_name(cls)))
                    inspect_class(base)
    for cls in [e for e in module.__dict__.values() if inspect.isclass(e)]:
        inspect_class(cls)
    return module_classes, inheritances

def classes_tree_to_mermaid(module_classes, inheritances):
    lines = ["graph TD;"]
    for c in sorted(module_classes):
        lines.append(f"{c}")
    for (a, b) in sorted(inheritances):
        lines.append(f"{a} --> {b}")
    return "\n".join(lines)

def get_first_paragraph(docstring):
    if not docstring:
        return ""
    lines = docstring.strip().split('\n')
    # First paragraph is up to a blank line or the end of lines
    paragraph = []
    for line in lines:
        if line.strip() == "":
            break
        paragraph.append(line.strip())
    return " ".join(paragraph)

def get_method_info(obj):
    """Return a dict with method name, docstring first paragraph, number of lines, and __version__."""
    info = {}
    if not (inspect.isfunction(obj) or inspect.ismethod(obj)):
        return info

    name = obj.__name__
    doc = inspect.getdoc(obj)
    first_par = get_first_paragraph(doc)
    try:
        lines, _ = inspect.getsourcelines(obj)
        nlines = len(lines)
    except (OSError, TypeError):
        nlines = 0

    # Attempt to find __version__ attribute:
    # Methods typically won't have their own __version__, but we can check their __globals__ if function
    __version__ = None
    if inspect.isfunction(obj):
        __version__ = obj.__globals__.get('__version__', None)
    # If not found there, we can try attributes
    if hasattr(obj, '__version__'):
        __version__ = getattr(obj, '__version__', __version__)

    info['name'] = name
    info['doc'] = first_par
    info['nlines'] = nlines
    info['__version__'] = __version__ if __version__ is not None else ""
    return info

def get_module_methods(module):
    """Get methods and functions defined at the top-level and within classes of the module."""
    # Top-level functions
    top_methods = []
    for name, obj in inspect.getmembers(module, inspect.isfunction):
        if obj.__module__ == module.__name__:
            top_methods.append(get_method_info(obj))

    # Class methods
    class_methods = []
    for name, cls in inspect.getmembers(module, inspect.isclass):
        if cls.__module__ == module.__name__:
            for mname, mobj in inspect.getmembers(cls, inspect.isfunction):
                # To filter out inherited methods from outside modules if desired:
                if mobj.__module__ == module.__name__:
                    method_info = get_method_info(mobj)
                    method_info['class'] = cls.__name__
                    class_methods.append(method_info)

    return top_methods, class_methods

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_mermaid.py output_markdown.md modules_details.json < modules_list.txt", file=sys.stderr)
        sys.exit(1)

    output_markdown = sys.argv[1]
    modules_details_json = sys.argv[2]

    if not os.path.isfile(modules_details_json):
        print(f"Error: JSON file '{modules_details_json}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # === NEW SECTION: Load Module Details from JSON ===
    # This section loads the modules_details.json to access module details
    with open(modules_details_json, 'r', encoding='utf-8') as f:
        modules_details = json.load(f)
    # === END OF NEW SECTION ===

    # Read module names from STDIN
    modules = [line.strip() for line in sys.stdin if line.strip()]

    with open(output_markdown, "w", encoding="utf-8") as f:
        # Write the main header and generation date
        f.write(f"# SFPPy Modules Documentation\n\n")
        f.write(f"Generated on: **{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}**\n\n")
        f.write('<hr style="border: none; height: 1px; background-color: #e0e0e0;" />\n\n')
        # Write the configuration preambule
        f.write(ConfigurationPreambule + "\n\n\n\n")
        f.write('<hr style="border: none; height: 3px; background-color: #4CAF50;" />\n\n\n\n')
        #f.write('<hr style="border: none; height: 1px; background-color: #e0e0e0;" />\n\n')

        # Add a styled, responsive Table of Contents with three columns
        f.write('<a id="table_of_contents" name="table_of_contents"></a>\n')
        f.write("## Main Classes\n\n")
        f.write('<div id="table_of_contents" style="display: flex; flex-wrap: wrap; gap: 20px; justify-content: space-between; overflow-x: auto; padding: 10px;">\n')
        for idx, mod_name in enumerate(sorted(modules)):
            f.write('<div style="flex: 1 1 calc(33.33% - 20px); min-width: 200px;">\n')
            f.write(f'<a href="#{mod_name.replace(".", "_")}" style="text-decoration: none; font-weight: bold;">\n')
            f.write(f'{idx + 1}. {mod_name}\n')
            f.write('</a>\n')
            f.write('</div>\n')
        f.write('</div>\n\n')

        # Iterate over each module
        for i, mod_name in enumerate(sorted(modules)):
            # Attempt to import the module
            try:
                mod = importlib.import_module(mod_name)
            except Exception as e:
                f.write(f"## Module `{mod_name}`\n\n")
                f.write(f"**Error importing module**: {e}\n\n")
                continue

            # anchor
            f.write(f'<a id="{mod_name.replace(".", "_")}" name="{mod_name.replace(".", "_")}"></a>\n')
            # Generate anchors for TOC, previous module, and next module
            prev_mod = sorted(modules)[i - 1] if i > 0 else None
            next_mod = sorted(modules)[i + 1] if i < len(modules) - 1 else None
            # Navigation links container with alignment styling
            f.write('<div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px; font-size: 0.8em;">')
            # Add link to the previous module if it exists (left-aligned)
            if prev_mod:
                f.write(f'<a href="#{prev_mod.replace(".", "_")}" title="Go to Previous Module: {prev_mod}" style="text-decoration: none;">⬅️ Previous</a>\n')
            else:
                f.write('<span></span>\n')  # Empty span to maintain layout consistency
            # Add link to TOC (centered)
            f.write('<a href="#table_of_contents" title="Back to Table of Contents" style="text-decoration: none;">⬆️ TOC</a>\n')
            # Add link to the next module if it exists (right-aligned)
            if next_mod:
                f.write(f'<a href="#{next_mod.replace(".", "_")}" title="Go to Next Module: {next_mod}" style="text-decoration: none;">➡️ Next</a>\n')
            else:
                f.write('<span></span>\n')  # Empty span to maintain layout consistency
            # Close the navigation container
            f.write('</div>\n\n')

            # Start the module section
            f.write(f'## Module `{mod_name}`\n\n')

            # Get class inheritance tree
            module_classes, inheritances = classes_tree(mod, base_module='patankar')
            # Generate Mermaid diagram if classes are found
            if module_classes:
                f.write("### Class Inheritance Diagram\n")
                f.write("```mermaid\n")
                f.write(classes_tree_to_mermaid(module_classes, inheritances))
                f.write("\n```\n\n")
            else:
                f.write("*No classes found in this module.*\n\n")

            # === NEW SECTION: Add Link to Class Examples ===
            # This section retrieves module details from the JSON and adds a link to class_examples.html
            module_detail = modules_details.get(mod_name)
            if module_detail:
                local_path = module_detail.get("local_path", "")
                number_of_examples = module_detail.get("number_of_examples", 0)
                url_anchor = module_detail.get("url_anchor", "")
                # Construct relative link (e.g., "class_examples.html#patankar_dscript")
                link = f"class_examples.html{url_anchor}"
                # Add the link in bold with module name and number of examples
                f.write(f"**[Class Examples for `{local_path}` ({number_of_examples})]({link})**\n\n")
            else:
                # If module details are missing, indicate that class examples are not available
                f.write(f"**Class Examples:** Not available.\n\n")
            # === END OF NEW SECTION ===

            # Get methods of the module
            top_methods, class_methods = get_module_methods(mod)

            if top_methods or class_methods:
                f.write("### Methods Table\n\n")
                # Write the table header
                f.write("| Class | Method | Docstring First Paragraph | # Lines | __version__ |\n")
                f.write("|-------|---------|---------------------------|---------|-------------|\n")
                # Write top-level methods
                for m in top_methods:
                    f.write(f"| (module-level) | `{m['name']}` | {m['doc']} | {m['nlines']} | {m['__version__']} |\n")
                # Write class methods
                for m in class_methods:
                    f.write(f"| `{m.get('class','')}` | `{m['name']}` | {m['doc']} | {m['nlines']} | {m['__version__']} |\n")
                f.write("\n")
            else:
                f.write("*No methods found in this module.*\n\n")

if __name__ == "__main__":
    main()
