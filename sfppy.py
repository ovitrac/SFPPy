"""
===============================================================================
SFPPy Command Line Interface (CLI)
===============================================================================

Provides a command-line tool (`sfppy`) to run SFPPy scripts from anywhere
without manually modifying `sys.path`. It automatically adds the SFPPy root
directory (which contains `patankar/` and example scripts) to `sys.path`.

**Features:**
- Detects and adds `SFPPy/` (the parent of `patankar/`) to `sys.path`.
- Runs Python scripts in the correct directory to handle relative imports.
- Enables running SFPPy scripts globally without needing to be in the project folder.

**Usage:**
```sh
sfppy example1.py
```
or with arguments:
```sh
sfppy example2.py --some-option value
```

**Installation:**
1. Ensure `cli.py` is placed inside `SFPPy/` (not inside `patankar/`).
2. Modify `setup.py` to register `sfppy` as a console script:
```python
entry_points={
    "console_scripts": [
        "sfppy=cli:main",
    ],
},
```
3. Install SFPPy:
```sh
pip install -e .
```
4. Now you can run SFPPy scripts from anywhere using `sfppy <script.py>`.

-------------------------------------------------------------------------------
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: Olivier Vitrac
@license: MIT
@date: 2025-03-07
===============================================================================
"""

import sys
import os

def main():
    """
    Entry point for the SFPPy CLI.

    This function:
    1. Checks if a script path is provided as an argument.
    2. Identifies `SFPPy/` (the parent of `patankar/`) as the project root.
    3. Adds `SFPPy/` to `sys.path` to allow imports like `from patankar.module import X`.
    4. Changes the working directory to the script's location for proper relative paths.
    5. Executes the script in a controlled namespace.

    Usage:
        sfppy <script.py> [args]

    Raises:
        SystemExit: If no script is provided.
    """
    if len(sys.argv) < 2:
        print("Usage: sfppy <script.py> [args]")
        sys.exit(1)

    script_path = os.path.abspath(sys.argv[1])  # Absolute path of the script

    if not os.path.isfile(script_path):
        print(f"Error: Script '{script_path}' not found.")
        sys.exit(1)

    # Locate SFPPy root directory (parent of patankar/)
    sfppy_root = os.path.dirname(os.path.abspath(__file__))  # This is SFPPy/

    # Ensure SFPPy contains patankar/
    if not os.path.isdir(os.path.join(sfppy_root, "patankar")):
        print(f"Error: The detected SFPPy root '{sfppy_root}' does not contain 'patankar/'.")
        sys.exit(1)

    # Ensure SFPPy is in sys.path
    if sfppy_root not in sys.path:
        sys.path.insert(0, sfppy_root)

    # Change working directory to the script's location
    script_dir = os.path.dirname(script_path)
    os.chdir(script_dir)

    # Execute the script
    script_globals = {}
    with open(script_path, "r", encoding="utf-8") as f:
        exec(f.read(), script_globals)

if __name__ == "__main__":
    main()
