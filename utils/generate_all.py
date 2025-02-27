#!/usr/bin/env python3

"""
Script to dynamically generate and update __all__ variables in Python modules.
Restricts __all__ to symbols defined within the module itself.
Generates backup files with __all__ definitions in utils/__all__/ for review.
"""

# Maintained by INRAE\\olivier.vitrac@agroparistech.fr
# Revision history: 2024-12-08
# 2025-01-07 fix extract_symbols by adding package_root
# 2025-01-17 add SFPPy.converted modules

import os
import inspect
import pkgutil
import sys

def extract_symbols(module_name,package_root):
    """
    Extracts all public symbols (classes, functions, variables) from a module,
    restricting to symbols defined in the module itself.
    """
    try:
        module = __import__(module_name, fromlist=["*"])
        symbols = []
        for name, obj in inspect.getmembers(module):
            # Skip private symbols and ensure the object is defined in the current module
            if name.startswith("_"):
                continue

            # Include only symbols where the __module__ attribute matches the module_name
            origin = getattr(obj, "__module__", None)
            if origin is not None and origin.startswith(package_root): # origin == module_name:
                symbols.append(name)

        return symbols
    except ImportError as e:
        print(f"Error importing module {module_name}: {e}")
        return []

def write_backup_all(module_name, symbols, all_dir):
    """
    Generates a backup file for __all__ in the utils/__all__/ directory.
    """
    # Prepare the backup directory
    if not os.path.exists(all_dir):
        os.makedirs(all_dir)

    # Write the __all__ variable to the backup file
    backup_path = os.path.join(all_dir, f"{module_name.split('.')[-1]}_all.py")
    with open(backup_path, "w") as f:
        f.write(f"# __all__ for {module_name}\n")
        f.write(f"__all__ = {symbols!r}\n")

    print(f"Generated backup __all__ for {module_name}: {backup_path}")

def update_module_with_all(module_name, symbols, base_path, package_root):
    """
    Updates the `__all__` variable in the target module.
    Adds or replaces the `__all__` definition in the module file.
    """
    # Resolve the correct path for the module file
    relative_path = module_name.replace(package_root + ".", "").replace(".", "/") + ".py"
    module_path = os.path.join(base_path, relative_path)
    if not os.path.exists(module_path):
        print(f"Module file not found: {module_path}")
        return

    with open(module_path, "r") as f:
        lines = f.readlines()

    # Detect encoding comment (e.g., # -*- coding: utf-8 -*-)
    encoding_line_index = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("# -*- coding:"):
            encoding_line_index = i
            break

    # Detect imports
    import_end_index = encoding_line_index + 1 if encoding_line_index >= 0 else 0
    for i, line in enumerate(lines[import_end_index:], start=import_end_index):
        if not line.startswith("import") and not line.startswith("from") and line.strip():
            break
        import_end_index = i + 1

    # Check if __all__ already exists
    updated = False
    for i, line in enumerate(lines):
        if line.strip().startswith("__all__"):
            lines[i] = f"__all__ = {symbols!r}\n"
            updated = True
            break

    # Insert __all__ after imports (or encoding line if no imports)
    if not updated:
        insert_index = import_end_index
        lines.insert(insert_index, f"\n__all__ = {symbols!r}\n")

    # Write back the updated module
    with open(module_path, "w") as f:
        f.writelines(lines)

    print(f"Updated {module_name} with __all__ = {symbols!r}")

def generate_all_for_package(package_name, base_path, package_root, all_dir):
    """
    Generates and updates the `__all__` variable for all modules in the given package.
    """
    try:
        package = __import__(package_name, fromlist=["*"])
        for loader, module_name, is_pkg in pkgutil.iter_modules(package.__path__):
            if is_pkg:
                continue  # Skip sub-packages
            module_full_name = f"{package_name}.{module_name}"
            print(f"Processing module: {module_full_name}")

            symbols = extract_symbols(module_full_name,package_root)
            if symbols:
                # Write a backup __all__ file in utils/__all__/
                write_backup_all(module_full_name, symbols, all_dir)

                # Update the module file with the __all__ variable
                update_module_with_all(module_full_name, symbols, base_path, package_root)
            else:
                print(f"No public symbols found in {module_full_name}")
    except ImportError as e:
        print(f"Error importing package {package_name}: {e}")

if __name__ == "__main__":
    # Ensure the script is run from utils/
    script_dir = os.path.dirname(os.path.realpath(__file__))
    if not os.path.basename(script_dir) == "utils":
        print("Error: This script must be run from the SFPPy/utils/ directory.")
        sys.exit(1)

    # Adjust the paths
    mainfolder = os.path.abspath(os.path.join(script_dir, "../"))
    SFPPy_base_path = os.path.join(mainfolder, "patankar")  # Base path for SFPPy package
    all_dir = os.path.join(script_dir, "__all__")  # Directory for backup __all__ files
    sys.path.insert(0, mainfolder)  # Add mainfolder as root
    sys.path.insert(0, SFPPy_base_path)  # Add SFPPy directory as package root

    # Packages to process
    packages = ["patankar", "patankar.private", "patankar.converted"]

    # Generate __all__ for each package
    for package in packages:
        generate_all_for_package(package, SFPPy_base_path, "patankar", all_dir)

    print("\n__all__ generation and updates completed.")
