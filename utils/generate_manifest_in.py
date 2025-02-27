#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates MANIFEST.in for setup.py from Pizza3.simple.manifest

The following project structure is assumed.

Pizza3/
│
├── utils/
│   ├── generate_requirements.py  # SETUP script
│   ├── generate_manifest_in.py   # SETUP script
│   └── generate_setup.py         # SETUP script
│
├── pizza/
│   ├── __init__.py
│   ├── private/
│   │   ├── __init__.py
│   │   └── PIL/
│   │       ├── __init__.py
│   │       └── ... (other PIL modules)
│   └── ... (other modules)
│
├── example2.py
├── tmp/
├── README.md
├── LICENSE
├── Pizza3.simple.manifest
├── requirements.txt             # run ./generate_requirements.py   from utils/
├── MANIFEST.in                  # run ./generate_manifest_in.py  from utils/
└── setup.py                     # run ./generate_setup.py   from utils/


Generating MANIFEST.in from Pizza3.simple.manifest
The MANIFEST.in file tells setuptools about the non-Python files to include in Pizza3 package. 

Author:
    INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2024-12-12

    """

import os
import sys
from pathlib import Path

def is_utils_directory(current_path):
    """
    Verify that the script is run from the 'utils/' directory.
    """
    return current_path.name == 'utils'

def prompt_overwrite(file_path):
    """
    Prompt the user to overwrite an existing MANIFEST.in file.
    """
    while True:
        choice = input(f"'{file_path}' already exists. Overwrite? [Y/n]: ").strip().lower()
        if choice in ['', 'y', 'yes']:
            return True
        elif choice in ['n', 'no']:
            return False
        else:
            print("Please enter 'Y' or 'N'.")

def generate_manifest_in(parent_dir, manifest_file, output_file):
    """
    Generate MANIFEST.in based on the contents of Pizza3.simple.manifest.
    """
    manifest_path = parent_dir / manifest_file
    manifest_in_path = parent_dir / output_file

    if not manifest_path.exists():
        print(f"Error: '{manifest_file}' does not exist in '{parent_dir}'.")
        sys.exit(1)

    if manifest_in_path.exists():
        if not prompt_overwrite(manifest_in_path):
            print("Operation cancelled. Existing 'MANIFEST.in' was not modified.")
            sys.exit(0)

    try:
        with open(manifest_path, 'r') as mf:
            lines = mf.readlines()

        with open(manifest_in_path, 'w') as mif:
            mif.write("# MANIFEST.in\n")
            mif.write("# Automatically generated from Pizza3.simple.manifest\n\n")
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip empty lines and comments
                # For each file, determine if it's inside a package
                if '/' in line:
                    # Assuming data files are inside packages or subdirectories
                    mif.write(f"recursive-include {line} *\n")
                else:
                    # Files in the root directory
                    mif.write(f"include {line}\n")

        print(f"'{output_file}' has been generated successfully in '{parent_dir}'.")
    except Exception as e:
        print(f"An error occurred while generating '{output_file}': {e}")
        sys.exit(1)

def main():
    # Get the current working directory
    current_dir = Path.cwd()

    # Ensure the script is run from the 'utils/' directory
    if not is_utils_directory(current_dir):
        print("Error: This script must be run from the 'utils/' directory of the Pizza3 project.")
        sys.exit(1)

    # Define the parent (root) directory
    parent_dir = current_dir.parent

    # Define the manifest file and the output MANIFEST.in
    manifest_file = 'Pizza3.simple.manifest'
    output_file = 'MANIFEST.in'

    # Generate MANIFEST.in
    generate_manifest_in(parent_dir, manifest_file, output_file)

if __name__ == '__main__':
    main()
