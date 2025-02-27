#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates requirements.txt for SETUP

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


Update dependencies to needs and future evolutions:
    dependencies = [
        "numpy>=1.21.0"
    ]

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

def check_project_structure(parent_dir):
    """
    Check if the parent directory contains essential project files.
    Modify the expected_files list based on your project's actual structure.
    """
    expected_files = ['README.md', 'setup.py']  # Customize as needed
    for file in expected_files:
        if not (parent_dir / file).exists():
            return False
    return True

def prompt_overwrite(file_path):
    """
    Prompt the user to overwrite an existing requirements.txt file.
    """
    while True:
        choice = input(f"'{file_path}' already exists. Overwrite? [Y/n]: ").strip().lower()
        if choice in ['', 'y', 'yes']:
            return True
        elif choice in ['n', 'no']:
            return False
        else:
            print("Please enter 'Y' or 'N'.")

def generate_requirements(parent_dir, output_file, dependencies):
    """
    Generate a requirements.txt file with the specified dependencies.
    """
    requirements_path = parent_dir / output_file

    if requirements_path.exists():
        if not prompt_overwrite(requirements_path):
            print("Operation cancelled. Existing 'requirements.txt' was not modified.")
            sys.exit(0)

    try:
        with open(requirements_path, 'w') as f:
            if dependencies:
                for dep in dependencies:
                    f.write(f"{dep}\n")
                print(f"'{output_file}' has been created with specified dependencies in '{parent_dir}'.")
            else:
                # If no external dependencies, create an empty requirements.txt or add a comment
                f.write("# No external dependencies required.\n")
                print(f"'{output_file}' has been created as an empty file in '{parent_dir}'.")
    except Exception as e:
        print(f"An error occurred while writing '{output_file}': {e}")
        sys.exit(1)

def main():
    # Get the current working directory
    current_dir = Path.cwd()

    # Ensure the script is run from the 'utils/' directory
    if not is_utils_directory(current_dir):
        print("Error: This script must be run from the 'utils/' directory of the Pizza project.")
        sys.exit(1)

    # Define the parent (root) directory
    parent_dir = current_dir.parent

    # Verify project structure (optional but recommended)
    if not check_project_structure(parent_dir):
        print(f"Error: The parent directory '{parent_dir}' does not contain the expected project files.")
        print("Please ensure you're running the script from the correct 'utils/' directory.")
        sys.exit(1)

    # Define the output file
    output_file = 'requirements.txt'

    # Define external dependencies
    # List all external packages your project depends on.
    # Since your project primarily uses standard libraries and an internal PIL,
    # you might not need to list any external dependencies.
    dependencies = [
        "numpy>=1.21.0"
    ]

    # Generate requirements.txt
    generate_requirements(parent_dir, output_file, dependencies)

if __name__ == '__main__':
    main()
