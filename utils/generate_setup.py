#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates setup.py for SETUP

The following project structure is assumed.

SFPPy/
│
├── utils/
│   ├── generate_requirements.py  # SETUP script
│   ├── generate_manifest_in.py   # SETUP script
│   └── generate_setup.py         # SETUP script
│
├── patankar/
│   ├── __init__.py (if any, by default none)
│   ├── private/
│   │   ├── __init__.py (if any, by default none)
│   │   ├── chemspipy/
│   │       ├── __init__.py (if it exists)
│   │   └── pint/
│   │       ├── __init__.py (if it exists)
│   └── ... (other modules)
│
├── example1.py
├── example2.py
├── example3.py
├── tmp/
├── README.md
├── LICENSE
├── SFPPy.simple.manifest
├── requirements.txt             # run ./generate_requirements.py   from utils/
├── MANIFEST.in                  # run ./generate_manifest_in.py  from utils/
└── setup.py                     # run ./generate_setup.py   fropm utils/


Update dependencies to needs and future evolutions:
    dependencies = [
        "numpy>=1.21.0"
        "openpyxl >= 3.0.10"
    ]

Author:
    INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2025-02-12

    """


import os
import sys
import re
from pathlib import Path

# Define dependencies (same for pip and Conda)
dependencies = [
    "numpy>=1.21.0",      # Matrix computations (finite-volume method)
    "matplotlib>=3.4.0",  # Modern plotting functions
    "scipy>=1.7.0",       # Numerical integration and solvers
    "pandas>=1.3.0",      # Data manipulation and Excel support
    "openpyxl>=3.0.10"    # Save to .xlsx format
]

conda_channels = ["conda-forge", "defaults"]  # Recommended Conda channels


def is_utils_directory(current_path):
    """Verify that the script is run from the 'utils/' directory."""
    return current_path.name == 'utils'


def get_version(parent_dir):
    """Extract the version number of SFPPy from VERSION.txt."""
    version_file = parent_dir / "utils" / "VERSION.txt"
    if not version_file.exists():
        sys.stderr.write(f"Error: {version_file} not found. Please create a file with content: version=\"XX.YY.ZZ\"\n")
        sys.exit(1)
    
    with open(version_file, "r") as f:
        for line in f:
            match = re.match(r'^version\s*=\s*"(.*?)"$', line.strip())
            if match:
                return match.group(1)

    sys.stderr.write(f"Error: No valid version string found in {version_file}. Ensure it contains: version=\"XX.YY.ZZ\"\n")
    sys.exit(1)


def prompt_overwrite(file_path):
    """Prompt the user to overwrite an existing file."""
    while True:
        choice = input(f"'{file_path}' already exists. Overwrite? [Y/n]: ").strip().lower()
        if choice in ['', 'y', 'yes']:
            return True
        elif choice in ['n', 'no']:
            return False
        else:
            print("Please enter 'Y' or 'N'.")


def generate_setup_py(parent_dir, dependencies):
    """Generate setup.py with specified dependencies."""
    setup_path = parent_dir / "setup.py"

    if setup_path.exists() and not prompt_overwrite(setup_path):
        print("Skipping setup.py generation.")
        return

    setup_content = f"""from setuptools import setup, find_packages

setup(
    name="SFPPy",
    version="{get_version(parent_dir)}",
    description="Software Simulating Mass Transfer from Food Packaging",
    author="Olivier Vitrac",
    author_email="olivier.vitrac@agroparistech.fr",
    url="https://github.com/ovitrac/SFPPy",
    packages=find_packages(include=['patankar', 'patankar.*']),
    install_requires=[
        {', '.join([f'"{dep}"' for dep in dependencies])}
    ],
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    include_package_data=True,
    zip_safe=True,
)
"""

    with open(setup_path, 'w') as f:
        f.write(setup_content)
    
    print(f"✔ setup.py created successfully in '{parent_dir}'.")


def generate_environment_yml(parent_dir, dependencies, conda_channels):
    """Generate environment.yml for Conda users."""
    env_path = parent_dir / "environment.yml"

    if env_path.exists() and not prompt_overwrite(env_path):
        print("Skipping environment.yml generation.")
        return

    # Convert dependencies for Conda (removes version constraints if necessary)
    conda_deps = [dep.replace(">=", "=") for dep in dependencies]

    env_content = f"""name: sfppy
    channels:
    {chr(10).join(f"  - {channel}" for channel in conda_channels)}
    dependencies:
    - python=3.10
    {chr(10).join(f"  - {dep}" for dep in conda_deps)}
    """


    with open(env_path, 'w') as f:
        f.write(env_content)

    print(f"✔ environment.yml created successfully in '{parent_dir}'.")


def main():
    """Main script execution."""
    current_dir = Path.cwd()

    # Ensure the script is run from 'utils/'
    if not is_utils_directory(current_dir):
        print("Error: This script must be run from the 'utils/' directory of the SFPPy project.")
        sys.exit(1)

    # Define parent project directory
    parent_dir = current_dir.parent

    # Generate setup.py
    generate_setup_py(parent_dir, dependencies)

    # Generate environment.yml
    generate_environment_yml(parent_dir, dependencies, conda_channels)


if __name__ == '__main__':
    main()
