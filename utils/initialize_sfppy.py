#!/usr/bin/env python
"""
Temporary Environment Setup for SFPPy (initialize_sfppy.py)

This script is intended for use in read-only environments (e.g., Colab) or
launched from index.ipynb/command line. It temporarily sets up the SFPPy
environment as follows:

1. Determines the SFPPy root folder from its own location (i.e., the parent folder
   of the utils directory).
2. Checks if that folder is already in sys.path; if so, it prints an informative message.
3. If not, it adds the SFPPy root to sys.path and, optionally, changes the current working directory
   to that folder (controlled by the flag do_cd).
4. It then prints instructions for a full installation using pip or conda.  
   These instructions remind the user to **cd** into the SFPPy root folder before running:
      - pip install commands (e.g., "pip install -e .")
      - conda commands (e.g., "conda env create -f environment.yml")
5. Finally, it suggests running an example (e.g., example1.py) or launching index.ipynb
   in a Jupyter Notebook or JupyterLab environment.

Usage:
    Run this script from index.ipynb or directly from the command line.
    
Flags:
    --no-verbose  : Turn off verbose output (default is verbose output).
    --no-cd       : Prevent changing the working directory (default is to change it).

**INRAE\Olivier Vitrac**  
Email: [olivier.vitrac@agroparistech.fr](mailto:olivier.vitrac@agroparistech.fr)
    
"""

import os
import sys
import argparse

def main(verbosity=True, do_cd=True):
    # Determine the directory of this script (utils/)
    utils_dir = os.path.dirname(os.path.abspath(__file__))
    # The SFPPy root folder is the parent of the utils folder
    main_folder = os.path.dirname(utils_dir)
    
    if verbosity:
        print("== SFPPy Temporary Environment Setup ==")
        print(f"Detected SFPPy root folder: {main_folder}\n")
    
    # Check if the SFPPy root is already in sys.path
    if main_folder in sys.path:
        if verbosity:
            print("INFO: SFPPy is already in sys.path.")
            print(f"Current sys.path entry: {main_folder}\n")
    else:
        if verbosity:
            print("INFO: SFPPy is not in sys.path. Adding it now...")
        sys.path.insert(0, main_folder)
        if do_cd:
            try:
                os.chdir(main_folder)
                if verbosity:
                    print(f"Changed current working directory to: {main_folder}")
            except Exception as e:
                if verbosity:
                    print(f"WARNING: Could not change directory to {main_folder}: {e}")
        else:
            if verbosity:
                print("Skipping changing directory due to do_cd flag.")
        if verbosity:
            print(f"Added {main_folder} to sys.path.\n")
    
    if verbosity:
        print("== Full Installation Instructions ==")
        print("To fully install SFPPy (for persistent and standard use), please follow these steps:")
        print("1. Make sure you are in the SFPPy root folder (current directory should be:)")
        print(f"      {main_folder}")
        print("   If you are not, run:")
        print("      cd <path_to_SFPPy_root>\n")
        print("2. For pip users, install SFPPy in editable mode:")
        print("      pip install -e .\n")
        print("3. For conda users, create the environment using the provided environment.yml:")
        print("      conda env create -f environment.yml")
        print("   (Then activate the environment with 'conda activate SFPPy')\n")
        print("4. Once installed, you can run an example, e.g.:")
        print("      python example1.py")
        print("   Or launch the index.ipynb in a Jupyter Notebook/JupyterLab environment.\n")
        print("Thank you for using SFPPy!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-verbose', action='store_false', dest='verbosity',
                        default=True, help="Turn off verbose output")
    parser.add_argument('--no-cd', action='store_false', dest='do_cd',
                        default=True, help="Prevent changing the working directory")
    args = parser.parse_args()
    main(verbosity=args.verbosity, do_cd=args.do_cd)