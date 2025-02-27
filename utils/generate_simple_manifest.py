#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate the SFPPy.simple.manifest

This script traverses the project directory structure to create a manifest file 
listing all relevant files, excluding specified directories and file types. 
It ensures execution from the 'utils/' directory and outputs the manifest 
to the designated main folder with a timestamp.

Usage:
    python3 manifest_generator.py

Outputs:
    - SFPPy.simple.manifest: A text file listing all relevant project files.

Configuration:
    - Output Folder: Specified as the parent directory of 'utils/'.
    - Output Filename: 'SFPPy.simple.manifest'.
    - Excluded Folders and File Types: Defined in the 'ignore' list.

Author:
    INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2025-02-03
"""

import os
import sys
import datetime

def main():
    # Define the manifest filename
    manifest_filename = "SFPPy.simple.manifest"

    # Define the list of file extensions and directories to ignore
    ignore_extensions = [
        '.manifest', '.py~', '.sh~', '.pyc',
        '.sample', '.cache', '.xml', '.iml', '.zip', '.js',
        '.pdf', '.png', '.mp4', '.avi', '.html~', '.gitignore~',
        '.gif','.pptx','.egg-info'
    ]
    
    ignore_dirs = [
        ".env",
        ".venv",
        ".git",
        "fork",
        "history",
        "help",
        "debug",
        "sandbox",
        "draft",
        "src",
        "tmp",
        "trash",
        "PIL.egg-info",
        "restore",
        "__all__",
        "windowsONLY",
        "obsolete",
        ".spyproject",
        ".vscode",
        "videos",
        "converted",
        "release",
        "first_doc",
        "old-doc",
        "FIles_moved_to_Post",
        ".ipynb_checkpoints",
        ".virtual_documents",
        "books",
        "external",
        "tmp",
        "html",
        ".spyproject",
        ".ipynb_checkpoints"
    ]

    # Determine the current working directory
    current_dir = os.getcwd()

    # Verify that the script is executed from the 'utils/' directory
    if os.path.basename(current_dir) != "utils":
        print("Error: This script must be executed from the 'utils/' directory.")
        sys.exit(1)

    # Define the main folder as the parent directory of 'utils/'
    mainfolder = os.path.abspath(os.path.join(current_dir, os.pardir))

    # Define the full path for the manifest file (placed inside mainfolder)
    manifest_path = os.path.join(mainfolder, manifest_filename)

    try:
        with open(manifest_path, 'w', encoding='utf-8') as manifest_file:
            # Write the header with the date of generation
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            manifest_file.write("# SFPPy.simple.manifest\n")
            manifest_file.write(f"# Generated on: {current_time}\n\n")
            
            # Traverse the directory structure
            for root, dirs, files in os.walk(mainfolder):
                # Exclude specified directories
                dirs[:] = [d for d in dirs if d not in ignore_dirs]
                
                for name in files:
                    file_path = os.path.join(root, name)
                    _, ext = os.path.splitext(file_path)
                    
                    # Exclude files with ignored extensions
                    if ext in ignore_extensions:
                        continue
                    
                    # Write the relative file path to the manifest
                    relative_path = os.path.relpath(file_path, mainfolder)
                    manifest_file.write(f"{relative_path}\n")
        
        print(f"Manifest generated successfully at: {manifest_path}")
    
    except Exception as e:
        print(f"An error occurred while generating the manifest: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

