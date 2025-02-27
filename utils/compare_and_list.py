#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
    compare_and_list.py

----------------------------------------------------------------------
This script compares files in the local directory structure (source)
against a destination folder (remote copy). It checks for:
  - missing files
  - "obsolete" files (based on file size, last modification date, or both)

NEW FEATURES:
  - For MISSING files, display (YYYY-MM-DD) and size of the source file.
  - Third argument "missing" => only show missing files
                "update"  => show missing and obsolete files (default).

Usage:
  ./compare_and_list.py <destination_folder> [comparison_mode] [operation_mode] [create_missing_folders_flag]

Examples:
  ./compare_and_list.py /path/to/remote  # Default: comparison_mode=date, operation_mode=update, no folder creation
  ./compare_and_list.py /path/to/remote size missing yes
    => Compare by size only; display only missing files; create missing folders.

Production Examples:
    /compare_and_list.py $HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3

OUTPUT FORMAT (examples):
  (size diff: +2.3 KB, time diff: +45.0 min) OBSOLETE: /path/to/source -> /path/to/dest
  (source: 2025-01-05, 12 KB) MISSING: /path/to/source -> /path/to/dest

----------------------------------------------------------------------

    # Documentation

    ## Overview

    This script, *compare_and_list.py*, inspects a local *Pizza3* project directory (which we call the **source**) and a remote copy (which we call the **destination**). It reports files that are:
    - **Missing** in the destination.
    - **Obsolete** in the destination, based on one of three modes:
    - **date**: The source file has a more recent modification date.
    - **size**: The source file is larger than the destination file.
    - **both**: Either the source is newer or the source is larger.

    By default, *compare_and_list.py* compares using the **date** mode.

    It also checks the **inclusion** and **exclusion** patterns defined in the script, mirroring a backup configuration. Only files that match an inclusion pattern and do **not** match an exclusion pattern will be considered for comparison.

    ## Prerequisites

    1. **Script location**: The script must reside in `$mainfolder/utils/` and must be run from there.
    2. **Name of the current folder**: The script checks that the current directory’s name is `utils`. If not, it fails.
    3. **Presence of pdocme.sh**: The script requires a file named `pdocme.sh` in the `utils/` folder to confirm the environment is correct.
    4. **Destination folder**: The script requires a mandatory first argument indicating the destination folder. It checks that this folder has a `utils/` subfolder containing `pdocme.sh`, ensuring it is a remote copy of the same structure.

    ## Command-Line Arguments

    The script takes up to four arguments:

    1. **destination_folder** (mandatory)  
    The path to the remote/copy of the main folder. If relative, the script automatically converts it to an absolute path.  
    Example: `$HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3`

    2. **comparison_mode** (optional)  
    Specifies the comparison mode:  
    - **date** : Compare the modification times (source newer => replace needed).  
    - **size** : Compare the sizes (source bigger => replace needed).  
    - **both** : Compare both criteria (if source is either newer or bigger => replace needed).  
    If not provided, the default is **date**.

    3. **operation_mode** (optional)  
    - **missing**: Only list missing files.  
    - **update** (default): List missing and obsolete.  

    4. **create_missing_folders_flag** (optional)  
    - **yes** : If set to 'yes', the script will automatically create any missing folders in the destination when a missing file is found.  
    - **no**  : No folder creation if the file’s destination path doesn’t exist.  
    The default is **no**.

    ## Usage Examples

    ### Example 1: Compare only by date, no folder creation

    ```bash
    cd $HOME/han/dev/Pizza3/utils/
    ./compare_and_list.py $HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3
    ```

    1. The script changes into the *utils* directory:  
    *$HOME/han/dev/Pizza3/utils/*
    2. Runs the script with only one argument: the destination folder.  
    3. Uses the default comparison mode (**date**).  
    4. Does **not** create missing folders.

    Any files missing or “obsolete by date” in the destination are listed in the output, which can be redirected as needed, for example:

    ```bash
    ./compare_and_list.py $HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3 > missing_or_obsolete.txt
    ```

    ### Example 2: Compare by size and create missing folders

    ```bash
    cd $HOME/han/dev/Pizza3/utils/
    ./compare_and_list.py $HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3 size yes
    ```

    1. The script again is run from the same local directory (*utils*).
    2. **destination_folder** is `$HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3`.  
    3. **comparison_mode** is `size`.  
    4. **create_missing_folders_flag** is `yes`, so any required subfolders that do not exist in the destination will be created.

    The script outputs lines like:

    - `MISSING: source_file_path -> destination_file_path`  
    - `OBSOLETE: source_file_path -> destination_file_path`  

    You can filter or redirect this output to logs or another program.

    ## Explanation of Output Lines

    - **MISSING:** A file was found in the source but is absent in the destination.  
    - **OBSOLETE:** A file exists in both locations, but the source is considered “newer” or “larger” (depending on the mode), meaning a replacement is appropriate.

    Each line contains the **full path** of the source file and the **corresponding path** in the destination. This format is convenient for other scripts or tooling to parse and handle necessary copies.

    ## Inclusion/Exclusion Patterns

    - **include_patterns**: List of filename patterns (wildcards) that must match for the file to be considered (e.g., `*.py`, `*.m`, `*.txt`, etc.).  
    - **exclude_files**: Specific file patterns to exclude (e.g., `*.log`, `backupme.README.md`, etc.).  
    - **exclude_folders_rel**: Relative folders to exclude (e.g., `./old`, `./tmp`).  
    - **exclude_folders_abs**: Absolute folders to exclude, relevant to `$mainfolder` (e.g., `$mainfolder/release`).  

    The script automatically ignores files and folders that match any exclusion pattern.

    ## Logging and Creation of Folders

    - When *create_missing_folders_flag* is set to **yes**, the script attempts to create the destination folder structure if it is missing. Upon successful creation, it logs a line such as:  
    `Created missing folder: /path/to/new/folder`

    - In all cases, the script **never** modifies the source files or folders. It only prints actions that might be needed in the destination directory.

    ## Error Handling

    1. **Script not in utils**: The script will refuse to run if the current directory name is not “utils”.  
    2. **pdocme.sh missing**: The script checks for this file as a marker that `$mainfolder/utils` is correct.  
    3. **Destination folder not found**: The script exits with an error if the user-supplied destination folder path does not exist.  
    4. **Invalid remote copy**: If the destination path lacks `utils/pdocme.sh`, the script fails.  
    5. **Invalid comparison mode**: If the user supplies an unrecognized mode, the script warns and reverts to the default “date” mode.

    --- 

**Author:**
---------
INRAE\Olivier Vitrac  
Email: olivier.vitrac@agroparistech.fr  
Last Revised: 2025-01-09
'''

import os
import sys
import fnmatch
import time
import datetime

# ----------------------------------------------------------------------
# Global Patterns for Inclusions and Exclusions
# ----------------------------------------------------------------------
include_patterns = [
    "*.m",       # MATLAB files
    "*.asv",     # Auto-saved files
    "*.m~",      # Backup files
    "*.pynb",    # Jupyter notebooks
    "*.py",      # Python scripts
    "*.sh",      # Shell scripts
    "*.txt",     # Text files
    "*.md",      # Markdown files
    "*.html",    # HTML files
    "*.json",    # JSON files
    "*.css",     # CSS files
    "*.manifest" # Manifest files
]

exclude_folders_rel = [
    "./old",
    "./tmp",
    "./sandbox",
    "./debug",
    "./obsolete",
    "./.git",
    "./.vscode",
    "./.spyproject",
    "./__all__",
    "./__pycache__"
]

exclude_folders_abs = [
    # e.g. "/absolute/path/to/history",
    # e.g. "/absolute/path/to/release",
]

exclude_files = [
    "*.log",
    "*.zip",
    "backupme.README.md"
]

# ----------------------------------------------------------------------
# Additional constants
# ----------------------------------------------------------------------
TIME_DIFF_TOLERANCE = 3600  # e.g., 1 hour

def is_included(filename):
    """
    Checks if a file should be included based on the include_patterns.
    """
    return any(fnmatch.fnmatch(filename, pat) for pat in include_patterns)

def is_excluded_file(filename):
    """
    Checks if a file is explicitly excluded based on exclude_files.
    """
    return any(fnmatch.fnmatch(filename, pat) for pat in exclude_files)

def is_excluded_folder(path, mainfolder):
    """
    Checks if a folder should be excluded based on:
      1) The relative folder exclusions (exclude_folders_rel)
      2) The absolute folder exclusions (exclude_folders_abs)
    """
    abs_path = os.path.abspath(path)
    
    # Check relative folder patterns
    for rel_excl in exclude_folders_rel:
        if path.endswith(rel_excl.lstrip("./")):
            return True
    
    # Check absolute folder patterns
    for abs_excl in exclude_folders_abs:
        if abs_path.startswith(abs_excl):
            return True
    
    return False

def usage():
    """
    Prints usage instructions and exits.
    """
    print("Usage:")
    print("  ./compare_and_list.py <destination_folder> [comparison_mode] [operation_mode] [create_missing_folders_flag]")
    print("")
    print("Where:")
    print("  <destination_folder> : Mandatory. Path to remote/copy of mainfolder.")
    print("  [comparison_mode]    : Optional. 'date', 'size', or 'both'. Default is 'date'.")
    print("  [operation_mode]     : Optional. 'missing' or 'update'. Default is 'update'.")
    print("      'missing' => only show missing files")
    print("      'update'  => show missing and obsolete files")
    print("  [create_missing_folders_flag] : Optional. 'yes' or 'no'. Default is 'no'.")
    sys.exit(1)

# ----------------------------------------------------------------------
# Helper functions to format size/time differences
# ----------------------------------------------------------------------
def format_size(bytes_val):
    """
    Return a short string for a file size in bytes, e.g. '12 B', '1.2 KB', '3.4 MB', etc.
    """
    if bytes_val < 1024:
        return f"{bytes_val} B"
    elif bytes_val < 1024**2:
        return f"{bytes_val/1024:.1f} KB"
    else:
        return f"{bytes_val/(1024**2):.1f} MB"

def format_size_diff(size_diff_bytes):
    sign = "+" if size_diff_bytes >= 0 else "-"
    abs_val = abs(size_diff_bytes)
    if abs_val < 1024:
        return f"{sign}{abs_val} B"
    elif abs_val < 1024**2:
        return f"{sign}{abs_val/1024:.1f} KB"
    else:
        return f"{sign}{abs_val/(1024**2):.1f} MB"

def format_time_diff(time_diff_seconds):
    sign = "+" if time_diff_seconds >= 0 else "-"
    abs_diff = abs(time_diff_seconds)
    # if < 60 seconds
    if abs_diff < 60:
        return f"{sign}{int(abs_diff)} s"
    elif abs_diff < 3600:
        mins = abs_diff / 60.0
        return f"{sign}{mins:.1f} min"
    elif abs_diff < 86400:
        hours = abs_diff / 3600.0
        return f"{sign}{hours:.1f} h"
    else:
        days = abs_diff / 86400.0
        return f"{sign}{days:.1f} days"

def format_date(timestamp):
    """
    Convert an epoch timestamp (float or int) to YYYY-MM-DD.
    """
    dt = datetime.datetime.fromtimestamp(timestamp)
    return dt.strftime("%Y-%m-%d")

def main():
    # -------------------------------------------------
    # Step 0: Preliminary checks on script location
    # -------------------------------------------------
    current_dir = os.path.basename(os.getcwd())
    if current_dir != "utils":
        print("Error: You must run this script from the 'utils' folder.")
        sys.exit(1)
    # Check presence of pdocme.sh
    if not os.path.exists("pdocme.sh"):
        print("Error: 'pdocme.sh' is not found in the current folder. Execution aborted.")
        sys.exit(1)
    
    # -------------------------------------------------
    # Step 1: Parse arguments
    # -------------------------------------------------
    if len(sys.argv) < 2:
        usage()
    destination_folder_input = sys.argv[1]
    comparison_mode = "date"  # default
    operation_mode = "update" # default
    create_missing_folders = False
    
    # Arg2: comparison_mode
    if len(sys.argv) >= 3:
        arg2 = sys.argv[2].lower()
        if arg2 in ["date", "size", "both"]:
            comparison_mode = arg2
        elif arg2 in ["missing", "update"]:
            # If user provided 'missing' or 'update' as second arg, treat that as operation_mode
            operation_mode = arg2
        else:
            # not recognized => keep default comparison_mode
            pass
    
    # Arg3: either operation_mode or create_missing_folders_flag
    if len(sys.argv) >= 4:
        arg3 = sys.argv[3].lower()
        if arg3 in ["missing", "update"]:
            operation_mode = arg3
        elif arg3 in ["yes", "no"]:
            create_missing_folders = (arg3 == "yes")
        else:
            # not recognized
            pass
    
    # Arg4: possibly leftover for create_missing_folders_flag
    if len(sys.argv) >= 5:
        arg4 = sys.argv[4].lower()
        if arg4 in ["yes", "no"]:
            create_missing_folders = (arg4 == "yes")
    
    # -------------------------------------------------
    # Step 2: Validate destination folder
    # -------------------------------------------------
    destination_folder = os.path.abspath(destination_folder_input)
    
    if not os.path.exists(destination_folder):
        print(f"Error: Destination folder '{destination_folder}' does not exist.")
        sys.exit(1)
    
    # Check that destination_folder includes 'utils' subfolder and 'pdocme.sh'
    utils_subfolder = os.path.join(destination_folder, "utils")
    pdocme_file = os.path.join(utils_subfolder, "pdocme.sh")
    if not os.path.exists(utils_subfolder) or not os.path.exists(pdocme_file):
        print("Error: Destination folder is not a valid remote copy of the mainfolder. "
              "It must contain 'utils/pdocme.sh'.")
        sys.exit(1)
    
    # -------------------------------------------------
    # Step 3: Walk through the source folder ($mainfolder)
    # -------------------------------------------------
    mainfolder = os.path.dirname(os.path.abspath(os.getcwd()))
    
    for root, dirs, files in os.walk(mainfolder):
        rel_path = os.path.relpath(root, mainfolder)
        if rel_path == ".":
            rel_path = ""
        
        if is_excluded_folder(root, mainfolder):
            dirs[:] = []
            continue
        
        for filename in files:
            source_file = os.path.join(root, filename)
            
            if not is_included(filename):
                continue
            if is_excluded_file(filename):
                continue
            
            dest_file = os.path.join(destination_folder, rel_path, filename)
            source_stat = os.stat(source_file)
            
            # Check if file is missing
            if not os.path.exists(dest_file):
                if operation_mode == "update" or operation_mode == "missing":
                    # We do display missing if mode is 'missing' or 'update'
                    # Show date + size of the source
                    source_date = format_date(source_stat.st_mtime)
                    source_size = format_size(source_stat.st_size)
                    msg = f"(source: {source_date}, {source_size}) MISSING: {source_file} -> {dest_file}"
                    print(msg)
                    if create_missing_folders:
                        dest_dir = os.path.dirname(dest_file)
                        if not os.path.exists(dest_dir):
                            os.makedirs(dest_dir, exist_ok=True)
                            print(f"Created missing folder: {dest_dir}")
                # No need to check obsolete logic because it's missing
                continue
            
            # If operation_mode == "missing", we skip checking obsolete
            if operation_mode == "missing":
                continue
            
            # operation_mode == "update": check if it is obsolete
            dest_stat = os.stat(dest_file)
            replace_required = False
            
            # Evaluate time difference
            time_diff = source_stat.st_mtime - dest_stat.st_mtime
            # Evaluate size difference
            size_diff = source_stat.st_size - dest_stat.st_size
            
            # Check date
            if comparison_mode in ["date", "both"]:
                if time_diff > TIME_DIFF_TOLERANCE:
                    replace_required = True
            
            # Check size
            if comparison_mode in ["size", "both"]:
                if source_stat.st_size > dest_stat.st_size:
                    replace_required = True
            
            if replace_required:
                size_diff_str = format_size_diff(size_diff)
                time_diff_str = format_time_diff(time_diff)
                print(f"({size_diff_str}, {time_diff_str}) OBSOLETE: {source_file} -> {dest_file}")
    
    print("Comparison completed.")

if __name__ == "__main__":
    main()
