#!/usr/bin/env bash

###############################################################################
#                                                                             #
#    create_Pizza3_manifest.sh                                                #
#                                                                             #
#    Script to create the default manifest for the Pizza3 project              #
#                                                                             #
#    Maintained by INRAE\olivier.vitrac@agroparistech.fr                       #
#    Revision history: 2024-12-10                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

# Expected file structure
# Pizza3/
# ├── utils/
# │   ├── pdocme.sh
# │   └── manifestManager.py
# ├── pizza/
# │   ├── __init__.py
# │   ├── other_python_files.py
# │   └── ... (other directories and files)
# └── html/
#     ├── index.html
#     ├── pizza/
#     │   ├── __init__.html
#     │   ├── other_python_files.html
#     │   └── ... (other HTML files)
#     └── Pizza3.manifest


# Exit immediately if a command exits with a non-zero status
set -e

# Function to display error messages
error_exit() {
    echo "Error: $1"
    exit 1
}

# Identify the main project directory (parent directory)
mainfolder="$(realpath ../)"

# Verify that mainfolder exists and is a directory
[[ -d "$mainfolder" ]] || error_exit "Main project directory '$mainfolder' does not exist."

# Path to generate_manifest.sh
generate_manifest_script="$(pwd)/generate_manifest.sh"

# Check if generate_manifest.sh exists and is executable
if [[ ! -x "$generate_manifest_script" ]]; then
    error_exit "generate_manifest.sh not found or not executable at '$generate_manifest_script'."
fi

# Define the manifest name
manifest_name="Pizza3.manifest"

# Define excluded directories and files
excluded_dirs=(
    ".git"
    "fork"
    "history"
    "help"
    "debug"
    "tmp"
    "PIL"
    "restore"
    "__all__"
    "windowsONLY"
)

excluded_files=(
    "__init__.py"
    "__main__.py"
    "manifest.py"
    "debug.py"
)

# Define file extensions to include
extensions=(
    ".py"
    ".sh"
    ".md"
    ".html"
)

# Create the manifest using generate_manifest.sh with exclusions and extensions
echo "Creating default manifest '$manifest_name' for Pizza3 project with exclusions and extensions..."
"$generate_manifest_script" create -p "$mainfolder" -m "$manifest_name" --exclude-dirs "${excluded_dirs[@]}" --exclude-files "${excluded_files[@]}" --extensions "${extensions[@]}"

# Verify that the manifest was created
manifest_path="$mainfolder/$manifest_name"
if [[ -f "$manifest_path" ]]; then
    echo "Manifest '$manifest_name' successfully created at '$manifest_path'."
else
    error_exit "Failed to create manifest '$manifest_name' at '$manifest_path'."
fi
