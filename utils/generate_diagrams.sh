#!/usr/bin/env bash

# Generate a global Markdown documentation with Mermaid diagrams
# for all modules found in $mainfolder/patankar
#
# Prerequisites:
#  - Python 3.x
#  - generate_mermaid.py (the Python script) must be in the same directory
#  - generate_examples.py (the Python script) must be in the same directory
#
# Usage:
#  ./generate_diagrams.sh

# INRAE\\Olivier Vitrac
# Email: olivier.vitrac@agroparistech.fr
# Last Revised: 2025-01-17

# ------------------------
# Ensure Script is Run in Proper Context
# ------------------------

# Check if the '../patankar' directory exists
if [[ ! -d "../patankar" ]]; then
    echo "Error: '../patankar' directory not found. Adjust mainfolder path accordingly."
    exit 1
fi

# ------------------------
# Configuration Section
# ------------------------

# Main folder and paths
mainfolder="$(realpath ../)"
output_markdown="$mainfolder/html/SFPPy_classes_documentation.md"

# Python version (adjust if needed)
PYTHON_VERSION="3.10"

# Output directories and variables
output_dir="$(dirname "$output_markdown")"

# Ensure output directory exists
mkdir -p "$output_dir"

# Paths to include in PYTHONPATH
additional_paths=(
    "$mainfolder"
    "$mainfolder/patankar"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION/lib-dynload"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION/site-packages"
    "$HOME/.ipython"
)

# Set PYTHONPATH dynamically
export PYTHONPATH=$(IFS=:; echo "${additional_paths[*]}")
echo "PYTHONPATH set to: $PYTHONPATH"

# Directories to exclude (e.g., PIL)
excluded_dirs=(
    ".git"
    "cache.PubChem"
    "chemspipy"
    "pint"
    "demo"
    "examples"
    "debug"
    "sandbox"
    "windowsONLY"
    "fork"
    "tmp"
    "history"
)

# ------------------------
# Generate Modules List
# ------------------------

# Build the find command with exclusions
find_cmd="find \"$mainfolder/patankar\" -type f -name \"*.py\""
for dir in "${excluded_dirs[@]}"; do
    find_cmd+=" ! -path \"*$dir/*\""
done

# Execute the find command and transform paths for Python import
modules_list=$(eval "$find_cmd" | sed "s|$mainfolder/||;s|\.py$||" | sed 's|/|.|g')

# Check if any modules were found
if [[ -z "$modules_list" ]]; then
    echo "Error: No Python files found in the specified directory."
    exit 1
fi

# ------------------------
# Generate Class Examples and JSON Details
# ------------------------

# Define filenames for examples list and JSON details
moduleexamplesList="modules_withexamples_list.txt"
moduleexamplesDetails="class_examples_details.json"
moduleexamplesHTML="class_examples.html"
lookfolders=(
    "$mainfolder/patankar"
    "$mainfolder/patankar/private"
    #"$mainfolder/patankar/converted"
)

# Find .py files and generate the modules list
find "${lookfolders[@]}" -maxdepth 1 -type f -name '*.py' | sed "s|$mainfolder/||" | sort > "$output_dir/$moduleexamplesList"

# Verify if the modules list was created successfully
if [[ -f "$output_dir/$moduleexamplesList" ]]; then
    echo "Module examples list generated successfully at '$output_dir/$moduleexamplesList'"
else
    echo "Failed to generate module examples list."
    exit 1
fi

# Run generate_examples.py to create HTML and JSON files
./generate_examples.py "$output_dir/$moduleexamplesHTML" "$output_dir/$moduleexamplesDetails" < "$output_dir/$moduleexamplesList"

# Verify if both HTML and JSON files were created successfully
if [[ -f "$output_dir/$moduleexamplesHTML" && -f "$output_dir/$moduleexamplesDetails" ]]; then
    echo "Module examples page generated successfully at '$output_dir/$moduleexamplesHTML'"
    echo "Module details JSON generated successfully at '$output_dir/$moduleexamplesDetails'"
else
    echo "Failed to generate module examples page or details."
    exit 1
fi

# ------------------------
# Generate Documentation with Mermaid
# ------------------------

# Path to the JSON details file
mermaid_details_json="$output_dir/$moduleexamplesDetails"

# Ensure generate_mermaid.py is executable
chmod +x generate_mermaid.py

# Run generate_mermaid.py with the JSON details
echo "$modules_list" | python3 generate_mermaid.py "$output_markdown" "$mermaid_details_json"

# Check if the Markdown documentation was created successfully
if [[ -f "$output_markdown" ]]; then
    echo "Markdown documentation generated successfully at '$output_markdown'"
else
    echo "Failed to generate Markdown documentation."
    exit 1
fi

# ------------------------
# Final Message
# ------------------------

echo "Documentation (Mermaid+Markdown) generation completed successfully."
