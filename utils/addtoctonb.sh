#!/bin/bash
# addtoctonb.sh
# This script adds a TOC to all notebook HTML files located in $mainfolder/../notebooks/export/
# by converting them with autotoc.py (located in $mainfolder/utils/) and saving the result in
# $mainfolder/../html/notebooks/ with the same filename.
#
# Maintained by INRAE\olivier.vitrac@agroparistech.fr
# Revision history: 2024-12-20

# -----------------------------------
#           INITIALIZATION
# -----------------------------------

# Ensure the script is run from the correct directory (SFPPy/utils/)
if [[ ! -f "pdocme.sh" ]]; then
    echo "Error: This script must be run from the SFPPy/utils/ directory."
    exit 1
fi

# Identify mainfolder dynamically
mainfolder="$(realpath ../)"

# Read __version__ from VERSION.txt
version_file="$mainfolder/utils/VERSION.txt"
if [[ ! -f "$version_file" ]]; then
  echo "Error: $version_file not found. Please create a file with content: version=\"XX.YY.ZZ\"" >&2
  exit 1
fi
__version__=$(grep -m 1 '^version=' "$version_file" | sed -E 's/version\s*=\s*"([^"]+)"/\1/')
if [[ -z "$__version__" ]]; then
  echo "Error: No valid version string found in $version_file. Ensure it contains: version=\"XX.YY.ZZ\"" >&2
  exit 1
fi
echo "SFPPy Version: $__version__"

# Set directories for input and output
input_dir="$mainfolder/notebooks/export"
output_dir="$mainfolder/wiki/notebooks"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# -----------------------------------
#           CONVERSION LOOP
# -----------------------------------

# Loop through each HTML file in the input directory
for file in "$input_dir"/*.html; do
    if [[ -f "$file" ]]; then
        base_name=$(basename "$file")
        output_file="$output_dir/$base_name"
        # Convert the HTML file using autotoc.py
        python3 "$mainfolder/utils/autotoc.py" "$file" "$output_file"
        if [[ $? -eq 0 ]]; then
            # Get the file size in bytes
            size=$(stat -c%s "$output_file")
            echo "Converted $base_name: ${size} bytes"
        else
            echo "Error converting $base_name"
        fi
    fi
done
