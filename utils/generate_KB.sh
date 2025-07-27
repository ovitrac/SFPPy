#!/bin/bash
# generate_KB.sh
# This script creates a Markdown knowledge base in $mainfolder/docs/KB by collecting all .md files
# under $mainfolder, preserving subfolder structure and excluding backup/history/temp folders.
#
# Maintained by INRAE\olivier.vitrac@agroparistech.fr
# Revision history: 2025-07-25 (array-based folder exclusion)

# -----------------------------------
#         INITIALIZATION
# -----------------------------------

if [[ ! -f "pdocme.sh" ]]; then
    echo "Error: This script must be run from the SFPPy/utils/ directory."
    exit 1
fi

mainfolder="$(realpath ../)"
srcdir="$mainfolder"
kbdir="$mainfolder/docs/KB"

mkdir -p "$kbdir"

# -----------------------------------
#         CONFIGURATION
# -----------------------------------

# Define excluded folders (case-insensitive, name-based match)
excluded_folders=("history" "backup" "trash" "temp" "tmp" "sandbox" "release" "src" "books" ".*")

# Construct find exclusion expression
find_exclude=()
for dir in "${excluded_folders[@]}"; do
  find_exclude+=(-iname "$dir" -o)
done
# Remove last -o
unset 'find_exclude[${#find_exclude[@]}-1]'

# -----------------------------------
#         FILE COPYING
# -----------------------------------

echo "Creating knowledge base in: $kbdir"
echo "Searching for .md files under: $srcdir"
echo "Excluding folders: ${excluded_folders[*]}"

find "$srcdir" \
    -type d \( "${find_exclude[@]}" \) -prune -o \
    -type f -name "*.md" ! -path "$kbdir/*" -print | while read -r filepath; do
        relpath="${filepath#$srcdir/}"
        targetpath="$kbdir/$relpath"
        targetdir="$(dirname "$targetpath")"
        mkdir -p "$targetdir"
        cp "$filepath" "$targetpath"
        echo "Copied: $relpath"
done

echo "Knowledge base generation complete."

