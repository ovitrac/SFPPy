#!/bin/bash
# generate_release.sh
# This script creates a release package for Pizza3, based on a manifest file.
#
# Use ./generate_simple_manifest.py to refresh the manifest

# Maintained by INRAE\olivier.vitrac@agroparistech.fr
# Revision history: 2024-12-20


# -----------------------------------
#           INITIALIZATION
# -----------------------------------

# Ensure the script is run from the correct directory (Pizza3/utils/)
if [[ ! -f "pdocme.sh" ]]; then
    echo "Error: This script must be run from the Pizza3/utils/ directory."
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
echo "Pizza3 Version: $__version__"


# -----------------------------------
# Configuration Section
# -----------------------------------

VERSION="${__version__}min"
MANIFEST_FILE="$(realpath ../Pizza3.mini.manifest)"
RELEASE_FOLDER="$(realpath ../release)"
OUTPUT_FILE="Pizza3_v${VERSION}.zip"
INFO_FILE="release.info.txt"
RELEASE_NAME="Pizza3_v${VERSION}"


# -----------------------------------
# Validation Checks
# -----------------------------------

# Verify that the manifest file exists
if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "Error: Manifest file '$MANIFEST_FILE' not found."
    exit 1
fi

# Create the release folder if it doesn't exist
mkdir -p "$RELEASE_FOLDER"

# Full path of the output archive
output_archive="${RELEASE_FOLDER}/${OUTPUT_FILE}"

# Check if the release archive already exists to prevent overwriting
if [[ -f "$output_archive" ]]; then
    echo "Error: Release archive for version ${VERSION} already exists at:"
    echo "       $output_archive"
    exit 1
fi

# -----------------------------------
# Release Generation Process
# -----------------------------------

echo "Generating release package for Pizza3 (v${VERSION})..."

# Create a temporary directory to collect files
temp_dir=$(mktemp -d)
if [[ ! -d "$temp_dir" ]]; then
    echo "Error: Failed to create a temporary directory."
    exit 1
fi

# Ensure temporary directory is cleaned up on exit
trap "rm -rf '$temp_dir'" EXIT

# Define the release directory inside the temporary directory
release_dir="${temp_dir}/${RELEASE_NAME}"
mkdir -p "$release_dir"

echo "Using manifest file: $MANIFEST_FILE"

# Read the manifest file and copy listed files/directories into the release directory
while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip comments and empty lines
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue

    # Define the absolute path of the source file or directory
    abs_path="${mainfolder}/${line}"

    if [[ -e "$abs_path" ]]; then
        # Determine the destination directory structure
        dest_dir="${release_dir}/$(dirname "$line")"
        mkdir -p "$dest_dir"

        # Copy files and directories
        if [[ -d "$abs_path" ]]; then
            cp -r "$abs_path" "$dest_dir/"
        else
            cp -p "$abs_path" "$dest_dir/"
        fi
    else
        echo "Warning: File or directory not found in manifest: $line"
    fi
done < "$MANIFEST_FILE"

# -----------------------------------
# Adding Release Information
# -----------------------------------

# Create the release information file
info_file_path="${release_dir}/${INFO_FILE}"
{
    echo "Pizza3 Release v${VERSION}"
    echo "Generated on: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Contact: INRAE\\olivier.vitrac@agroparistech.fr"
    echo
    echo "### Included Files:"
    echo
    tree "$release_dir" --noreport | sed 's/^/    /'
} > "$info_file_path"

# -----------------------------------
# Creating the ZIP Archive
# -----------------------------------

echo "Compressing files..."

# Navigate to the temporary directory to ensure correct ZIP structure
pushd "$temp_dir" > /dev/null

# Create the ZIP archive with maximum compression
zip -rTp -9 "$output_archive" "$RELEASE_NAME" > /dev/null

# Capture the exit status of the ZIP command
zip_status=$?

# Return to the original directory
popd > /dev/null

# Check if ZIP was successful
if [[ $zip_status -eq 0 ]]; then
    echo "Release package created successfully: $output_archive"
else
    echo "Error: Failed to create the release archive."
    exit 1
fi

# -----------------------------------
# Displaying the ZIP File Size
# -----------------------------------

# Get the size of the ZIP file in a human-readable format
zip_size=$(du -h "$output_archive" | cut -f1)
echo "Release package size: $zip_size"

# -----------------------------------
# Cleanup and Completion
# -----------------------------------

# The trap set earlier will automatically remove the temporary directory
echo "Release generation complete."
