#!/bin/bash
:<<'END_DOC'
Script Name: refresh_allreleases.sh
Purpose: Refresh all releases of the SFPPy project by archiving previous releases and generating
new full and minimalist releases for GitHub.

Description:
This script automates the process of managing SFPPy releases by performing the following steps:
1. **Execution Directory Validation**:
   - Ensures the script is executed from the `utils` directory.
2. **Archiving Old Releases**:
   - Moves existing release zip files from the `release` directory to the `release/history` directory.
   - Appends a timestamp (`YYYY-MM-DD_HH-MM`) to the filename for version control.
3. **Generating New Releases**:
   - Executes scripts to generate both full and minimalist releases:
       - Full release: Includes all features and files.
       - Minimalist release: Lightweight version with essential files.
4. **Displaying Release Summary**:
   - Summarizes the newly generated release zip files and their sizes.
   - Lists all archived releases.

Dependencies:
- **Required Scripts**:
  - `generate_simple_manifest.py` (Python): Generates the manifest for the full release.
  - `generate_release.sh` (Bash): Creates the full release zip file.
  - `generate_mini_manifest.py` (Python): Generates the manifest for the minimalist release.
  - `generate_mini_release.sh` (Bash): Creates the minimalist release zip file.
- **Commands**:
  - `find`, `du`, `mv`, `bash`, `python3`.

Execution Order:
1. `generate_simple_manifest.py` (Python)
2. `generate_release.sh` (Bash)
3. `generate_mini_manifest.py` (Python)
4. `generate_mini_release.sh` (Bash)

Exclusions:
- The script skips subdirectories when looking for zip files in the `release` directory.
- Does not generate additional release types beyond full and minimalist.

Output:
- **Archived Releases**:
  - Old releases are moved to `release/history` with a timestamp appended to their filenames.
- **New Releases**:
  - Full and minimalist releases are saved in the `release` directory.
  - The script displays a summary of all new and archived release files, including their sizes.

Usage:
- Ensure the script is run from the `utils` directory:
  ```bash
  cd /path/to/SFPPy/utils
  ./refresh_allreleases.sh


Contact:
    Author: INRAE\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr

Last Revised: 2025-01-09
END_DOC

# Exit immediately if a command exits with a non-zero status
set -e

# Function to print headers
print_header() {
    echo ""
    echo "=============================================="
    echo "$1"
    echo "=============================================="
    echo ""
}

# Function to check if a script exists and is executable
check_script() {
    if [[ ! -x "$1" ]]; then
        echo "Error: '$1' does not exist or is not executable."
        exit 1
    fi
}

# Step 1: Verify Execution Directory
print_header "Step 1/3: Verifying Execution Directory"

current_dir=$(basename "$PWD")
if [[ "$current_dir" != "utils" ]]; then
    echo "Error: This script must be run from the 'utils' directory."
    exit 1
else
    echo "Verified: Running from 'utils' directory."
fi

# Define absolute paths
mainfolder="$(realpath ../)"
releases_dir="${mainfolder}/release"
history_dir="${releases_dir}/history"

# Step 2: Archive Existing Releases
print_header "Step 2/3: Archiving Existing Releases"

# Create history directory if it doesn't exist
if [[ ! -d "$history_dir" ]]; then
    echo "Creating history directory at '$history_dir'."
    mkdir -p "$history_dir"
    echo "History directory created."
fi

# Find zip files in releases_dir (excluding subdirectories)
mapfile -t zip_files < <(find "$releases_dir" -maxdepth 1 -type f -name "*.zip")

if [[ ${#zip_files[@]} -eq 0 ]]; then
    echo "No release zip files found in '$releases_dir'. Skipping archiving."
else
    echo "Archiving the following release zip files:"
    printf "%-50s %10s\n" "Filename" "Size"
    printf "%-50s %10s\n" "--------" "----"
    
    for zip_file in "${zip_files[@]}"; do
        # Get filename without path
        filename=$(basename "$zip_file")
        
        # Get the modification date in YYYY-MM-DD_HH-MM format
        mod_date=$(date -r "$zip_file" '+%Y-%m-%d_%H-%M')
        
        # Define new filename with date
        new_filename="${filename%.zip}_$mod_date.zip"
        
        # Move and rename the zip file to history_dir
        mv "$zip_file" "${history_dir}/${new_filename}"
        
        # Get the size of the moved file in human-readable format
        filesize=$(du -h "${history_dir}/${new_filename}" | cut -f1)
        
        # Print the details
        printf "%-50s %10s\n" "$new_filename" "$filesize"
    done
fi

echo ""

# Step 3: Generating New Releases
print_header "Step 3/3: Generating New Releases"

# Define paths to generation scripts
generate_simple_manifest="${PWD}/generate_simple_manifest.py"
generate_release="${PWD}/generate_release.sh"
generate_mini_manifest="${PWD}/generate_mini_manifest.py"
generate_mini_release="${PWD}/generate_mini_release.sh"

# Check if the required scripts exist and are executable
check_script "$generate_simple_manifest"
check_script "$generate_release"
check_script "$generate_mini_manifest"
check_script "$generate_mini_release"

# Function to generate a release
generate_release_type() {
    local manifest_script="$1"
    local release_script="$2"
    local release_type="$3"
    
    print_header "Generating ${release_type} Release"
    
    echo "Running '${manifest_script##*/}'..."
    python3 "$manifest_script"
    echo "Manifest generated successfully."
    
    echo "Running '${release_script##*/}'..."
    bash "$release_script"
    echo "Completed '${release_script##*/}'."
}

# Generate Full Release
generate_release_type "$generate_simple_manifest" "$generate_release" "Full"

# Generate Minimalist Release
#generate_release_type "$generate_mini_manifest" "$generate_mini_release" "Minimalist"

echo ""

# Step 4: Display Details of New Releases
print_header "Release Generation Summary"

# Find new zip files in releases_dir (excluding history/)
mapfile -t new_zip_files < <(find "$releases_dir" -maxdepth 1 -type f -name "*.zip")

if [[ ${#new_zip_files[@]} -eq 0 ]]; then
    echo "No new release zip files were generated."
else
    echo "Newly generated release zip files:"
    printf "%-50s %10s\n" "Filename" "Size"
    printf "%-50s %10s\n" "--------" "----"
    
    for new_zip in "${new_zip_files[@]}"; do
        filename=$(basename "$new_zip")
        filesize=$(du -h "$new_zip" | cut -f1)
        printf "%-50s %10s\n" "$filename" "$filesize"
    done
fi

echo ""
echo "All release generation steps have been executed successfully."
echo "Archived releases can be found in the '$history_dir' directory."
