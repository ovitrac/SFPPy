#!/bin/bash
:<<'END_DOC'
Script Name: refresh_alldocs.sh
Purpose: Refresh the entire documentation for the Pizza3 project, including updating setup configuration, archiving
          previous documentation, and regenerating all associated documentation files.
Description:
This script ensures the setup configuration and documentation for the Pizza3 project are refreshed. It performs the
following steps in order:

1. **Environment Check**:
   - Verifies if a Conda environment is active. 
   - Provides a warning if no environment is active, with an option for the user to proceed without one.

2. **Working Directory Validation**:
   - Ensures the script is executed from the `utils` directory of the project.

3. **Archive Existing Documentation**:
   - Archives the `html` directory, if it exists, into the `history` directory.
   - The archive is named with the format: `username@hostname_YYYY-MM-DD_HH-MM.zip`.

4. **Remove Existing Documentation**:
   - Deletes the `html` directory to allow for regeneration of the documentation.

5. **Run Documentation Generation Scripts**:
   - Executes the following scripts in the specified order:
     1. **`generate_setup.py`** (Python): Updates the setup configuration, including version information, used in the documentation.
     2. **`generate_matlab_docs.py`** (Python): Generates MATLAB-related documentation.
     3. **`generate_post_docs.py`** (Python): Generates POST-related documentation.
     4. **`generate_diagrams.sh`** (Bash): Creates project diagrams.
     5. **`pdocme.sh`** (Bash): Generates Python API documentation.
     6. **`convert_py_to_html.py`** (Python): Converts Python scripts to HTML.

6. **Open Generated Documentation**:
   - Opens relevant HTML files for review, including:
     - MATLAB documentation (`index_matlab.html`).
     - POST documentation (`index_post.html`).
     - Main documentation (`index.html`).

7. **Completion**:
   - Displays a final success message summarizing the actions performed and the location of archived documentation.

Dependencies:
- The following scripts must exist and be executable in the `utils` directory:
  - `generate_setup.py` (Python)
  - `generate_matlab_docs.py` (Python)
  - `generate_post_docs.py` (Python)
  - `generate_diagrams.sh` (Bash)
  - `pdocme.sh` (Bash)
  - `convert_py_to_html.py` (Python)
- The script assumes `zip` and `xdg-open` commands are available on the system.

Exclusions:
- The script does not attempt to regenerate excluded documentation components not listed above.
- Skips archiving if the `html` directory does not exist.

Usage:
Run the script from the `utils` directory of the Pizza3 project:
```bash
cd /path/to/Pizza3/utils
./refresh_alldocs.sh

Contact:
    Author: INRAE\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr

Last Revised: 2025-01-17
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

# Function: Check if a Conda environment is active
check_conda_env() {
    # Check if the CONDA_DEFAULT_ENV variable is set (indicates an active Conda environment)
    if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
        echo "────────────────────────────────────────────────────────────────────"
        echo "⚠️  Warning: No Conda environment is currently activated."
        echo "It is highly recommended to run this script within a Conda environment"
        echo "to ensure the proper dependencies are available and compatible."
        echo "────────────────────────────────────────────────────────────────────"
        read -p "❓ Do you want to proceed without a Conda environment? (yes/no) [no]: " user_input
        
        # Default to "no" if input is empty
        user_input=${user_input:-no}

        if [[ "$user_input" != "yes" ]]; then
            echo "────────────────────────────────────────────────────────────────────"
            echo "❌ Exiting: No Conda environment is active. Activate an environment"
            echo "before rerunning this script to ensure a smooth execution."
            echo "────────────────────────────────────────────────────────────────────"
            return 1
        else
            echo "────────────────────────────────────────────────────────────────────"
            echo "⚠️  Proceeding without a Conda environment. Be cautious of potential"
            echo "dependency issues or conflicts."
            echo "────────────────────────────────────────────────────────────────────"
        fi
    else
        echo "────────────────────────────────────────────────────────────────────"
        echo "✅ Conda environment detected: '$CONDA_DEFAULT_ENV'."
        echo "────────────────────────────────────────────────────────────────────"
    fi
    return 0
}

# Step 0: Ensure a Conda environment is active before running the script
# It is not mandatory, but running the script outside of Conda might cause issues
# with dependencies, especially with older Python versions.
if check_conda_env; then
    echo "✨ Environment check passed. Proceeding with the full script execution..."
else
    echo "❌ Environment check failed. Script execution halted." >&2
    exit 1
fi

# Step 1: Ensure the script is run from the 'utils' directory
current_dir=$(basename "$PWD")
if [[ "$current_dir" != "utils" ]]; then
    echo "Error: This script must be run from the 'utils' directory."
    exit 1
fi

# Define paths
HTML_DIR="../html"
HISTORY_DIR="../history"

# Step 2: Archive the existing 'html' directory
print_header "Step 1/3: Archiving the existing 'html' directory"

# Create history directory if it doesn't exist
if [[ ! -d "$HISTORY_DIR" ]]; then
    echo "Creating history directory at '$HISTORY_DIR'."
    mkdir -p "$HISTORY_DIR"
    echo "History directory created."
fi

# Check if 'html' directory exists
if [[ -d "$HTML_DIR" ]]; then
    # Get username and hostname
    USERNAME=$(whoami)
    HOSTNAME=$(hostname)
    
    # Get current timestamp in 'yyyy-mm-dd_HH-MM' format
    TIMESTAMP=$(date '+%Y-%m-%d_%H-%M')
    
    # Define zip filename
    ZIP_FILENAME="${USERNAME}@${HOSTNAME}_${TIMESTAMP}.zip"
    
    # Define full path for the zip file
    ZIP_PATH="${HISTORY_DIR}/${ZIP_FILENAME}"
    
    # Archive the 'html' directory with maximum compression
    echo "Archiving '$HTML_DIR' to '$ZIP_PATH' with maximum compression (-9)..."
    zip -r -9 "$ZIP_PATH" "$HTML_DIR" > /dev/null
    echo "Archiving completed successfully."
else
    echo "Warning: '$HTML_DIR' does not exist. Skipping archiving."
fi

# Step 3: Remove the existing 'html' directory and regenerate documentation
print_header "Step 2/3: Removing the existing 'html' directory and regenerating documentation"

# Remove the existing 'html' directory
if [[ -d "$HTML_DIR" ]]; then
    echo "Removing the existing '$HTML_DIR' directory..."
    rm -rf "$HTML_DIR"
    echo "Removal completed successfully."
else
    echo "No existing '$HTML_DIR' directory found. Proceeding to generate new documentation."
fi

# Check if the required scripts exist and are executable
check_script "./generate_setup.py"
check_script "./generate_all.py"
check_script "./generate_matlab_docs.py"
check_script "./generate_post_docs.py"
check_script "./generate_diagrams.sh"
check_script "./convert_py_to_html.py"
check_script "./pdocme.sh"

# Run the documentation generation scripts with status messages
echo "Launch 'generate_setup.py'to update the version (repored also in documentatation)..."
./generate_matlab_docs.py
echo "Completed 'generate_setup.py'."

echo "Launch 'generate_all.py'to refresh all variables __all__..."
./generate_all.py
echo "Completed 'generate_all.py'."

# Commented (not used in SFPPy)
# echo "Running 'generate_matlab_docs.py'..."
# ./generate_matlab_docs.py
# echo "Completed 'generate_matlab_docs.py'."
# xdg-open ../html/index_matlab.html

# Commented (not used in SFPPy)
# echo "Running 'generate_post_docs.py'..."
# ./generate_post_docs.py
# echo "Completed 'generate_matlab_docs.py'."
# xdg-open ../html/post/index_post.html

echo "Running 'generate_diagrams.sh'..."
./generate_diagrams.sh
echo "Completed 'generate_diagrams.sh'."

echo "Running 'pdocme.sh'..."
./pdocme.sh
echo "Completed 'pdocme.sh'."
xdg-open ../html/index.html

# Step 4: Final Completion Message
print_header "Step 3/3: Documentation refresh completed successfully"

echo "All documentation steps have been executed successfully."
echo "The new 'html' documentation has been generated."
echo "Previous documentation archives can be found in the '$HISTORY_DIR' directory."
