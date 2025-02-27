#!/usr/bin/env bash

###############################################################################
#                                                                             #
#    generate_manifest.sh                                                     #
#                                                                             #
#    Script to manage project manifests for Pizza3 (create, diff, update)     #
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

### **Example Usage**

# 1. **Creating a New Manifest:**
   
#    ```bash
#    ./generate_manifest.sh create -p /path/to/project -m project.manifest
#    ./generate_manifest.sh create -p /path/to/project -m Pizza3.manifest --exclude-dirs .git tmp build --exclude-files __init__.py debug.py
#    ```
   
#    **Explanation:**
#    - **`create`:** Command to generate a new manifest.
#    - **`-p /path/to/project`:** Specifies the path to the project directory.
#    - **`-m project.manifest`:** Names the manifest file as `project.manifest`.

# 2. **Diffing Two Manifests:**
   
#    ```bash
#    ./generate_manifest.sh diff -l /local/manifest/dir -s https://example.com/source/manifest --print
#    ```
   
#    **Explanation:**
#    - **`diff`:** Command to compare manifests.
#    - **`-l /local/manifest/dir`:** Specifies the local manifest directory.
#    - **`-s https://example.com/source/manifest`:** Provides the source manifest location (can be a URL or a local path).
#    - **`--print`:** Option to print the differences identified.

# 3. **Updating the Local Manifest:**
   
#    ```bash
#    ./generate_manifest.sh update -l /local/manifest/dir -s https://example.com/source/manifest --prompt
#    ```
   
#    **Explanation:**
#    - **`update`:** Command to synchronize the local manifest with the source.
#    - **`-l /local/manifest/dir`:** Specifies the local manifest directory.
#    - **`-s https://example.com/source/manifest`:** Provides the source manifest location.
#    - **`--prompt`:** Option to prompt the user before downloading and applying updates.

# Exit immediately if a command exits with a non-zero status
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 [create|diff|update] [options]"
    echo
    echo "Commands:"
    echo "  create     Create a new manifest."
    echo "  diff       Compare local and source manifests."
    echo "  update     Update local manifest based on source."
    echo
    echo "Options for create:"
    echo "  -p, --path <path>           Path to the project directory."
    echo "  -m, --manifest <manifest>   Name of the manifest file. Default: .dmanifest"
    echo "  --exclude-dirs <dirs>       List of directories to exclude (space-separated)."
    echo "  --exclude-files <files>     List of files to exclude (space-separated)."
    echo "  --extensions <exts>         List of file extensions to include (space-separated). Example: .py .sh .md .html"
    echo
    echo "Options for diff:"
    echo "  -l, --local <local_dir>     Local manifest directory."
    echo "  -s, --source <source>       Source manifest location (URL or path)."
    echo "  --local-manifest <manifest> Name of the local manifest file. Default: .dmanifest"
    echo "  --source-manifest <manifest> Name of the source manifest file. Default: .dmanifest"
    echo "  --print                      Print differences."
    echo
    echo "Options for update:"
    echo "  -l, --local <local_dir>     Local manifest directory."
    echo "  -s, --source <source>       Source manifest location (URL or path)."
    echo "  --local-manifest <manifest> Name of the local manifest file. Default: .dmanifest"
    echo "  --source-manifest <manifest> Name of the source manifest file. Default: .dmanifest"
    echo "  --prompt                     Prompt before downloading updates."
    echo
    echo "Example Commands:"
    echo "  $0 create -p /path/to/project -m Pizza3.manifest --exclude-dirs .git tmp build --exclude-files __init__.py debug.py --extensions .py .sh .md .html"
    echo "  $0 diff -l /local/manifest/dir -s https://example.com/source/manifest --print"
    echo "  $0 update -l /local/manifest/dir -s https://example.com/source/manifest --prompt"
    exit 1
}

# Check if the script is run from the utils/ directory
if [[ ! -f "generate_manifest.sh" ]]; then
    echo "Error: This script must be run from the Pizza3/utils/ directory."
    exit 1
fi

# Identify the main project directory (parent directory)
mainfolder="$(realpath ../)"

# Configuration
PYTHON_VERSION="3.10"
additional_paths=(
    "$mainfolder"
    "$mainfolder/pizza"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION/lib-dynload"
    "$HOME/anaconda3/lib/python$PYTHON_VERSION/site-packages"
    "$HOME/.ipython"
    "$(pwd)"  # Include the current directory to access generate_manifest.py
)

# Set PYTHONPATH dynamically
export PYTHONPATH=$(IFS=:; echo "${additional_paths[*]}")
echo "PYTHONPATH set to: $PYTHONPATH"

# Path to generate_manifest.py
generate_manifest_script="$(pwd)/generate_manifest.py"

# Check if generate_manifest.py exists and is executable
if [[ ! -x "$generate_manifest_script" ]]; then
    echo "Error: generate_manifest.py not found or not executable at $generate_manifest_script"
    exit 1
fi

# Parse the first argument as the command
command="$1"
shift

# Initialize variables for options
create_path=""
create_manifest=".dmanifest"
create_exclude_dirs=()
create_exclude_files=()
create_extensions=()

diff_local=""
diff_source=""
diff_local_manifest=".dmanifest"
diff_source_manifest=".dmanifest"
diff_print=false

update_local=""
update_source=""
update_local_manifest=".dmanifest"
update_source_manifest=".dmanifest"
update_prompt=false

# Function to create a new manifest
create_manifest() {
    # Prepare exclusion arguments
    exclude_dirs_args=()
    exclude_files_args=()
    extensions_args=()
    if [[ ${#create_exclude_dirs[@]} -gt 0 ]]; then
        exclude_dirs_args=(--exclude-dirs "${create_exclude_dirs[@]}")
    fi
    if [[ ${#create_exclude_files[@]} -gt 0 ]]; then
        exclude_files_args=(--exclude-files "${create_exclude_files[@]}")
    fi
    if [[ ${#create_extensions[@]} -gt 0 ]]; then
        extensions_args=(--extensions "${create_extensions[@]}")
    fi

    python3 "$generate_manifest_script" create -p "$create_path" -m "$create_manifest" "${exclude_dirs_args[@]}" "${exclude_files_args[@]}" "${extensions_args[@]}"
    if [[ $? -eq 0 ]]; then
        echo "Manifest '$create_manifest' created at '$create_path'."
    else
        echo "Error: Failed to create manifest '$create_manifest'."
    fi
}

# Function to diff manifests
diff_manifests() {
    python3 "$generate_manifest_script" diff -l "$diff_local" -s "$diff_source" --local-manifest "$diff_local_manifest" --source-manifest "$diff_source_manifest" $( [[ $diff_print == true ]] && echo "--print")
    if [[ $? -eq 0 ]]; then
        echo "Diff operation completed."
    else
        echo "Error: Diff operation failed."
    fi
}

# Function to update manifest
update_manifest() {
    python3 "$generate_manifest_script" update -l "$update_local" -s "$update_source" --local-manifest "$update_local_manifest" --source-manifest "$update_source_manifest" $( [[ $update_prompt == true ]] && echo "--prompt")
    if [[ $? -eq 0 ]]; then
        echo "Manifest update completed successfully."
    else
        echo "Error: Manifest update failed."
    fi
}

# Parse options based on command
case "$command" in
    create)
        while [[ "$#" -gt 0 ]]; do
            case "$1" in
                -p|--path)
                    create_path="$2"
                    shift 2
                    ;;
                -m|--manifest)
                    create_manifest="$2"
                    shift 2
                    ;;
                --exclude-dirs)
                    shift
                    while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
                        create_exclude_dirs+=("$1")
                        shift
                    done
                    ;;
                --exclude-files)
                    shift
                    while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
                        create_exclude_files+=("$1")
                        shift
                    done
                    ;;
                --extensions)
                    shift
                    while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
                        create_extensions+=("$1")
                        shift
                    done
                    ;;
                *)
                    echo "Unknown option for create: $1"
                    usage
                    ;;
            esac
        done
        if [[ -z "$create_path" ]]; then
            echo "Error: Path to the project directory is required for create."
            usage
        fi
        create_manifest
        ;;
    diff)
        while [[ "$#" -gt 0 ]]; do
            case "$1" in
                -l|--local)
                    diff_local="$2"
                    shift 2
                    ;;
                -s|--source)
                    diff_source="$2"
                    shift 2
                    ;;
                --local-manifest)
                    diff_local_manifest="$2"
                    shift 2
                    ;;
                --source-manifest)
                    diff_source_manifest="$2"
                    shift 2
                    ;;
                --print)
                    diff_print=true
                    shift
                    ;;
                *)
                    echo "Unknown option for diff: $1"
                    usage
                    ;;
            esac
        done
        if [[ -z "$diff_local" || -z "$diff_source" ]]; then
            echo "Error: Both local and source manifest locations are required for diff."
            usage
        fi
        diff_manifests
        ;;
    update)
        while [[ "$#" -gt 0 ]]; do
            case "$1" in
                -l|--local)
                    update_local="$2"
                    shift 2
                    ;;
                -s|--source)
                    update_source="$2"
                    shift 2
                    ;;
                --local-manifest)
                    update_local_manifest="$2"
                    shift 2
                    ;;
                --source-manifest)
                    update_source_manifest="$2"
                    shift 2
                    ;;
                --prompt)
                    update_prompt=true
                    shift
                    ;;
                *)
                    echo "Unknown option for update: $1"
                    usage
                    ;;
            esac
        done
        if [[ -z "$update_local" || -z "$update_source" ]]; then
            echo "Error: Both local and source manifest locations are required for update."
            usage
        fi
        update_manifest
        ;;
    *)
        echo "Error: Unknown command '$command'"
        usage
        ;;
esac
