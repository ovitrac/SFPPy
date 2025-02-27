#!/usr/bin/env bash

# Backup Script for SFPPy Project
# INRAE\Olivier Vitrac - olivier.vitrac@agroparistech.fr
# Revision: 2024-01-09

# ----------------------------------------------------------------------
# Script Overview:
# `backupme.sh` is designed to back up code files and other related files
# in the SFPPy project directory. It identifies and archives files based
# on defined inclusion patterns while excluding specific folders and files.
# The output is a compressed `.zip` file containing the backed-up files,
# along with detailed Markdown and HTML reports.
#
# Reports:
# - Markdown Report: `history/backupme.README.md`
# - HTML Report: `history/backupme.README.html` (embedded CSS)
#
# ----------------------------------------------------------------------

# Exit immediately if a command exits with a non-zero status
set -e
# Treat unset variables as an error
set -u

# ----------------------------------------------------------------------
# Usage:
# 1. Run the script from the `SFPPy/utils/` directory:
#    ./backupme.sh [ -y ]
#    -y: Optional flag to force yes without confirmation.
# 2. The backup archive will be created in the `history` folder under
#    the SFPPy project directory.
# 3. After execution, the script will display:
#    - Confirmation of the backup.
#    - The full path to the backup file.
#    - The file size in a human-readable format (e.g., B, kB, MB).
#    - Detailed backup reports included in `backupme.README.md` and `backupme.README.html` inside the zip.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Features:
# 1. Automatic Directory Validation:
#    Ensures the script is executed from the correct location, avoiding unintended behavior.
# 2. Dynamic File Inclusion:
#    Includes files matching patterns defined in the `include_patterns` array.
# 3. Folder and File Exclusions:
#    Skips folders and files specified in the `exclude_folders_rel`, `exclude_folders_abs`, and `exclude_files` arrays.
# 4. Human-Readable File Size:
#    Displays the size of the generated backup file in a user-friendly format.
# 5. Custom Backup Name:
#    The backup file name includes the current directory name, username, hostname, and timestamp for easy identification.
# 6. Backup Reports:
#    Generates both `backupme.README.md` and `backupme.README.html` markdown files with backup details, included in the zip.
# 7. Progress Indicators and Timing:
#    Provides textual indicators and timing for different script steps.
# 8. Pandoc Installation Check:
#    Verifies whether Pandoc is installed before attempting to generate the HTML report.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Configuration:
# Modify the following arrays to customize the script:
# - `include_patterns`: Add or remove file extensions to include in the backup.
# - `exclude_folders_rel`: Add or remove relative folders to exclude.
# - `exclude_folders_abs`: Add or remove absolute folders to exclude.
# - `exclude_files`: Add or remove specific file patterns to exclude.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Function Definitions
# ----------------------------------------------------------------------

# Function to display usage information
usage() {
    echo "Usage: $0 [ -y ]"
    echo "  -y    Force yes to proceed with the backup without confirmation."
    exit 1
}

# Function to check if Pandoc is installed and version is 3.x or later
check_pandoc_installed() {
    # Check if pandoc command exists
    if ! command -v pandoc >/dev/null 2>&1; then
        echo "Error: Pandoc is not installed. Please install Pandoc to generate HTML reports."
        echo "Visit https://pandoc.org/installing.html for installation instructions."
        return 1
    fi
    # Check pandoc version
    local pandoc_version
    pandoc_version=$(pandoc --version 2>/dev/null | head -n 1 | awk '{print $2}')
    if [[ -z "$pandoc_version" ]]; then
        echo "Error: Unable to determine Pandoc version. Ensure Pandoc is correctly installed."
        return 1
    fi
    # Parse the major version
    local major_version
    major_version=$(echo "$pandoc_version" | cut -d. -f1)
    # Check if the major version is 3 or later
    if [[ "$major_version" -lt 3 ]]; then
        echo "Error: Pandoc version is $pandoc_version (earlier than 3.x). Please upgrade to Pandoc 3.x or later."
        echo "Visit https://pandoc.org/installing.html for upgrade instructions."
        return 1
    fi
    echo "Pandoc is installed and version is $pandoc_version (3.x or later)."
    return 0
}

# Function to convert size to best unit
convert_size() {
    local size=$1
    if (( size < 1024 )); then
        echo "${size} B"
    elif (( size < 1024 * 1024 )); then
        printf "%.2f kB" "$(echo "scale=2; $size / 1024" | bc)"
    else
        printf "%.2f MB" "$(echo "scale=2; $size / (1024 * 1024)" | bc)"
    fi
}

# Function to calculate elapsed time with two decimal precision
calculate_elapsed_time() {
    local start=$1
    local end=$2
    elapsed=$(echo "$end - $start" | bc)
    # Format to two decimal places
    printf "%.2f" "$elapsed"
}

# Function to truncate strings exceeding a certain length
truncate_string() {
    local str=$1
    local max_length=$2
    if [[ ${#str} -gt $max_length ]]; then
        echo "${str:0:$((max_length-3))}..."
    else
        echo "$str"
    fi
}

# ----------------------------------------------------------------------
# Parse Optional Parameters
# ----------------------------------------------------------------------

force_yes=false
while getopts "y" opt; do
    case $opt in
        y)
            force_yes=true
            ;;
        *)
            usage
            ;;
    esac
done

# ----------------------------------------------------------------------
# Ensure Script is Run from the Correct Directory
# ----------------------------------------------------------------------

# This script must be run from SFPPy/utils/
if [[ ! -f "pdocme.sh" ]]; then
    echo "Error: This script must be run from the SFPPy/utils/ directory."
    exit 1
fi

# ----------------------------------------------------------------------
# Define Paths
# ----------------------------------------------------------------------

mainfolder="$(realpath ../)" # Main project directory
backfolder="$mainfolder/history" # Backup folder
mkdir -p "$backfolder" # Create backup folder if it doesn't exist

# Temporary files for storing find results
tempfile=$(mktemp)
tempfile_with_sizes=$(mktemp)
tempfile_zip=$(mktemp)

# Other definitions
CSSfile="$mainfolder/utils/css/sfppy.css"

# ----------------------------------------------------------------------
# Define Inclusion and Exclusion Patterns
# ----------------------------------------------------------------------

# Patterns for including files in the backup
include_patterns=(
    "*.m"         # MATLAB files
    "*.asv"       # Auto-saved files
    "*.m~"        # Backup files
    "*.pynb"      # Jupyter notebooks
    "*.py"        # Python scripts
    "*.sh"        # Shell scripts
    "*.txt"       # Text files
    "*.md"        # Markdown files
    "*.html"      # HTML files
    "*.json"      # JSON files
    "*.css"       # CSS files
    "*.manifest"  # Manifest files
)

# Relative folders to exclude from the backup
exclude_folders_rel=(
    "./old"       # Folder for old files (local)
    "./tmp"       # Temporary folder (local)
    "./sandbox"   # Sandbox folder (local)
    "./debug"     # Debug folder (local)
    "./obsolete"  # Obsolete folder (local)
    "./.git"      # Git folder (local)
    "./.vscode"   # Visual Studio Code folder (local)
    "./.spyproject" # Spyder project folder (local)
    "./__all__"   # all (internal to SFPPy)
    "./__pycache__" # Python folder (local)
)

# Absolute folders to exclude from the backup
exclude_folders_abs=(
    "$mainfolder/history"  # History folder
    "$mainfolder/release"  # Release folder
)

# Specific file patterns to exclude from the backup
exclude_files=(
    "*.log"                # Log files
    "*.zip"                # ZIP files
    "backupme.README.md"   # Exclude README.md from backup
)

# ----------------------------------------------------------------------
# Start Total Timer
# ----------------------------------------------------------------------

total_start_time=$(date +%s.%N)

# ----------------------------------------------------------------------
# Step 1: Find Files to Include in the Backup
# ----------------------------------------------------------------------

echo "Step 1: Searching for files to include in the backup..."
step1_start_time=$(date +%s.%N)

# Construct the find command
find_cmd="find \"$mainfolder\" \( -type f \( "
for pattern in "${include_patterns[@]}"; do
    find_cmd+=" -iname \"$pattern\" -o"
done
find_cmd="${find_cmd% -o}" # Remove trailing -o
find_cmd+=" \) \)"

# Exclude absolute folders
for folder in "${exclude_folders_abs[@]}"; do
    find_cmd+=" ! -path \"$folder/*\""
done

# Exclude relative folders (converted to glob patterns)
for folder in "${exclude_folders_rel[@]}"; do
    folder_basename=$(basename "$folder")
    find_cmd+=" ! -path \"*/$folder_basename/*\""
done

# Exclude specific files
for file in "${exclude_files[@]}"; do
    find_cmd+=" ! -iname \"$file\""
done

# Execute find command and save results to tempfile
eval "$find_cmd" > "$tempfile"

# Preprocess find output to include file sizes and modification dates
echo "Processing file details..."
while IFS= read -r file; do
    size=$(stat --printf="%s" "$file" 2>/dev/null || echo 0) # Get file size or 0 if stat fails
    mod_date=$(stat --printf="%y" "$file" 2>/dev/null || echo "1970-01-01 00:00:00")
    mod_date=$(date -d "$mod_date" +"%Y-%m-%d %H:%M:%S") # Limit to seconds
    relative_path=$(realpath --relative-to="$mainfolder" "$file")
    printf "%s\t%s\t%s\n" "$relative_path" "$size" "$mod_date" >> "$tempfile_with_sizes"
done < "$tempfile"

# Calculate time spent on Step 1
step1_end_time=$(date +%s.%N)
step1_elapsed=$(calculate_elapsed_time "$step1_start_time" "$step1_end_time")
echo "Step 1 completed in $step1_elapsed seconds."

# ----------------------------------------------------------------------
# Step 2: Generate README.md Report
# ----------------------------------------------------------------------

echo "Step 2: Generating backupme.README.md..."
step2_start_time=$(date +%s.%N)

# Get user@host and current date
user_host="$(whoami)@$(hostname)"
backup_date="$(date +"%Y-%m-%d %H:%M:%S")"

# Calculate total number of files and total size
total_files=$(wc -l < "$tempfile")
total_size_bytes=$(awk -F'\t' '{sum += $2} END {print sum}' "$tempfile_with_sizes")
total_size_formatted=$(convert_size "$total_size_bytes")

# Generate Rules Applied table in README.md
readme_file="$backfolder/backupme.README.md"

# Initialize README.md with header and Rules Applied table
{
    echo "# Backup Report"
    echo ""
    echo "**User:** $user_host  "
    echo "**Backup Path:** $mainfolder  "
    echo "**Backup Date:** $backup_date  "
    echo ""
    echo "## Rules Applied"
    echo "| Inclusion Rules   | Exclusion Rules                         |"
    echo "|-------------------|-----------------------------------------|"
} > "$readme_file"

# Prepare inclusion and exclusion lists
inclusion_list=("${include_patterns[@]}")
exclusion_list=("${exclude_folders_rel[@]}" "${exclude_folders_abs[@]}" "${exclude_files[@]}")

# Determine the maximum number of rows needed
max_rules=${#inclusion_list[@]}
if (( ${#exclusion_list[@]} > max_rules )); then
    max_rules=${#exclusion_list[@]}
fi

# Iterate and add rules to the table
for ((i=0; i<max_rules; i++)); do
    include="${inclusion_list[i]:-}"
    exclude="${exclusion_list[i]:-}"
    
    # Truncate exclusion if necessary for folder paths
    if [[ -n "$exclude" ]]; then
        exclude_display=$(truncate_string "$exclude" 40)
    else
        exclude_display=""
    fi

    if [[ -n "$include" ]]; then
        # Write to README.md with backticks for Markdown formatting
        printf "| __BT__%-15s__BT__ | __BT__%-37s__BT__ |\n" "$include" "$exclude_display" >> "$readme_file"
        # Write to screen without backticks and with fixed width
        printf "| %-17s | %-40s |\n" "$include" "$exclude_display"
    else
        # Handle cases with empty inclusions
        printf "| %-17s | __BT__%-37s__BT__ |\n" "" "$exclude_display" >> "$readme_file"
        # Write to screen without backticks and with fixed width
        printf "| %-17s | %-40s |\n" "" "$exclude_display"
    fi
done

# Replace placeholders with backticks in the README.md file
sed -i 's/__BT__/`/g' "$readme_file"

# Summarize by extension and write to README.md and display on screen
{
    echo ""
    echo "### Summary of Files by Extension (Total Files: $total_files, Total Size: $total_size_formatted)"
    echo "| Extension      | Count      |  Total Size  |"
    echo "|----------------|------------|--------------|"
} >> "$readme_file"

extension_data=$(awk -F'\t' '
{
    file = $1;
    size = $2;
    n = split(file, parts, "/");
    file_name = parts[n];
    split(file_name, name_parts, ".");
    if(length(name_parts) > 1){
        ext = name_parts[length(name_parts)];
    }
    else{
        ext = "unknown";
    }
    sizes_ext[ext] += size;
    counts_ext[ext]++;
}
END {
    for(ext in sizes_ext){
        total_size_kb = sizes_ext[ext] / 1024;
        if(total_size_kb < 1){
            size = sprintf("%.0f B", sizes_ext[ext]);
        }
        else if(total_size_kb < 1024){
            size_kb = total_size_kb
            size = sprintf("%.2f kB", size_kb)
        }
        else{
            size_mb = total_size_kb / 1024
            size = sprintf("%.2f MB", size_mb)
        }
        printf "| __BT__%-12s__BT__ | %-10d | %-12s |\n", ext, counts_ext[ext], size;
    }
}' "$tempfile_with_sizes" | sort -k2 -nr)

# Write extension summary to README.md with backticks
echo "$extension_data" >> "$readme_file"

# Display extension summary on screen without backticks
{
    echo ""
    echo "### Summary of Files by Extension (Total Files: $total_files, Total Size: $total_size_formatted)"
    echo "|  Extension   |    Count   |  Total Size  |"
    echo "|--------------|------------|--------------|"
}
echo "$extension_data" | sed 's/__BT__//g' | sort -k2 -nr

# Summarize by folder and write to README.md and display on screen
{
    echo ""
    echo "### Summary of Files by Folder (Total Files: $total_files, Total Size: $total_size_formatted)"
    echo "|                                      Folder                                      |    Count   |  Total Size  |"
    echo "|------------------------------------------------------------------------------------|------------|--------------|"
} >> "$readme_file"

folder_data=$(awk -F'\t' '
{
    file = $1;
    size = $2;
    n = split(file, parts, "/");
    if(n > 1){
        folder_path = substr($1, 1, length($1) - length(parts[n]) -1);
    }
    else{
        folder_path = ".";
    }
    sizes_folder[folder_path] += size;
    counts_folder[folder_path]++;
}
END {
    for(folder in sizes_folder){
        total_size_kb = sizes_folder[folder] / 1024
        if(total_size_kb < 1){
            size = sprintf("%.0f B", sizes_folder[folder])
        }
        else if(total_size_kb < 1024){
            size_kb = total_size_kb
            size = sprintf("%.2f kB", size_kb)
        }
        else{
            size_mb = total_size_kb / 1024
            size = sprintf("%.2f MB", size_mb)
        }
        printf "| __BT__%-80s__BT__ | %-10d | %-12s |\n", folder, counts_folder[folder], size;
    }
}' "$tempfile_with_sizes" | sort -k2 -nr)

# Write folder summary to README.md with backticks
echo "$folder_data" >> "$readme_file"

# Display folder summary on screen without backticks
{
    echo ""
    echo "### Summary of Files by Folder (Total Files: $total_files, Total Size: $total_size_formatted)"
    echo "|                                      Folder                                      |    Count   |  Total Size  |"
    echo "|------------------------------------------------------------------------------------|------------|--------------|"
}
echo "$folder_data" | sed 's/__BT__//g' | sort -k2 -nr

# Calculate time spent on Step 2
step2_end_time=$(date +%s.%N)
step2_elapsed=$(calculate_elapsed_time "$step2_start_time" "$step2_end_time")
echo "Step 2 completed in $step2_elapsed seconds."

# ----------------------------------------------------------------------
# Step 3: Generate Detailed List of Files by Folder
# ----------------------------------------------------------------------

echo "Step 3: Generating detailed list of files by folder..."
step3_start_time=$(date +%s.%N)

# Add Detailed List section to README.md
{
    echo ""
    echo "## Detailed List of Files by Folder" >> "$readme_file"
} 

# Get unique list of folders, sorted
folders=$(awk -F'\t' '{
    n = split($1, parts, "/");
    if(n > 1){
        folder_path = substr($1, 1, length($1) - length(parts[n]) -1);
    }
    else{
        folder_path = ".";
    }
    folders[folder_path]++;
}
END {
    for (f in folders) print f;
}' "$tempfile_with_sizes" | sort)

# Iterate over each folder to generate detailed lists
for folder in $folders; do
    # Calculate number of files and total size for the folder
    folder_info=$(awk -F'\t' -v folder="$folder" '
    {
        file = $1;
        size = $2;
        n = split(file, parts, "/");
        if(n > 1){
            current_folder_path = substr($1, 1, length($1) - length(parts[n]) -1);
        }
        else{
            current_folder_path = ".";
        }
        if(current_folder_path == folder){
            size_b += size;
            count++;
        }
    }
    END {
        if(count > 0){
            if(size_b < 1024){
                size = sprintf("%.0f B", size_b)
            }
            else if(size_b < 1024 * 1024){
                size_kb = size_b / 1024
                size = sprintf("%.2f kB", size_kb)
            }
            else{
                size_mb = size_b / (1024 * 1024)
                size = sprintf("%.2f MB", size_mb)
            }
            printf " (%d files, %s)", count, size
        }
        else{
            printf " (0 files, 0 B)"
        }
    }' "$tempfile_with_sizes")
    
    # Truncate folder name if it exceeds 80 characters
    truncated_folder=$(truncate_string "$folder" 77) # 80 - 3 for "..."
    
    # Write folder title with counts and size to README.md
    echo "" >> "$readme_file"
    if [[ "$folder" == "." ]]; then
        echo "### Folder: root$folder_info" >> "$readme_file"
    else
        echo "### Folder: $truncated_folder$folder_info" >> "$readme_file"
    fi
    echo "| File                                     | Size        | Last Modified       |" >> "$readme_file"
    echo "|------------------------------------------|-------------|---------------------|" >> "$readme_file"

    # Collect files for the folder and format as Markdown links
    folder_files=$(awk -F'\t' -v folder="$folder" '
    BEGIN {
        if(folder == "."){
            folder_prefix = "";
        }
        else{
            folder_prefix = folder "/";
        }
    }
    {
        file = $1;
        size = $2;
        mod_date = $3 " " $4;
        if(folder == "." && index(file, "/") == 0){
            file_name = file;
            size_b = size;
        }
        else if(folder != "." && index(file, folder_prefix) == 1){
            relative_file = substr(file, length(folder_prefix)+1);
            n = split(relative_file, parts, "/");
            file_name = parts[n];
            size_b = size;
        }
        else{
            next;
        }
        # Choose best unit
        if(size_b < 1024){
            size_formatted = sprintf("%-10s", size_b " B")
        }
        else if(size_b < 1024 * 1024){
            size_kb = size_b / 1024
            size_formatted = sprintf("%.2f kB", size_kb)
        }
        else{
            size_mb = size_b / (1024 * 1024)
            size_formatted = sprintf("%.2f MB", size_mb)
        }
        # Truncate filename to 30 characters, append "..." if truncated
        if(length(file_name) > 30){
            display_name = substr(file_name, 1, 27) "..."
        }
        else{
            display_name = file_name
        }
        # Create the relative link path from history/
        link_target = "../" file
        # Enclose filename in Markdown link syntax
        printf "| [%s](%s) | %-10s | %-19s |\n", display_name, link_target, size_formatted, mod_date
    }' "$tempfile_with_sizes" | sort)

    # Write detailed file list to README.md
    echo "$folder_files" >> "$readme_file"
done

# Replace any remaining placeholders with backticks in the README.md file
sed -i 's/__BT__/`/g' "$readme_file"

# Calculate time spent on Step 3
step3_end_time=$(date +%s.%N)
step3_elapsed=$(calculate_elapsed_time "$step3_start_time" "$step3_end_time")
echo "Step 3 completed in $step3_elapsed seconds."

# ----------------------------------------------------------------------
# Step 3a: Convert README.md to HTML using Pandoc
# ----------------------------------------------------------------------

echo "Step 3a: Converting README.md to HTML..."
step3a_start_time=$(date +%s.%N)

# Define output HTML file path
html_file="$backfolder/backupme.README.html"

# Check if Pandoc is installed
if check_pandoc_installed; then
    # Check if the CSS file exists
    if [ ! -f "$CSSfile" ]; then
        echo "Error: CSS file '$CSSfile' does not exist. Please ensure the file path is correct."
        exit 1
    fi
    # Execute Pandoc conversion with embedded resources and metadata
    if pandoc "$readme_file" \
              --embed-resources \
              --standalone \
              --metadata title="SFPPy" \
              -c "$CSSfile" \
              -o "$html_file"; then
        echo "HTML report generated successfully at $html_file"
    else
        echo "Error: Failed to generate HTML report with Pandoc."
        exit 1
    fi
else
    echo "Skipping HTML report generation."
    # Optionally, you can choose to exit the script instead
    # exit 1
fi

# Calculate time spent on Step 3a
step3a_end_time=$(date +%s.%N)
step3a_elapsed=$(calculate_elapsed_time "$step3a_start_time" "$step3a_end_time")
echo "Step 3a completed in $step3a_elapsed seconds."

# ----------------------------------------------------------------------
# Step 4: Prompt User for Confirmation Before Zipping
# ----------------------------------------------------------------------

echo "Step 4: Preparing to create the backup ZIP file..."
step4_start_time=$(date +%s.%N)

if ! $force_yes; then
    echo ""
    read -p "Do you want to proceed with the backup? (y/n): " response
    if [[ "$response" != "y" ]]; then
        echo "Backup canceled."
        rm "$tempfile" "$tempfile_with_sizes" "$tempfile_zip"
        exit 0
    fi
fi

# Calculate time spent on Step 4
step4_end_time=$(date +%s.%N)
step4_elapsed=$(calculate_elapsed_time "$step4_start_time" "$step4_end_time")
echo "Step 4 completed in $step4_elapsed seconds."

# ----------------------------------------------------------------------
# Step 5: Prepare List of Files to Backup, Including Reports
# ----------------------------------------------------------------------

echo "Step 5: Preparing the list of files to include in the backup ZIP..."
step5_start_time=$(date +%s.%N)

{
    # List files relative to mainfolder
    awk -F'\t' '{print $1}' "$tempfile_with_sizes"
    # Include README.md with relative path
    echo "history/backupme.README.md"
    # Include README.html with relative path if it was generated
    if [[ -f "$html_file" ]]; then
        echo "history/backupme.README.html"
    fi
} > "$tempfile_zip"

# Optional: Verify the contents of tempfile_zip (Uncomment for debugging)
# echo "Files to be zipped:"
# cat "$tempfile_zip"
# echo ""

# Calculate time spent on Step 5
step5_end_time=$(date +%s.%N)
step5_elapsed=$(calculate_elapsed_time "$step5_start_time" "$step5_end_time")
echo "Step 5 completed in $step5_elapsed seconds."

# ----------------------------------------------------------------------
# Step 6: Create the ZIP File with Paths Relative to mainfolder
# ----------------------------------------------------------------------

echo "Step 6: Creating the backup ZIP file..."
step6_start_time=$(date +%s.%N)

# Define backup filename before changing directory
backup_filename="${PWD##*/}_backup_$(whoami)@$(hostname)_$(date +"%Y_%m_%d__%H-%M").zip"
backup_path="$backfolder/$backup_filename"

# Create the zip file with paths relative to mainfolder
# Change directory to mainfolder to ensure relative paths in zip
(
    cd "$mainfolder" || exit 1
    zip -@ "$backup_path" < "$tempfile_zip"
)

# Remove temporary zip list
rm "$tempfile_zip"

# Verify the ZIP file by testing it
if zip -T "$backup_path" >/dev/null 2>&1; then
    echo "Test of $backup_path OK"
else
    echo "Error: ZIP file validation failed."
    exit 1
fi

# Calculate time spent on Step 6
step6_end_time=$(date +%s.%N)
step6_elapsed=$(calculate_elapsed_time "$step6_start_time" "$step6_end_time")
echo "Step 6 completed in $step6_elapsed seconds."

# ----------------------------------------------------------------------
# Step 7: Display Backup File Details
# ----------------------------------------------------------------------

echo "Step 7: Verifying and displaying backup details..."
step7_start_time=$(date +%s.%N)

if [[ -f "$backup_path" ]]; then
    backup_size=$(du -h "$backup_path" | cut -f1)
    echo -e "\nBackup successfully created!"
    echo "Location: $backup_path"
    echo "Size: $backup_size"
else
    echo "Error: Backup file was not created."
    exit 1
fi

# Calculate time spent on Step 7
step7_end_time=$(date +%s.%N)
step7_elapsed=$(calculate_elapsed_time "$step7_start_time" "$step7_end_time")
echo "Step 7 completed in $step7_elapsed seconds."

# ----------------------------------------------------------------------
# Step 8: Clean Up Temporary Files
# ----------------------------------------------------------------------

echo "Step 8: Cleaning up temporary files..."
step8_start_time=$(date +%s.%N)

rm "$tempfile" "$tempfile_with_sizes"

# Calculate time spent on Step 8
step8_end_time=$(date +%s.%N)
step8_elapsed=$(calculate_elapsed_time "$step8_start_time" "$step8_end_time")
echo "Step 8 completed in $step8_elapsed seconds."

# ----------------------------------------------------------------------
# Total Time Calculation
# ----------------------------------------------------------------------

total_end_time=$(date +%s.%N)
total_elapsed=$(calculate_elapsed_time "$total_start_time" "$total_end_time")
echo -e "\nTotal Time Spent: $total_elapsed seconds."

exit 0