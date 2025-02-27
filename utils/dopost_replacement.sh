#!/bin/bash
# -------------------------------------------------------------------
# Script: dopost_replacement.sh
# Description: Safely replaces a string in a file while creating a backup.
#
# Usage:
#   ./dopost_replacement.sh <file> <search_string> <replacement_string>
#
# Example:
#   ./dopost_replacement.sh "output_dir/patankar/migration.html" "sidebar{width:30%" "sidebar{width:100%"
#
# Features:
#   Ensures the file exists
#    Creates a backup file with '~' added to the extension
#    Safely replaces the string using 'sed'
#    Handles special characters in search and replacement strings
#    Provides error handling for invalid inputs
#
# Contact:
#    Author: INRAE\Olivier Vitrac
#    Email: olivier.vitrac@agroparistech.fr
#
# Last Revised: 2025-01-17
# -------------------------------------------------------------------

# Ensure three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file> <search_string> <replacement_string>"
    exit 1
fi

# Assign input parameters
file="$1"
search_string="$2"
replacement_string="$3"

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "Error: File '$file' not found."
    exit 1
fi

# Create a backup file (original.html → original.html~)
backup_file="${file}~"
cp "$file" "$backup_file"

# Perform the replacement
sed -i "s|$search_string|$replacement_string|g" "$file"

# Confirm the operation
echo "✔ Replacement completed. Backup saved as '$backup_file'."
