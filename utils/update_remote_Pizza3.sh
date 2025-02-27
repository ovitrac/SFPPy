#!/usr/bin/env bash
:<<'END_DOC'
update_remote_Pizza3.sh
----------------------------------------------------------------------
This script wraps around compare_and_list.py to actually perform the
creation or replacement of files in the remote folder. It replicates
the same usage and arguments, plus supports a "-y" flag for "yes to all".

Usage:
  ./update_remote_Pizza3.sh <destination_folder> [comparison_mode] [operation_mode] [create_missing_folders_flag] [-y]

Example:
  ./update_remote_Pizza3.sh /path/to/remote both update yes

Production example:
    ./update_remote_Pizza3.sh $HOME/onedriveOV_AgroParisTech/Han/dev/Pizza3

Workflow:
  1) Calls compare_and_list.py to list missing and/or obsolete files.
  2) Reads each line from that output.
  3) For each line, determines if it's "MISSING" or "OBSOLETE".
  4) If "MISSING", optionally confirm creation + copying, or auto-confirm if -y.
  5) If "OBSOLETE", display file metadata on source & destination, ask for justification,
     then proceed if user confirms (or skip if refused).

----------------------------------------------------------------------

### Usage

```bash
cd /path/to/Pizza3/utils
./update_remote_Pizza3.sh <destination_folder> [comparison_mode] [operation_mode] [create_missing_folders_flag] [-y]
```

- **destination_folder**: path to remote (must contain `utils/pdocme.sh`)  
- **comparison_mode** (optional): `date`, `size`, or `both` (default: `date`).  
- **operation_mode** (optional): `missing` or `update` (default: `update`).  
- **create_missing_folders_flag** (optional): `yes` or `no` (default: `no`).  
- **-y** (optional): auto-confirm all missing/obsolete file actions (skips prompt for justification, logs "Auto-accepted with -y").  

### Examples:

1. **List only missing files** (do not overwrite anything), interactively deciding whether to copy them or not, *no folder creation* by default:

   ```bash
   ./update_remote_Pizza3.sh /path/to/remote missing
   ```

2. **List missing & obsolete** (default is `update`, with date-based comparison), and create needed folders automatically, *without asking for confirmation*:

   ```bash
   ./update_remote_Pizza3.sh /path/to/remote both update yes -y
   ```

   This runs:
   - `compare_and_list.py /path/to/remote both update yes`  
   - Then automatically copies all missing/obsolete files, providing the justification “Auto-accepted with -y”.

----------------------------------------------------------------------
Contact:
    Author: INRAE\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr

Last Revised: 2025-01-09
END_DOC

ALL_ARGS=("$@")
LAST_ARG="${ALL_ARGS[${#ALL_ARGS[@]}-1]}"

YES_TO_ALL="no"
if [ "$LAST_ARG" == "-y" ]; then
  YES_TO_ALL="yes"
  # Remove -y from argument list
  unset 'ALL_ARGS[${#ALL_ARGS[@]}-1]'
fi

# Rebuild argument list (for the python script) without the last -y if present
PY_ARGS=("${ALL_ARGS[@]}")

# Check we are in 'utils' folder
CURRENT_DIR=$(basename "$(pwd)")
if [ "$CURRENT_DIR" != "utils" ]; then
  echo "Error: You must run this script from the 'utils' folder." > /dev/tty
  exit 1
fi

if [ ! -f pdocme.sh ]; then
  echo "Error: 'pdocme.sh' not found in the current folder. Execution aborted." > /dev/tty
  exit 1
fi

# 1) Run compare_and_list.py, store output in a temp file
TMPFILE=$(mktemp /tmp/compare_and_list_output.XXXXXX)
./compare_and_list.py "${PY_ARGS[@]}" > "$TMPFILE"

# 2) Read the temp file line by line
while IFS= read -r line; do

  # We pass normal lines back to standard output if they don't indicate MISSING or OBSOLETE.
  # But we direct interactive stuff to /dev/tty so the user sees it in real time.
  
  # If this line is the final "Comparison completed.", just echo it to standard output
  if [[ "$line" == *"Comparison completed."* ]]; then
    echo "$line"
    continue
  fi

  # If line is neither "MISSING:" nor "OBSOLETE:", just echo it normally.
  if [[ "$line" != *"MISSING:"* && "$line" != *"OBSOLETE:"* ]]; then
    echo "$line"
    continue
  fi

  # Determine whether the action is MISSING or OBSOLETE
  if [[ "$line" == *"MISSING:"* ]]; then
    ACTION_TYPE="MISSING"
  else
    ACTION_TYPE="OBSOLETE"
  fi

  # Typically lines look like:
  #   (source: 2025-01-09, 12 KB) MISSING: /source/path -> /dest/path
  #   (+2.3 KB, +45.0 min) OBSOLETE: /source/path -> /dest/path
  # Let's parse out /source/path and /dest/path
  TMP_STR=$(echo "$line" | sed -E 's/.* (MISSING|OBSOLETE): //')
  SOURCE_PATH=$(echo "$TMP_STR" | awk -F ' -> ' '{print $1}')
  DEST_PATH=$(echo "$TMP_STR"   | awk -F ' -> ' '{print $2}')

  if [ "$YES_TO_ALL" == "no" ]; then
    # Interactive mode
    echo "" > /dev/tty
    echo "--------------------------------------------------------------------" > /dev/tty
    echo "Detected $ACTION_TYPE action:" > /dev/tty
    echo "$line" > /dev/tty

    # Show details to help user decide
    if [ "$ACTION_TYPE" == "OBSOLETE" ]; then
      echo "" > /dev/tty
      echo "Source file details:" > /dev/tty
      ls -l "$SOURCE_PATH" 2>/dev/null > /dev/tty
      echo "" > /dev/tty
      echo "Destination file details:" > /dev/tty
      ls -l "$DEST_PATH" 2>/dev/null > /dev/tty || echo "(destination missing?)" > /dev/tty
    else
      echo "" > /dev/tty
      echo "Source file details:" > /dev/tty
      ls -l "$SOURCE_PATH" 2>/dev/null > /dev/tty
      echo "" > /dev/tty
    fi

    # Prompt for yes/no
    echo -n "Do you want to proceed with this $ACTION_TYPE file? [y/N] " > /dev/tty
    read -r CONFIRM < /dev/tty
    if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
      echo "Skipping..." > /dev/tty
      continue
    else
        JUSTIFICATION="validated by user"
    fi

    # Prompt for justification
    # echo "Please provide a justification for this update:" > /dev/tty
    # read -r JUSTIFICATION < /dev/tty
    echo "You said: '$JUSTIFICATION'" > /dev/tty
    echo "" > /dev/tty

  else
    # -y => auto-confirm everything
    JUSTIFICATION="Auto-accepted with -y"
  fi

  # 3) Perform the copy/creation
  echo "Performing $ACTION_TYPE action for: $SOURCE_PATH" > /dev/tty
  DEST_DIR=$(dirname "$DEST_PATH")
  if [ ! -d "$DEST_DIR" ]; then
    mkdir -p "$DEST_DIR"
    echo "Created missing folder: $DEST_DIR" > /dev/tty
  fi

  cp -p "$SOURCE_PATH" "$DEST_PATH"
  if [ $? -eq 0 ]; then
    echo "Copied/Overwritten successfully: $DEST_PATH" > /dev/tty
    echo "Justification: $JUSTIFICATION" > /dev/tty
  else
    echo "Error copying file $SOURCE_PATH -> $DEST_PATH" > /dev/tty
  fi

done < "$TMPFILE"

# Cleanup
rm -f "$TMPFILE"
