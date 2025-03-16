#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_wiki.py

Generates an index HTML page for the SFPPy project that uses an iframe
to display original Typora-generated HTML files located in the wiki folder.

This script scans the SFPPy project's HTML files (converted from Markdown) in the wiki folder,
copies them to a destination folder (html/wikipages) while preserving the subfolder structure,
renaming files if duplicates occur (appending a "~" to the filename), and builds a hierarchical,
collapsible left menu. It also parses each HTML file to look for asset references (from "assets/" or "./assets/")
and copies these files preserving their local folder structure.
The index file is written in the destination folder and refers to the HTML files from this folder.
Initially, the iframe loads wiki.html.

Usage:
    Ensure that this script is run from the `SFPPy/utils/` directory.

    ```bash
    python generate_wiki.py
    ```

Output:
    - The `html/wikipages` directory will contain the copied HTML files, any referenced assets,
      and the generated `index.html` file.

Author:
    - **INRAE\Olivier Vitrac**
    - **Email:** olivier.vitrac@agroparistech.fr
    - **Last Revised:** 2025-02-12

Version:
    SFPPy v.1.30
"""

import os
import re
import sys
import shutil
from datetime import datetime
import html

# Ensure the script is run from SFPPy/utils/
if not os.path.isfile("pdocme.sh"):
    print("Error: This script must be run from the SFPPy/utils/ directory.")
    sys.exit(1)

# Root folder and version file
mainfolder = os.path.realpath(os.path.join(".."))
version_file = os.path.join(mainfolder, "utils", "VERSION.txt")
def get_version():
    """Extract the version number of SFPPy from version_file."""
    if not os.path.isfile(version_file):
        sys.stderr.write(f"Error: {version_file} not found. Please create a file with content: version=\"XX.YY.ZZ\"\n")
        sys.exit(1)
    with open(version_file, "r") as f:
        for line in f:
            line = line.strip()
            match = re.match(r'^version\s*=\s*"(.*?)"$', line)
            if match:
                return match.group(1)
    sys.stderr.write(f"Error: No valid version string found in {version_file}. Ensure it contains: version=\"XX.YY.ZZ\"\n")
    sys.exit(1)

# Configuration
# Source folder: mainfolder/wiki
source_dir = os.path.join(mainfolder, "wiki")
# Destination folder: mainfolder/html/wikipages
dest_dir = os.path.join(mainfolder, "html", "wikipages")
output_file = "index.html"
SFPPy_VERSION = f"SFPPy v.{get_version()}"
CONTACT = "INRAE\\olivier.vitrac@agroparistech.fr"

# Create destination folder if it doesn't exist
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# CSS Style for layout, iframe responsiveness, and collapsible menu styling
CSS_STYLE = """
body {
    font-family: 'Segoe UI', Arial, sans-serif;
    margin: 0;
    padding: 0;
    background-color: #f9f9f9;
    color: #333;
}
header {
    background: #4CAF50;
    color: #fff;
    padding: 10px;
    position: relative;
}
header h1 {
    margin: 0;
    font-size: 1.5em;
    padding-left: 50px;
}
#content {
    display: flex;
    height: calc(100vh - 50px);
}
#nav {
    width: 300px;
    background: #fff;
    border-right: 1px solid #ddd;
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
    flex-shrink: 0;
}
#nav.collapsed {
    width: 0;
    padding: 20px 0;
}
#main {
    flex: 1;
    padding: 0;
    overflow-y: auto;
    box-sizing: border-box;
}
#contentFrame {
    width: 100%;
    height: 100%;
    border: none;
}
.toggle-btn {
    position: absolute;
    top: 50%;
    left: 10px;
    transform: translateY(-50%);
    background-color: #4CAF50;
    border: none;
    color: white;
    padding: 10px 12px;
    cursor: pointer;
    font-size: 1.2em;
    border-radius: 4px;
    z-index: 1001;
}
.toggle-btn:hover {
    background-color: #45a049;
}
.folder-title {
    font-weight: bold;
    cursor: pointer;
    margin-bottom: 5px;
}
ul {
    list-style-type: none;
    padding-left: 15px;
    margin: 0;
}
li {
    margin: 5px 0;
}
a {
    text-decoration: none;
    color: #007BFF;
}
a:hover {
    text-decoration: underline;
}
@media screen and (max-width: 768px) {
    #nav {
        position: absolute;
        left: 0;
        top: 50px;
        height: calc(100% - 50px);
        z-index: 1000;
    }
}
"""

# JavaScript for toggling sidebar, folders, and loading the iframe
JS_SCRIPT = """
// Function to toggle a folder's visibility
function toggleFolder(el) {
    var content = el.nextElementSibling;
    if (content.style.display === "none" || content.style.display === "") {
        content.style.display = "block";
    } else {
        content.style.display = "none";
    }
}

// Toggle Sidebar Functionality
const toggleButton = document.getElementById('toggleSidebar');
const nav = document.getElementById('nav');
toggleButton.addEventListener('click', () => {
    nav.classList.toggle('collapsed');
    if(nav.classList.contains('collapsed')) {
        toggleButton.innerHTML = '<kbd>&#9776;</kbd>';
        toggleButton.setAttribute('aria-expanded', 'false');
    } else {
        toggleButton.innerHTML = '<kbd>&#10005;</kbd>';
        toggleButton.setAttribute('aria-expanded', 'true');
    }
});

// Load file into the iframe
function loadDoc(url) {
    document.getElementById('contentFrame').src = url;
    // Optionally update the URL hash for bookmarking
    window.location.hash = url;
}

// On page load, if a hash is present, load that file in the iframe; otherwise load wiki.html
window.addEventListener('load', function() {
    const hash = window.location.hash.substring(1);
    if(hash) {
        document.getElementById('contentFrame').src = hash;
    } else {
        document.getElementById('contentFrame').src = "wiki.html";
    }
});
"""

# Exclusion lists for copying files
excluded_dirs = [
    "__pycache__",
    ".ipynb_checkpoints",
    "anaconda_projects",
    "assets",        # Note: we do not copy an entire assets directory unless referenced by an HTML.
    "history",
    "help",
    "debug",
    "tmp",
    "docs"
]

excluded_files = [
    "index.html",
    "README.html"
    # Other exclusions can be added as needed.
    # Note: We are now copying wiki.html.
]

def is_excluded(path):
    for d in excluded_dirs:
        if os.path.sep + d + os.path.sep in path:
            return True
    if os.path.basename(path) in excluded_files:
        return True
    return False

# Set to track copied assets (destination paths)
copied_assets = set()

# Copy HTML files from source_dir to dest_dir while preserving subfolder structure.
# Rename files if duplicates occur.
# After copying an HTML file, search its content for asset references (assets/ or ./assets/) and copy them.
asset_pattern = re.compile(r'(?:src|href)\s*=\s*["\'](\.?\/?assets\/[^"\']+)["\']', re.IGNORECASE)

for root, dirs, files in os.walk(source_dir):
    # Remove directories in the exclusion list
    dirs[:] = [d for d in dirs if d not in excluded_dirs]
    for f in files:
        if f.endswith(".html"):
            source_file = os.path.join(root, f)
            if is_excluded(source_file):
                continue
            rel_path = os.path.relpath(source_file, source_dir)
            dest_file = os.path.join(dest_dir, rel_path)
            dest_folder = os.path.dirname(dest_file)
            if not os.path.exists(dest_folder):
                os.makedirs(dest_folder)
            # If destination file already exists, rename the new file by appending "~"
            if os.path.exists(dest_file):
                dest_file = dest_file + "~"
            shutil.copy2(source_file, dest_file)
            
            # Now, open the source HTML file and search for asset references
            with open(source_file, "r", encoding="utf-8", errors="replace") as fin:
                content = fin.read()
            for match in asset_pattern.finditer(content):
                asset_rel_path = match.group(1)
                # Normalize the asset path (remove leading "./" if any)
                asset_rel_path = asset_rel_path.lstrip("./")
                # Determine the source asset file path relative to the current HTML file
                html_dir = os.path.dirname(source_file)
                asset_source = os.path.join(html_dir, asset_rel_path)
                if not os.path.exists(asset_source):
                    continue  # asset not found; could log a warning here if desired
                # Determine destination: relative to the destination HTML file's folder.
                dest_asset = os.path.join(os.path.dirname(dest_file), asset_rel_path)
                # Create destination folder for the asset if needed
                asset_dest_folder = os.path.dirname(dest_asset)
                if not os.path.exists(asset_dest_folder):
                    os.makedirs(asset_dest_folder)
                # Apply duplicate renaming rule if asset already exists and hasn't been copied
                if dest_asset in copied_assets or os.path.exists(dest_asset):
                    dest_asset = dest_asset + "~"
                shutil.copy2(asset_source, dest_asset)
                copied_assets.add(dest_asset)

# Now, build a nested tree representing the folder structure in dest_dir.
tree = {}
for root, dirs, files in os.walk(dest_dir):
    dirs[:] = [d for d in dirs]
    for f in files:
        if f.endswith(".html"):
            if f == output_file:
                continue  # skip the generated index
            full_path = os.path.join(root, f)
            rel = os.path.relpath(full_path, dest_dir)
            parts = rel.split(os.path.sep)
            current = tree
            for part in parts[:-1]:
                current = current.setdefault(part, {})
            current.setdefault("_files", []).append((os.path.splitext(parts[-1])[0], rel))

def generate_nav_html(tree):
    """Recursively generate HTML for the navigation tree."""
    html_nav = "<ul>\n"
    if "_files" in tree:
        for fname, relpath in sorted(tree["_files"], key=lambda x: x[0].lower()):
            url = html.escape(relpath)
            html_nav += f"<li><a href='javascript:loadDoc(\"{url}\")'>{html.escape(fname)}</a></li>\n"
    for key in sorted(tree.keys()):
        if key == "_files":
            continue
        html_nav += f"<li><div class='folder-title' onclick='toggleFolder(this)'>{html.escape(key)}</div>\n"
        html_nav += generate_nav_html(tree[key])
        html_nav += "</li>\n"
    html_nav += "</ul>\n"
    return html_nav

nav_html = generate_nav_html(tree)
nav_html = "<div class='folder-title' onclick='toggleFolder(this)'>root</div>\n" + generate_nav_html(tree)

# Generate index.html with an iframe in the main panel in the destination folder.
index_file = os.path.join(dest_dir, output_file)
with open(index_file, "w", encoding="utf-8") as fout:
    fout.write("<!DOCTYPE html>\n<html lang='en'>\n<head>\n")
    fout.write("<meta charset='UTF-8'>\n")
    fout.write("<meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
    fout.write("<title>SFPPy Wiki Pages</title>\n")
    fout.write(f"<style>{CSS_STYLE}</style>\n")
    fout.write("</head>\n<body>\n")
    fout.write("<header>\n")
    fout.write("  <button class='toggle-btn' id='toggleSidebar' aria-label='Toggle Sidebar' aria-expanded='false'>\n")
    fout.write("    <kbd>&#9776;</kbd>\n")
    fout.write("  </button>\n")
    fout.write("  <h1>SFPPy Wiki Pages</h1>\n")
    fout.write("</header>\n")
    fout.write("<div id='content'>\n")
    fout.write("  <div id='nav'>\n")
    fout.write(f"    <p><strong>Version:</strong> {html.escape(SFPPy_VERSION)}</p>\n")
    fout.write(f"    <p><strong>Maintained by:</strong> {html.escape(CONTACT)}</p>\n")
    fout.write("    <hr>\n")
    fout.write(nav_html)
    fout.write("  </div>\n")
    fout.write("  <div id='main'>\n")
    fout.write("    <iframe id='contentFrame' src='wiki.html' title='Content'></iframe>\n")
    fout.write("  </div>\n")
    fout.write("</div>\n")
    fout.write(f"<script>{JS_SCRIPT}</script>\n")
    fout.write("</body>\n</html>")

print(f"Wiki pages generation completed. Output in {dest_dir}")
print(f"Index created at {index_file}")
