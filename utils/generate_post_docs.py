#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_post_docs.py
=====================

**Purpose:**
-----------
Automate the generation of the POST examples documentation page for Pizza3. The script performs the following tasks:
1. Archives the existing `post` folder by zipping it with maximum compression and appending the creation/modification date to its name. The ZIP is stored in `history/`.
2. Copies and renames HTML and PDF files from the source directory to the destination directory (`html/post/`).
3. Generates a JSON manifest pairing the renamed HTML and PDF files.
4. Extracts synopses from the HTML files.
5. Generates a comprehensive HTML documentation page with a left navigation menu and a right content panel.
6. Incorporates the provided CSS styles to maintain a consistent layout.
7. Appends contact information and the documentation creation date at the end of the right panel.

- **`FILE_MAPPING`**: Dictionary mapping the modified names to their original HTML and PDF filenames.

```python
FILE_MAPPING = {
    "POST_example1": ("example1.html", "Pizza3 - WORKSHOP - PostTreatment - Part 1.pdf"),
    "POST_example2": ("example2.html", "Pizza3 - WORKSHOP - PostTreatment - Part 2.pdf"),
    "POST_example3": ("example2bis.html", "Pizza3 - WORKSHOP - PostTreatment - Part 2bis.pdf"),
}
```

| Modified Name    | Original HTML Filename         | Original PDF Filename                                    |
|------------------|--------------------------------|----------------------------------------------------------|
| `POST_example1`  | `example1.html`                | `Pizza3 - WORKSHOP - PostTreatment - Part 1.pdf`         |
| `POST_example2`  | `example2.html`                | `Pizza3 - WORKSHOP - PostTreatment - Part 2.pdf`         |
| `POST_example3`  | `example2bis.html`             | `Pizza3 - WORKSHOP - PostTreatment - Part 2bis.pdf`      |

**Usage:**
----------
Run the script from the `utils/` directory:
```bash
python3 generate_post_docs.py
```

**Dependencies:**
---------------
- Python 3.x
- `beautifulsoup4` library

**Author:**
---------
INRAE\Olivier Vitrac  
Email: olivier.vitrac@agroparistech.fr  
Last Revised: 2025-01-09

"""

import os
import shutil
import json
from pathlib import Path
from datetime import datetime
from bs4 import BeautifulSoup
import zipfile
import sys
import re

# Root folder and version file
MAINFOLDER = Path("../").resolve()
version_file = os.path.join(MAINFOLDER, "utils", "VERSION.txt")

# Configuration
SOURCE_DIR = MAINFOLDER / "post" / "html"
DEST_DIR = MAINFOLDER / "html" / "post"
HISTORY_DIR = MAINFOLDER / "history"
HTML_POST_FOLDER = MAINFOLDER / "html" / "post"

# Mapping: Modified Name -> (Original HTML, Original PDF)
FILE_MAPPING = {
    "POST_example1": ("example1.html", "Pizza3 - WORKSHOP - PostTreatment - Part 1.pdf"),
    "POST_example2": ("example2.html", "Pizza3 - WORKSHOP - PostTreatment - Part 2.pdf"),
    "POST_example3": ("example2bis.html", "Pizza3 - WORKSHOP - PostTreatment - Part 2bis.pdf"),
}

# Provided CSS
BASE_CSS = """
body {
    font-family: 'Segoe UI', Arial, sans-serif; 
    margin: 0; 
    padding: 0; 
    line-height: 1.6; 
    background-color: #f9f9f9; 
    color: #333;
}
header {
    background: #4CAF50; 
    color: #fff; 
    padding: 10px;
    position: relative; /* For positioning the toggle button */
}
header h1 {
    margin: 0; 
    font-size: 1.5em;
    color: #fff; /* Explicitly set to white */
    padding-left: 50px; /* Space for the toggle button */
}
#content {
    display: flex;
    height: calc(100vh - 50px); /* Adjusted for header height */
    transition: all 0.3s ease;
}
#nav {
    min-width: 250px;
    max-width: 250px; /* Fixed width */
    background: #fff;
    border-right: 1px solid #ddd;
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
    transition: all 0.3s ease;
}
#nav.collapsed {
    margin-left: -250px; /* Hide the sidebar */
}
#main {
    flex: 1;
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
    transition: all 0.3s ease;
}
header .toggle-btn {
    position: absolute;
    top: 50%;
    transform: translateY(-50%); /* Center the button vertically */
    left: 10px; /* Place the button on the left */
    background-color: #4CAF50; /* Green background */
    border: none;
    color: white; /* Ensure the hamburger icon is white */
    padding: 10px 12px; /* Adjust padding for larger button */
    cursor: pointer;
    font-size: 1.2em; /* Increase font size for better visibility */
    border-radius: 4px;
    z-index: 1001; /* Ensure the button is above other elements */
}
header .toggle-btn:hover {
    background-color: #45a049;
}
header .toggle-btn kbd {
    font-family: 'Arial', sans-serif; /* Match the header font */
    color: white; /* Ensure the hamburger icon is white */
    font-size: 1.2em; /* Same size as the button text */
    background: none; /* Remove any background styling from <kbd> */
    border: none; /* Remove any borders from <kbd> */
}
h1 {
    font-size: 1.8em;
    color: #333;
}
h2 {
    color: #333; 
    border-bottom: 2px solid #4CAF50; 
    padding-bottom: 5px;
}
a {
    text-decoration: none; 
    color: #007BFF;
}
a:hover {
    text-decoration: underline;
}
ul {
    list-style-type: none; 
    padding-left: 0;
    margin: 0;
}
li {
    margin: 5px 0;
}
.folder-title {
    font-weight: bold; 
    color: #333; 
    padding: 5px 0;
    cursor: pointer;
}
.folder-content {
    margin-left: 20px;
    display: none; /* start collapsed */
}
.file {
    margin-left: 20px;
}
hr {
    margin: 20px 0; 
    border: 1px solid #ddd;
}
footer {
    font-size: 0.9em; 
    color: #666; 
    margin-top: 20px; 
    text-align: center;
}
/* Enhanced Table Styling with Banded Colors */
table {
    border-collapse: collapse;
    width: 100%;
}
th, td {
    border: 1px solid #ddd;
    padding: 8px;
}
th {
    background-color: #4CAF50;
    color: white;
}
tr:nth-child(even) {
    background-color: #f2f2f2; /* Light gray for even rows */
}
tr:nth-child(odd) {
    background-color: rgba(76, 175, 80, 0.1); /* Light green for odd rows */
}
/* Collapsible Code Section */
.collapsible {
    background-color: #f1f1f1;
    color: #333;
    cursor: pointer;
    padding: 10px;
    width: 100%;
    border: none;
    text-align: left;
    outline: none;
    font-size: 1em;
}
.active, .collapsible:hover {
    background-color: #ddd;
}
.content {
    padding: 0 18px;
    display: none;
    overflow: hidden;
    background-color: #f9f9f9;
}
/* Minimal Syntax Highlighting */
.code {
    background-color: #f4f4f4;
    padding: 10px;
    border: 1px solid #ddd;
    overflow-x: auto;
    font-family: 'Courier New', Courier, monospace;
    color: #333;
}
.keyword {
    color: #007BFF;
    font-weight: bold;
}
.comment {
    color: #6a9955;
    font-style: italic;
}
.string {
    color: #a31515;
}

/* Responsive Design */
@media screen and (max-width: 768px) {
    #nav {
        position: absolute;
        left: 0;
        top: 50px; /* Height of the header */
        height: calc(100% - 50px);
        z-index: 1000;
    }
    #nav.collapsed {
        margin-left: -250px;
    }
    #main {
        flex: 1;
    }
}
"""

def get_version():
    """Extract the version number of Pizza3 from version_file."""
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

def print_header(title):
    """Prints a formatted header for better readability in the terminal."""
    print("\n" + "="*60)
    print(title)
    print("="*60 + "\n")

def verify_execution_directory():
    """Verifies that the script is being run from the 'utils' directory."""
    current_dir = Path.cwd().name
    if current_dir != "utils":
        print(f"Error: This script must be run from the 'utils' directory. Currently in '{current_dir}'.")
        sys.exit(1)
    else:
        print(f"Verified: Running from '{current_dir}' directory.")

def archive_existing_html_post_folder():
    """
    Archives the existing 'html/post' folder by zipping it with maximum compression,
    appending the creation/modification date to the ZIP filename, and saving it in 'history/'.
    """
    if not HTML_POST_FOLDER.exists():
        print("No existing 'html/post' folder found. Skipping archiving.")
        return

    # Ensure history directory exists
    if not HISTORY_DIR.exists():
        print(f"Creating history directory at '{HISTORY_DIR}'.")
        HISTORY_DIR.mkdir(parents=True)
        print("History directory created.")
    else:
        print(f"History directory '{HISTORY_DIR}' already exists.")

    # Define ZIP filename with current date
    current_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    zip_filename = f"post_backup_{current_date}.zip"
    zip_path = HISTORY_DIR / zip_filename

    # Create ZIP archive with maximum compression
    print(f"Archiving 'html/post' folder to '{zip_path}' with maximum compression.")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(HTML_POST_FOLDER):
            for file in files:
                file_path = Path(root) / file
                # Calculate the relative path to maintain folder structure in ZIP
                relative_path = file_path.relative_to(MAINFOLDER)
                zipf.write(file_path, arcname=relative_path)
    print(f"Archived 'html/post' folder successfully as '{zip_filename}'.")

    # Remove the existing 'html/post' folder
    print(f"Removing the existing 'html/post' folder at '{HTML_POST_FOLDER}'.")
    shutil.rmtree(HTML_POST_FOLDER)
    print("'html/post' folder removed successfully.")

    # Recreate a fresh 'html/post' folder
    print(f"Recreating a fresh 'html/post' folder at '{HTML_POST_FOLDER}'.")
    HTML_POST_FOLDER.mkdir(parents=True)
    print("'html/post' folder recreated successfully.")

def copy_and_rename_files():
    """
    Copies and renames HTML and PDF files from the source directory to the destination directory.
    Generates a JSON manifest pairing the renamed HTML and PDF files.
    """
    print_header("Copying and Renaming Files")

    # Create destination directory if it doesn't exist
    if not DEST_DIR.exists():
        print(f"Creating destination directory at '{DEST_DIR}'.")
        DEST_DIR.mkdir(parents=True)
        print("Destination directory created.")
    else:
        print(f"Destination directory '{DEST_DIR}' already exists.")

    manifest = []

    for mod_name, (html_orig, pdf_orig) in FILE_MAPPING.items():
        src_html = SOURCE_DIR / html_orig
        src_pdf = SOURCE_DIR / pdf_orig

        # Define new filenames
        new_html = DEST_DIR / f"{mod_name}.html"
        new_pdf = DEST_DIR / f"{mod_name}.pdf"

        # Check if source files exist
        if not src_html.exists():
            print(f"Warning: Source HTML file '{html_orig}' does not exist. Skipping '{mod_name}'.")
            continue
        if not src_pdf.exists():
            print(f"Warning: Source PDF file '{pdf_orig}' does not exist. Skipping '{mod_name}'.")
            continue

        # Copy and rename HTML
        shutil.copy2(src_html, new_html)
        print(f"Copied and renamed '{html_orig}' to '{new_html.name}'.")

        # Copy and rename PDF
        shutil.copy2(src_pdf, new_pdf)
        print(f"Copied and renamed '{pdf_orig}' to '{new_pdf.name}'.")

        # Append to manifest
        manifest.append({
            "modified_name": mod_name,
            "html_filename": new_html.name,
            "pdf_filename": new_pdf.name
        })

    # Save manifest to JSON
    manifest_path = DEST_DIR / "manifest.json"
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=4)
    print(f"Manifest file created at '{manifest_path}'.")

    return manifest

def extract_synopsis(html_file):
    """
    Extracts the synopsis from an HTML file.
    The synopsis is defined as the content between <h2>Synopsis</h2> and the next <h2> tag.
    """
    with open(html_file, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')

    # Find the <h2>Synopsis</h2> tag
    synopsis_header = soup.find('h2', string=lambda text: text and 'Synopsis' in text)
    if not synopsis_header:
        print(f"Warning: 'Synopsis' section not found in '{html_file.name}'.")
        return "<p>No synopsis available.</p>"

    # Initialize an empty BeautifulSoup object for the synopsis content
    synopsis_content = BeautifulSoup("", 'html.parser')

    # Iterate through the siblings after the synopsis header until the next <h2>
    for sibling in synopsis_header.find_next_siblings():
        if sibling.name == 'h2':
            break
        synopsis_content.append(sibling)

    return str(synopsis_content)

def generate_documentation(manifest):
    """
    Generates the final HTML documentation page with a left navigation menu and a right content panel.
    Incorporates the provided CSS and includes introductory text and contact details.
    """
    print_header("Generating Documentation HTML Page")

    # Start building the navigation menu
    nav_html = '<ul>\n'
    for item in manifest:
        mod_name = item["modified_name"]
        nav_html += f'    <li><a href="#{mod_name}">{mod_name}</a></li>\n'
    nav_html += '</ul>\n'

    # Start building the main content
    main_content = ''
    for item in manifest:
        mod_name = item["modified_name"]
        html_filename = item["html_filename"]
        pdf_filename = item["pdf_filename"]
        html_file = DEST_DIR / html_filename

        # Extract synopsis
        synopsis = extract_synopsis(html_file)

        # Build the HTML block for this example
        example_html = f'''
        <div id="{mod_name}" class="example-content" style="display:none;">
            <h2>{mod_name}</h2>
            <div class="synopsis">
                {synopsis}
            </div>
            <div class="links">
                <a href="{html_filename}" target="_blank">View HTML</a> | 
                <a href="{pdf_filename}" target="_blank">Download PDF</a>
            </div>
        </div>
        '''
        main_content += example_html

    # Define the introductory text with a link to ../index_matlab.html
    introductory_text = f"""
    <h1>Pizza3 POST Documentation Index</h1>
    <p>Select a POST example from the left panel to view its documentation. You can have access to the <a href="../index_matlab.html">POST tools documentation here</a>.</p>
    <p><strong>Version:</strong> Pizza3 v.{get_version()}</p>
    <p><strong>Maintained by:</strong> INRAE\\olivier.vitrac@agroparistech.fr</p>
    """

    # Get the current date and time for the footer
    creation_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Define the full HTML structure with the provided CSS
    full_html = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Pizza3 POST Documentation Examples</title>
        <style>
            {BASE_CSS}
        </style>
    </head>
    <body>
        <header>
            <!-- Toggle Sidebar Button -->
            <button class="toggle-btn" id="toggleSidebar" aria-label="Toggle Sidebar">
                <kbd>&#9776;</kbd>
            </button>
            <h1>Pizza3 POST Documentation Index</h1>
        </header>
        <div id="content">
            <div id="nav">
                {nav_html}
            </div>
            <div id="main">
                {introductory_text}
                {main_content}
                <hr>
                <p>Contact: INRAE\\olivier.vitrac@agroparistech.fr</p>
                <p>Documentation Created on: {creation_datetime}</p>
            </div>
        </div>
        <footer>
            <p>&copy; {datetime.now().year} Pizza3 Project. All rights reserved.</p>
        </footer>
        <script>
            // Function to handle navigation clicks
            const navLinks = document.querySelectorAll('#nav ul li a');
            navLinks.forEach(link => {{
                link.addEventListener('click', function(e) {{
                    e.preventDefault();
                    const targetId = this.getAttribute('href').substring(1);
                    const examples = document.querySelectorAll('.example-content');
                    examples.forEach(example => {{
                        example.style.display = 'none';
                    }});
                    const target = document.getElementById(targetId);
                    if(target) {{
                        target.style.display = 'block';
                        target.scrollIntoView({{behavior: "smooth"}});
                    }}
                    // Collapse sidebar on small screens after selection
                    if (window.innerWidth <= 768) {{
                        nav.classList.add('collapsed');
                        toggleButton.innerHTML = '<kbd>&#9776;</kbd>';
                    }}
                }});
            }});

            // Toggle Sidebar Functionality
            const toggleButton = document.getElementById('toggleSidebar');
            const nav = document.getElementById('nav');

            toggleButton.addEventListener('click', () => {{
                nav.classList.toggle('collapsed');
                // Change icon based on sidebar state
                if(nav.classList.contains('collapsed')) {{
                    toggleButton.innerHTML = '<kbd>&#9776;</kbd>'; // Hamburger icon
                }} else {{
                    toggleButton.innerHTML = '<kbd>&#10005;</kbd>'; // Close icon (X)
                }}
            }});

            // Handle URL hash on page load to display the corresponding module examples or welcome page
            window.addEventListener('load', function() {{
                const hash = window.location.hash.substring(1);
                const examples = document.querySelectorAll('.example-content');
                if(hash) {{
                    examples.forEach(example => {{
                        if(example.id === hash) {{
                            example.style.display = 'block';
                        }} else {{
                            example.style.display = 'none';
                        }}
                    }});
                }} else {{
                    // Show welcome content if no hash is present
                    document.querySelector('#main').scrollTop = 0;
                }}
            }});
        </script>
    </body>
    </html>
    '''

    # Write the HTML to the destination directory
    doc_path = DEST_DIR / "index_post.html"
    with open(doc_path, 'w', encoding='utf-8') as f:
        f.write(full_html)
    print(f"POST Documentation HTML page created at '{doc_path}'.")

def main():
    """Main function to orchestrate the documentation generation."""
    print_header("Starting POST Examples Documentation Generation")

    # Step 1: Verify Execution Directory
    verify_execution_directory()

    # Step 2: Archive Existing 'html/post' Folder if it exists
    archive_existing_html_post_folder()

    # Step 3: Copy and Rename Files, Generate Manifest
    manifest = copy_and_rename_files()

    if not manifest:
        print("No files were copied and renamed. Exiting.")
        sys.exit(1)

    # Step 4: Generate Documentation HTML Page
    generate_documentation(manifest)

    print_header("Documentation Generation Completed Successfully")

if __name__ == "__main__":
    main()