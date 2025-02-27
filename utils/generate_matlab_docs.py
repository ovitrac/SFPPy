#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_matlab_docs.py

Generates self-contained HTML documentation for the MATLAB functions of the Pizza3 project.

This script scans the Pizza3 project's MATLAB (.m) files, extracts documentation from their comments,
and compiles them into a structured, navigable HTML format. The generated documentation includes
function descriptions, examples, notes, and collapsible sections for MATLAB code with minimal
syntax highlighting.

Usage:
    Ensure that this script is run from the `Pizza3/utils/` directory. If not, the script will exit
    with an error message.

    ```bash
    python generate_matlab_docs.py
    ```

Output:
    - An `html` directory is created in the main project folder.
    - An `index_matlab.html` file is generated within the `html` directory, containing the compiled
      documentation.

Features:
    - **Directory and File Filtering**: Excludes specified directories and files to avoid processing
      unwanted content.
    - **Documentation Extraction**: Extracts and formats help documentation from MATLAB comment blocks.
    - **HTML Generation**: Creates a navigable HTML interface with a sidebar for function navigation
      and a main panel for content display.
    - **Styling and Interactivity**: Includes CSS for styling and JavaScript for interactive elements
      such as collapsible code sections and dynamic content loading.
    - **Syntax Highlighting**: Applies minimal syntax highlighting to MATLAB code snippets.

Configuration:
    - `mainfolder`: The main project directory path.
    - `output_dir`: Directory path for the generated HTML documentation.
    - `output_file`: Name of the main HTML file.
    - `PIZZA3_VERSION`: Version identifier for the documentation.
    - `CONTACT`: Contact information for the maintainer.
    - `CSS_STYLE`: CSS styles applied to the generated HTML.
    - `JS_SCRIPT`: JavaScript for interactive functionality in the HTML.

Dependencies:
    - Python 3.x
    - Standard Python libraries:
        - `os`
        - `re`
        - `sys`
        - `datetime`
        - `html`

Requirements:
    - The script must be executed from the `Pizza3/utils/` directory where `pdocme.sh` is present.

Example:
    After running the script, open the generated `html/index_matlab.html` in a web browser to view the
    documentation.

Notes:
    - The script ensures that all generated HTML elements have valid IDs to facilitate internal linking.
    - Collapsible sections enhance readability by allowing users to hide or show MATLAB code as needed.
    - The navigation menu mirrors the directory structure of the MATLAB files, providing an intuitive
      browsing experience.

Author:
    - **INRAE\Olivier Vitrac**
    - **Email:** olivier.vitrac@agroparistech.fr
    - **Last Revised:** 2025-01-17

Version:
    Pizza3 v.1.00

"""

import os
import re
import sys
from datetime import datetime
import html

# Ensure the script is run from Pizza3/utils/
if not os.path.isfile("pdocme.sh"):
    print("Error: This script must be run from the Pizza3/utils/ directory.")
    sys.exit(1)

# Root folder and version file
mainfolder = os.path.realpath(os.path.join(".."))
version_file = os.path.join(mainfolder, "utils", "VERSION.txt")
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

# Configuration
output_dir = os.path.join(mainfolder, "html")
output_file = "index_matlab.html"
PIZZA3_VERSION = f"Pizza3 v.{get_version()}"
CONTACT = "INRAE\\olivier.vitrac@agroparistech.fr"

# CSS Style with toggle button integration and dynamic sidebar collapse
CSS_STYLE = """
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
    transition: all 0.3s ease; /* Enable transitions for smooth animations */
}
#nav {
    width: 300px; /* Set a fixed width */
    background: #fff;
    border-right: 1px solid #ddd;
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
    transition: width 0.3s ease, padding 0.3s ease; /* Transition for smooth animations */
    flex-shrink: 0; /* Prevent flexbox from shrinking */
}
#nav.collapsed {
    width: 0; /* Hide the sidebar completely */
    padding: 20px 0; /* Optionally adjust padding */
}
#main {
    flex: 1;
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
    transition: all 0.3s ease; /* Enable transitions for smooth animations */
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
        width: 0; /* Hide the sidebar completely */
        padding: 20px 0; /* Adjust padding */
    }
    #main {
        flex: 1;
    }
    /* Add overlay when sidebar is open on mobile */
    body.nav-open::before {
        content: "";
        position: fixed;
        top: 50px;
        left: 0;
        width: 100%;
        height: calc(100% - 50px);
        background: rgba(0, 0, 0, 0.5);
        z-index: 999;
    }
}
"""

# JavaScript for toggle functionality and dynamic content display
JS_SCRIPT = """
// Toggle visibility of folder contents
function toggleFolder(el) {
    var content = el.nextElementSibling;
    if (content.style.display === "none" || content.style.display === "") {
        content.style.display = "block";
    } else {
        content.style.display = "none";
    }
};

// Load documentation into the main panel
function loadDoc(id) {
    // Hide all documentation sections
    var docs = document.getElementsByClassName('doc-content');
    for (var i = 0; i < docs.length; i++) {
        docs[i].style.display = 'none';
    }
    // Hide the welcome message
    var welcome = document.getElementById('welcome-message');
    if (welcome) {
        welcome.style.display = 'none';
    }
    // Show the selected documentation
    var selectedDoc = document.getElementById(id);
    if (selectedDoc) {
        selectedDoc.style.display = 'block';
    }
    
    // Reinitialize collapsible buttons within the displayed documentation
    var coll = selectedDoc.getElementsByClassName("collapsible");
    for (var i = 0; i < coll.length; i++) {
        // Remove existing event listeners to prevent multiple bindings
        coll[i].removeEventListener("click", toggleCollapsible);
        // Add event listener
        coll[i].addEventListener("click", toggleCollapsible);
    }
}

// Function to toggle collapsible sections
function toggleCollapsible() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
        content.style.display = "none";
    } else {
        content.style.display = "block";
    }
}

// Toggle Sidebar Functionality
const toggleButton = document.getElementById('toggleSidebar');
const nav = document.getElementById('nav');

toggleButton.addEventListener('click', () => {
    nav.classList.toggle('collapsed');
    document.body.classList.toggle('nav-open'); // Toggle overlay on small screens
    // Change icon based on sidebar state
    if(nav.classList.contains('collapsed')) {
        toggleButton.innerHTML = '<kbd>&#9776;</kbd>'; // Hamburger icon
        toggleButton.setAttribute('aria-expanded', 'false');
    } else {
        toggleButton.innerHTML = '<kbd>&#10005;</kbd>'; // Close icon (X)
        toggleButton.setAttribute('aria-expanded', 'true');
    }
});

// Handle URL hash on page load to display the corresponding documentation or welcome page
window.addEventListener('load', function() {
    const hash = window.location.hash.substring(1);
    const docs = document.querySelectorAll('.doc-content');
    if(hash) {
        docs.forEach(doc => {
            if(doc.id === hash) {
                doc.style.display = 'block';
            } else {
                doc.style.display = 'none';
            }
        });
        // If sidebar is open on small screens, ensure it's visible
        if (window.innerWidth <= 768 && !nav.classList.contains('collapsed')) {
            document.body.classList.add('nav-open');
        }
    } else {
        // Show welcome content if no hash is present
        const welcome = document.getElementById('welcome-message');
        if(welcome) {
            welcome.style.display = 'block';
        }
    }
});
"""

# Build a list of MATLAB files, excluding certain directories or files
excluded_dirs = [
    "fork",
    "history",
    "help",
    "debug",
    "fork",
    "tmp",
    "PIL",
    "restore",
    "__all__",
    "windowsONLY"
]

excluded_files = [
    "__init__.m",
    "__main__.m",
    "manifest.m",
    "debug.m"
]

def is_excluded(path):
    # Check if path contains any excluded directory
    for d in excluded_dirs:
        if os.path.sep + d + os.path.sep in path:
            return True
    # Check if file is excluded
    basename = os.path.basename(path)
    if basename in excluded_files:
        return True
    return False

matlab_files = []
for root, dirs, files in os.walk(mainfolder):
    # filter out excluded dirs
    dirs[:] = [d for d in dirs if d not in excluded_dirs]
    for f in files:
        if f.endswith(".m"):
            full_path = os.path.join(root, f)
            if not is_excluded(full_path):
                matlab_files.append(full_path)

matlab_files.sort()

# Function to extract documentation from a MATLAB file
def extract_matlab_doc(filepath):
    """
    Extracts MATLAB help from a .m file according to the standard MATLAB format.
    Enhancements:
    - Correctly captures the function name from the function definition.
    - If no function line is found, uses the filename as the function name.
    - Ensures that the entire help block is captured, not just until the first empty line.
    - Recognizes "NOTE:", "Note:", "note:", "EXAMPLES:", etc., and handles multi-line content.
    - Includes the MATLAB code in a collapsible section with minimal syntax highlighting.
    """
    doc_lines = []
    function_name = None
    inside_help = False
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    # First pass: Try to find the function line with improved regex
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("function"):
            # Improved regex to capture the function name correctly
            match = re.match(r'^function\s+(?:\[[^\]]*\]\s*=\s*)?(?:[A-Za-z0-9_]+\s*=\s*)?([A-Za-z0-9_]+)', stripped)
            if match:
                function_name = match.group(1)
            # Start looking for help after this line
            help_start_idx = idx + 1
            break
    else:
        # No function line found; use filename as function name
        function_name = os.path.splitext(os.path.basename(filepath))[0]
        help_start_idx = 0  # Start from the beginning of the file

    # Extract help lines
    for line in lines[help_start_idx:]:
        stripped = line.strip()
        if stripped.startswith("%"):
            # Remove the leading '%' and one space if present
            helpline = stripped.lstrip('%').lstrip()
            doc_lines.append(helpline)
            inside_help = True
        else:
            # If we were inside help and hit a non-% line, stop
            if inside_help:
                break

    # If no help was found after function line, and no function line exists, check for leading comments
    if not doc_lines and function_name == os.path.splitext(os.path.basename(filepath))[0]:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("%"):
                helpline = stripped.lstrip('%').lstrip()
                doc_lines.append(helpline)
            else:
                break  # Stop at first non-comment line

    # Now doc_lines contains the help block
    if not doc_lines:
        # No help available
        html_content = "<p>No help available.</p>"
    else:
        # Process help lines
        # The first line is considered the synopsis
        title_line = doc_lines[0]
        html_lines = []
        html_lines.append(f"<h1>{html.escape(title_line)}</h1>")

        # Define keywords and their corresponding HTML tags
        keywords = {
            "EXAMPLES:": "Examples",
            "EXAMPLE:": "Example",
            "SEE ALSO:": "See also",
            "CREDITS:": "Credits",
            "NOTE:": "Note",
            "NOTES:": "Notes",
            "AUTHOR:": "Author",
            "AUTHORS:": "Authors"
        }

        # Regular expression to detect keywords
        keyword_regex = re.compile(r'^(EXAMPLES?:|SEE ALSO:|CREDITS?:|NOTES?:|AUTHORS?:)', re.IGNORECASE)

        current_section = None
        section_content = []

        def flush_section():
            nonlocal current_section, section_content
            if current_section:
                if current_section.lower() == "see also":
                    # Process SEE ALSO links
                    links = []
                    for item in re.split(r'[,\s]+', ' '.join(section_content)):
                        if item in all_functions:
                            # Construct internal link
                            links.append(f'<a href="#{item}">{html.escape(item)}</a>')
                        else:
                            links.append(html.escape(item))
                    html_lines.append(f"<h2>{html.escape(current_section)}</h2>")
                    html_lines.append("<p>" + ", ".join(links) + "</p>")
                elif current_section.lower() in ["examples", "example"]:
                    html_lines.append(f"<h2>{html.escape(current_section)}</h2>")
                    # Preserve indentation and line breaks
                    html_lines.append("<pre class='code'><code class='language-matlab'>" + html.escape('\n'.join(section_content)) + "</code></pre>")
                elif current_section.lower() in ["note", "notes"]:
                    html_lines.append(f"<h2>{html.escape(current_section)}</h2>")
                    html_lines.append("<p>" + "<br/>".join(html.escape(line) for line in section_content) + "</p>")
                else:
                    html_lines.append(f"<h2>{html.escape(current_section)}</h2>")
                    html_lines.append("<p>" + "<br/>".join(html.escape(line) for line in section_content) + "</p>")
                section_content = []
            elif section_content:
                # Regular paragraph
                html_lines.append("<p>" + "<br/>".join(html.escape(line) for line in section_content) + "</p>")
                section_content = []

        # Collect all function names for linking in SEE ALSO
        all_functions = [os.path.splitext(os.path.basename(x))[0] for x in matlab_files]

        # Iterate over the help lines starting from the second line
        for line in doc_lines[1:]:
            if not line.strip():
                flush_section()
                continue
            keyword_match = keyword_regex.match(line)
            if keyword_match:
                flush_section()
                key = keyword_match.group(1).rstrip(':').capitalize()
                current_section = keywords.get(keyword_match.group(1).upper(), key)
                remainder = line[keyword_match.end():].strip()
                if remainder:
                    section_content.append(remainder)
            else:
                section_content.append(line)

        # Flush any remaining content
        flush_section()

        # Insert collapsible MATLAB code
        # Extract the MATLAB code from the file
        matlab_code = ''.join(lines)
        # Minimal Syntax Highlighting using CSS classes
        # Simple regex-based highlighting for keywords, comments, and strings
        def syntax_highlight(code):
            # Highlight comments
            code = re.sub(r'(%[^\n]*)', r'<span class="comment">\1</span>', code)
            # Highlight strings
            code = re.sub(r'(\'[^\']*\')', r'<span class="string">\1</span>', code)
            # Highlight keywords (a minimal set)
            keywords = ['function', 'end', 'if', 'else', 'elseif', 'for', 'while', 'switch', 'case', 'otherwise', 'return']
            pattern = r'\b(' + '|'.join(keywords) + r')\b'
            code = re.sub(pattern, r'<span class="keyword">\1</span>', code)
            return code

        highlighted_code = syntax_highlight(matlab_code)

        # Add collapsible section
        html_lines.append('<button class="collapsible">Show MATLAB Code</button>')
        html_lines.append('<div class="content"><pre class="code"><code class="language-matlab">' + highlighted_code + '</code></pre></div>')

        # Join all HTML lines
        html_content = "\n".join(html_lines)

    return function_name, html_content

# Extract docs from all .m files and generate HTML snippets
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Collect all function names for linking purposes
all_function_names = [os.path.splitext(os.path.basename(x))[0] for x in matlab_files]

# Create a structure (directory nav) similar to the Python version
file_structure = {}
for fpath in matlab_files:
    rel = os.path.relpath(fpath, mainfolder)
    d = os.path.dirname(rel)
    b = os.path.splitext(os.path.basename(rel))[0]
    if d not in file_structure:
        file_structure[d] = []
    file_structure[d].append(b)

# Sort directories and files
sorted_dirs = sorted(file_structure.keys())
for d in sorted_dirs:
    file_structure[d].sort()

# Generate individual HTML content for each function (hidden by default)
docs_html = ""
for fpath in matlab_files:
    function_name, html_doc = extract_matlab_doc(fpath)
    # Ensure that function_name is a valid HTML id (no spaces, special chars)
    valid_id = re.sub(r'\s+', '_', function_name)
    valid_id = re.sub(r'[^\w\-]', '', valid_id)
    # Create a div with id equal to valid_id
    docs_html += f"<div id='{html.escape(valid_id)}' class='doc-content' style='display: none;'>\n"
    docs_html += html_doc + "\n"
    docs_html += "</div>\n"

# Generate index_matlab.html
index_file = os.path.join(output_dir, output_file)
with open(index_file, "w", encoding="utf-8") as fout:
    fout.write("<!DOCTYPE html>\n<html lang='en'>\n<head>\n")
    fout.write("<meta charset='UTF-8'>\n")
    fout.write("<meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
    fout.write("<title>Pizza3 Matlab Documentation</title>\n")
    fout.write(f"<style>{CSS_STYLE}</style>\n")
    fout.write("</head>\n<body>\n")
    fout.write("<header>\n")
    fout.write("    <!-- Toggle Sidebar Button -->\n")
    fout.write("    <button class='toggle-btn' id='toggleSidebar' aria-label='Toggle Sidebar' aria-expanded='false'>\n")
    fout.write("        <kbd>&#9776;</kbd>\n")
    fout.write("    </button>\n")
    fout.write("    <h1>Pizza3 Documentation</h1>\n")
    fout.write("</header>\n")
    fout.write("<div id='content'>\n")
    fout.write("<div id='nav'>\n")
    fout.write(f"<p><strong>Version:</strong> {html.escape(PIZZA3_VERSION)}</p>\n")
    fout.write(f"<p><strong>Maintained by:</strong> {html.escape(CONTACT)}</p>\n")
    fout.write("<hr>\n")
    # Build the nav menu
    fout.write("<ul class='folder-list'>\n")
    for d in sorted_dirs:
        dirname = d if d != "." else "(root)"
        fout.write(f"<li class='folder'>\n")
        fout.write(f"<div class='folder-title' onclick='toggleFolder(this)'>{html.escape(dirname)}</div>\n")
        fout.write("<ul class='folder-content'>\n")
        for fname in file_structure[d]:
            # Ensure that fname corresponds to a valid_id
            valid_id = re.sub(r'\s+', '_', fname)
            valid_id = re.sub(r'[^\w\-]', '', valid_id)
            fout.write(f"<li class='file'><a href='#{html.escape(valid_id)}' onclick=\"loadDoc('{html.escape(valid_id)}')\">{html.escape(fname)}</a></li>\n")
        fout.write("</ul>\n</li>\n")
    fout.write("</ul>\n</div>\n")  # end nav
    fout.write("<div id='main'>\n")
    fout.write("<div id='welcome-message'>\n")
    fout.write("<h2>Welcome to Pizza3 Matlab Documentation</h2>\n")
    fout.write("<p>Select a function in the left menu to view its documentation.</p>\n")
    fout.write("<p>POST examples are fully detailed <a href='post/index_post.html' target='_blank'>here</a>.</p>\n")
    fout.write("<p>Back to the <a href='index.html'>Python'Pizza3 documentation</a>.</p>\n")
    fout.write("<hr>\n")
    fout.write("<p><i>When no function is selected, you see this welcome page.</i></p>\n")
    # Print the date
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
    fout.write(f"Generated on: {formatted_datetime}")
    fout.write("</div>\n")
    # Insert all documentation content inside #main
    fout.write(docs_html)
    fout.write("</div>\n")  # end main
    fout.write("</div>\n")  # end content
    # Insert JavaScript at the end of the body
    fout.write("<script>")
    fout.write(JS_SCRIPT)
    fout.write("</script>\n")
    fout.write("</body>\n</html>")

print(f"Documentation generation completed. Output in {output_dir}")
print(f"Index created at {index_file}")
