#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_examples.py
====================

**Purpose:**
-----------
This script processes a list of Python modules to extract usage examples from their `__main__` sections. It generates two primary outputs:
1. An HTML file (`class_examples.html`) containing formatted usage examples for each module.
2. A JSON file (`modules_details.json`) containing metadata about each processed module, including paths, anchors, and the number of examples.

This script is called by `utils/generate_diagrams.sh`.
The outputs are used by `utils/generate_mermaid.py`.
The whole ecosystem is called by `generate_diagrams.sh` (use this script to update the entire documentation).

**Functionality:**
------------------
- **Module Validation:** Ensures each module file exists and contains a `__main__` section.
- **Example Extraction:** Parses the `__main__` section to extract descriptive comments and corresponding code blocks.
- **HTML Generation:** Formats the extracted examples into collapsible HTML sections for easy viewing.
- **JSON Metadata:** Compiles module details into a structured JSON file for integration with other documentation tools.

**Usage:**
---------
Run the script from the `SFPPy/utils/` directory using the following command:

```bash
./generate_examples.py output_examples.html modules_details.json < modules_withexamples_list.txt
```

Use this command to generate an input file from `utils/`:

```bash
mainfolder="$(realpath ../)"
find "$mainfolder/patankar" "$mainfolder/patankar/private" -maxdepth 1 -type f -name '*.py' | sed "s|$mainfolder/||" | sort > ../html/modules_withexamples_list.txt
```

For production use this minimal example (see generate_diagrams.sh)
```bash
    mainfolder="$(realpath ../)"
    output_markdown="$mainfolder/html/patanakar_classes_documentation.md"
    output_dir="$(dirname "$output_markdown")"
    moduleexamplesList="modules_withexamples_list.txt"
    moduleexamplesDetails="class_examples_details.json"
    moduleexamplesHTML="class_examples.html"
    lookfolders=(
        "$mainfolder/patankar"
        "$mainfolder/patankar/private"
    )
    find "${lookfolders[@]}" -maxdepth 1 -type f -name '*.py' | sed "s|$mainfolder/||" | sort > "$output_dir/$moduleexamplesList"
    ./generate_examples.py "$output_dir/$moduleexamplesHTML" "$output_dir/$moduleexamplesDetails" < "$output_dir/$moduleexamplesList"
```

**Arguments:**
------------
1. `output_examples.html`: Path to the output HTML file that will contain the usage examples.
2. `modules_details.json`: Path to the output JSON file that will store module metadata.
3. `< modules_withexamples_list.txt`: A file containing a list of module paths (one per line) to be processed.

**Inputs:**
---------
- **Standard Input (`STDIN`):** A list of module paths in the format `patankar.module_name.py` (e.g., `patankar.migration.py`).

**Outputs:**
----------
1. **HTML File (`output_examples.html`):** Contains usage examples with descriptions and code snippets for each module.
2. **JSON File (`modules_details.json`):** Stores metadata for each module, including:
   - `full_path`: Absolute file path.
   - `local_path`: Relative file path as provided.
   - `url_anchor`: HTML anchor for linking.
   - `number_of_examples`: Count of extracted examples.

**Dependencies:**
---------------
- **Python 3.x**
- **Optional:** `pygments` library for syntax highlighting. If not installed, the script will use minimal syntax highlighting.

**Error Handling:**
------------------
- **Missing Files:** Warns and skips modules that do not exist.
- **Missing `__main__` Sections:** Skips modules without a `__main__` section.
- **No Examples Found:** Skips modules without extractable examples.

**Author:**
---------
INRAE\Olivier Vitrac  
Email: olivier.vitrac@agroparistech.fr  
Last Revised: 2024-12-23

"""

import sys
import os
import re
import json
from collections import defaultdict
from datetime import datetime
from html import escape
from pathlib import Path

try:
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    USE_PYGMENTS = True
except ImportError:
    USE_PYGMENTS = False

def ensure_run_from_utils():
    """Ensure the script is run from SFPPy/utils/ directory by checking for pdocme.sh."""
    if not os.path.isfile("pdocme.sh"):
        print("Error: This script must be run from the SFPPy/utils/ directory.", file=sys.stderr)
        sys.exit(1)

def get_version():
    """Extract the version number of SFPPy from version_file."""
    mainfolder = os.path.realpath(os.path.join(".."))
    version_file = os.path.join(mainfolder, "utils", "VERSION.txt")
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

def read_module_list():
    """Reads module paths from STDIN, one per line."""
    modules = [line.strip() for line in sys.stdin if line.strip()]
    return modules

def file_has_main_section(file_content):
    """Check if the file contains a '__main__' section."""
    return re.search(r"if\s+__name__\s*==\s*['\"]__main__['\"]\s*:", file_content) is not None

def extract_examples(file_content):
    """
    Extracts examples from the `if __name__ == '__main__':` section.
    Returns a list of tuples: (description, code)
    """
    examples = []
    # Match the if __name__ == '__main__': block
    pattern = re.compile(r"if\s+__name__\s*==\s*['\"]__main__['\"]\s*:\s*(.*)", re.DOTALL)
    match = pattern.search(file_content)
    if not match:
        return examples  # No examples found

    main_block = match.group(1)

    lines = main_block.splitlines()
    idx = 0
    n = len(lines)

    while idx < n:
        # Skip empty lines
        while idx < n and lines[idx].strip() == '':
            idx += 1
        if idx >= n:
            break

        # Check for description: contiguous comment lines or docstrings
        desc_lines = []
        if lines[idx].strip().startswith('#') or lines[idx].strip().startswith(('"""', "'''")):
            # Collect all leading comment lines or a docstring
            while idx < n:
                line = lines[idx]
                stripped = line.strip()
                if stripped.startswith('#'):
                    desc_lines.append(re.sub(r'^#+\s*', '', stripped))
                    idx += 1
                elif stripped.startswith(('"""', "'''")):
                    delimiter = stripped[:3]
                    if stripped.count(delimiter) >= 2:
                        # Single-line docstring
                        doc_content = stripped.strip(delimiter)
                        desc_lines.append(doc_content)
                        idx += 1
                    else:
                        # Multiline docstring
                        doc_content = stripped.lstrip(delimiter)
                        desc_lines.append(doc_content)
                        idx += 1
                        while idx < n:
                            doc_line = lines[idx]
                            if doc_line.strip().endswith(delimiter):
                                doc_content = doc_line.strip().rstrip(delimiter)
                                desc_lines.append(doc_content)
                                idx += 1
                                break
                            else:
                                desc_lines.append(doc_line.strip())
                                idx += 1
                else:
                    break
            description = '\n'.join(desc_lines).strip()
        else:
            description = ''

        # Collect code lines until the next description or end
        code_lines = []
        while idx < n:
            line = lines[idx]
            stripped = line.strip()
            if stripped == '':
                idx += 1
                continue
            if stripped.startswith('#') or stripped.startswith(('"""', "'''")):
                # Possible new description starts here
                break
            else:
                code_lines.append(line)
                idx += 1
        code = '\n'.join(code_lines).strip()

        if code:
            examples.append((description, code))

    return examples

def organize_modules(module_paths):
    """
    Organizes modules into a nested dictionary based on their directory structure.
    module_paths: List of module paths relative to mainfolder, e.g., 'patankar/geometry.py'
    Returns a nested dictionary representing the folder structure.
    """
    tree = {}
    for path in module_paths:
        parts = Path(path).parts  # e.g., ('patankar', 'migration.py') or ('patankar', 'private', 'geometry.py')
        current = tree
        for part in parts[:-1]:  # All parts except the last (file)
            current = current.setdefault(part, {})
        current.setdefault('_files', []).append(parts[-1])
    return tree

def generate_nav_html(tree, parent_path=''):
    """
    Recursively generates the HTML for the navigation menu.
    """
    html = '<ul>\n'
    for key in sorted(tree.keys()):
        if key == '_files':
            for file in sorted(tree['_files']):
                module_relative_path = f"{parent_path}/{file}" if parent_path else file
                # Remove the .py extension for the anchor without appending '_examples'
                module_anchor = module_relative_path.replace('/', '_').replace('.py', '')
                html += f'<li class="file"><a href="#{module_anchor}" onclick="loadDoc(\'{module_anchor}\')">{file}</a></li>\n'
        else:
            new_path = f"{parent_path}/{key}" if parent_path else key
            html += f'<li class="folder">\n'
            html += f'<div class="folder-title" onclick="toggleFolder(this)">{key}</div>\n'
            html += f'<ul class="folder-content">\n'
            html += generate_nav_html(tree[key], new_path)
            html += '</ul>\n</li>\n'
    html += '</ul>\n'
    return html

def syntax_highlight(code):
    """Applies syntax highlighting to the code block."""
    if USE_PYGMENTS:
        formatter = HtmlFormatter(nowrap=True)
        return highlight(code, PythonLexer(), formatter)
    else:
        # Minimal syntax highlighting based on CSS classes
        # This is a very simplistic approach and can be improved
        code = escape(code)
        # Highlight keywords
        keywords = r'\b(False|class|finally|is|return|None|continue|for|lambda|try|' \
                   r'True|def|from|nonlocal|while|and|del|global|not|with|' \
                   r'as|elif|if|or|yield|assert|else|import|pass|break|except|' \
                   r'in|raise)\b'
        code = re.sub(keywords, r'<span class="keyword">\1</span>', code)
        # Highlight strings
        code = re.sub(r'(\'[^\']*\'|"[^"]*")', r'<span class="string">\1</span>', code)
        # Highlight comments (only inline comments)
        code = re.sub(r'(?<!#)#.*', r'<span class="comment">\g<0></span>', code)
        return code

def generate_example_html(module_anchor, examples):
    """
    Generates the HTML for a single module's examples.
    """
    html = f'<div id="{module_anchor}" class="module-content" style="display:none;">\n'
    html += '<button onclick="toggleAllCode(this)" class="collapsible">Toggle All Code Sections</button>\n'
    for idx, (desc, code) in enumerate(examples):
        if desc:
            # Replace line breaks with <br> for HTML
            desc_html = escape(desc).replace('\n', '<br>')
            html += f'<p>{desc_html}</p>\n'
        if code:
            highlighted_code = syntax_highlight(code)
            html += f'''
    <button type="button" class="collapsible">Show Code Example {idx + 1}</button>
    <div class="content">
    <pre class="code">{highlighted_code}</pre>
    </div>
    '''
    html += '</div>\n'
    return html

def generate_welcome_html(generation_datetime):
    """
    Generates the HTML for the welcome page.
    """
    welcome_text = f"""
<h2>Welcome to Usage Class Examples</h2>
<p>Select a module in the left menu to view usage examples. These examples are not for production and are automatically extracted from the main section of each module.</p>
<p>Back to the <a href='index.html'>Python'SFPPy documentation</a>.</p>
<p>When no module is selected, you see this welcome page. They are used to test classes with typical codes. The main section often serves as a testing script, example usage block, or self-contained test block. It's a way to demonstrate how the module's functionality works or to run simple unit tests and examples inline.</p>
<p>Generated on: {generation_datetime}</p>
"""
    return f'<div id="welcome" class="module-content" style="display:block;">\n{welcome_text}\n</div>\n'

def generate_full_html(nav_html, content_html, base_css, pygments_css, generation_datetime):
    """
    Combines all parts into the full HTML document.
    """
    html = f'''<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'>
<meta name='viewport' content='width=device-width, initial-scale=1.0'>
<title>SFPPy Usage/Class Examples</title>
<style>
{base_css}

{pygments_css}
</style>
</head>
<body>
<header>
    <!-- Toggle Sidebar Button -->
    <button class='toggle-btn' id='toggleSidebar' aria-label='Toggle Sidebar' aria-expanded='false'>
        <kbd>&#9776;</kbd>
    </button>
    <h1>SFPPy Documentation - Usage Examples</h1>
</header>
<div id='content'>
<div id='nav'>
<p><strong>Version:</strong> SFPPy v.{get_version()}</p>
<p><strong>Maintained by:</strong> INRAE\\olivier.vitrac@agroparistech.fr</p>
<hr>
{nav_html}
</div>
<div id='main'>
{generate_welcome_html(generation_datetime)}
{content_html}
</div>
</div>
<footer>
<p>&copy; {datetime.now().year} SFPPy Project. All rights reserved.</p>
</footer>
<script>
// Toggle visibility of folder contents
function toggleFolder(element) {{
    var content = element.nextElementSibling;
    if (content.style.display === "block") {{
        content.style.display = "none";
    }} else {{
        content.style.display = "block";
    }}
}}

// Load documentation into the main panel
function loadDoc(moduleId) {{
    var modules = document.getElementsByClassName('module-content');
    for (var i = 0; i < modules.length; i++) {{
        modules[i].style.display = 'none';
    }}
    var welcome = document.getElementById('welcome');
    welcome.style.display = 'none';
    var selected = document.getElementById(moduleId);
    if (selected) {{
        selected.style.display = 'block';
    }}
}}

// Toggle All Code Sections
function toggleAllCode(button) {{
    var moduleContent = button.parentElement;
    var collapsibles = moduleContent.querySelectorAll('.collapsible');
    var contents = moduleContent.querySelectorAll('.content');
    var shouldExpand = true;
    // Determine if we should expand or collapse based on the first content's display
    if (contents.length > 0 && (contents[0].style.display === 'block')) {{
        shouldExpand = false;
    }}
    collapsibles.forEach(function(collapsible, index) {{
        var content = collapsibles[index].nextElementSibling;
        if (shouldExpand) {{
            content.style.display = 'block';
            collapsible.textContent = 'Hide Code Example ' + (index + 1);
        }} else {{
            content.style.display = 'none';
            collapsible.textContent = 'Show Code Example ' + (index + 1);
        }}
    }});
    button.textContent = shouldExpand ? 'Collapse All Code Sections' : 'Expand All Code Sections';
}}

// Initialize all collapsible buttons and handle URL hash on page load
document.addEventListener("DOMContentLoaded", function() {{
    var coll = document.getElementsByClassName("collapsible");
    for (var i = 0; i < coll.length; i++) {{
        coll[i].addEventListener("click", function() {{
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "block") {{
                content.style.display = "none";
                this.textContent = this.textContent.replace('Hide', 'Show');
            }} else {{
                content.style.display = "block";
                this.textContent = this.textContent.replace('Show', 'Hide');
            }}
        }});
    }}
    
    // Handle URL hash on page load to display the corresponding module examples
    var hash = window.location.hash.substring(1); // Remove the #
    if (hash) {{
        loadDoc(hash);
    }}
}});

// Toggle Sidebar Functionality
const toggleButton = document.getElementById('toggleSidebar');
const nav = document.getElementById('nav');

toggleButton.addEventListener('click', () => {{
    nav.classList.toggle('collapsed');
    document.body.classList.toggle('nav-open'); // Toggle overlay on small screens
    // Change icon based on sidebar state
    if(nav.classList.contains('collapsed')) {{
        toggleButton.innerHTML = '<kbd>&#9776;</kbd>'; // Hamburger icon
        toggleButton.setAttribute('aria-expanded', 'false');
    }} else {{
        toggleButton.innerHTML = '<kbd>&#10005;</kbd>'; // Close icon (X)
        toggleButton.setAttribute('aria-expanded', 'true');
    }}
}});
</script>
</body>
</html>
'''
    return html

def main():
    ensure_run_from_utils()

    if len(sys.argv) != 3:
        print("Usage: generate_examples.py output_examples.html modules_details.json < modules_withexamples_list.txt", file=sys.stderr)
        sys.exit(1)

    output_html = sys.argv[1]
    output_json = sys.argv[2]
    modules_with_paths = read_module_list()
    if not modules_with_paths:
        print("No modules provided.", file=sys.stderr)
        sys.exit(1)

    # Configuration
    mainfolder = os.path.realpath(os.path.join(".."))

    # Validate and filter modules
    valid_module_paths = []
    modules_details = {}
    for module_path in modules_with_paths:
        full_path = os.path.join(mainfolder, module_path)
        if not os.path.isfile(full_path):
            print(f"Warning: File '{module_path}' does not exist. Skipping.", file=sys.stderr)
            continue
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
        if not file_has_main_section(content):
            # Discard files without a '__main__' section
            print(f"Info: File '{module_path}' does not contain a '__main__' section. Skipping.", file=sys.stderr)
            continue
        examples = extract_examples(content)
        if not examples:
            print(f"Info: No examples found in '{module_path}'. Skipping.", file=sys.stderr)
            continue
        valid_module_paths.append(module_path)
        
        # === MODIFICATION START ===
        # Replace slashes with dots and remove .py to match module naming
        module_name = module_path.replace('/', '.').replace('.py', '')
        
        # Generate an anchor by replacing dots with underscores
        module_anchor = module_name.replace('.', '_')
        
        # Update the modules_details dictionary with the module name as the key
        modules_details[module_name] = {
            "full_path": full_path,
            "local_path": module_path,
            "url_anchor": f"#{module_anchor}",
            "number_of_examples": len(examples)
        }
        # === MODIFICATION END ===

    if not valid_module_paths:
        print("No valid modules with '__main__' sections and examples to process.", file=sys.stderr)
        sys.exit(1)

    # Organize modules into folder structure
    module_tree = organize_modules(valid_module_paths)

    # Generate navigation HTML
    nav_html = generate_nav_html(module_tree)

    # Process each module and generate content HTML
    content_html = ''
    for module_path in valid_module_paths:
        full_path = os.path.join(mainfolder, module_path)
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
        examples = extract_examples(content)
        if not examples:
            continue
        # === MODIFICATION START ===
        # Convert module_path to module_name to retrieve details
        module_name = module_path.replace('/', '.').replace('.py', '')
        module_detail = modules_details.get(module_name)
        if not module_detail:
            continue
        module_anchor = module_detail["url_anchor"].lstrip('#')
        # === MODIFICATION END ===
        
        content_html += generate_example_html(module_anchor, examples)

    # Get current date and time
    generation_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Define base CSS
    base_css = """
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

    # Initialize pygments_css
    pygments_css = ""
    if USE_PYGMENTS:
        # Generate PYGMENTS CSS styles
        formatter = HtmlFormatter()
        pygments_css = formatter.get_style_defs('.code')
    else:
        print("Warning: PYGMENTS is not installed. Code will have minimal syntax highlighting.", file=sys.stderr)

    # Generate the full HTML with both base CSS and Pygments CSS
    full_html = generate_full_html(nav_html, content_html, base_css, pygments_css, generation_datetime)

    # Write the HTML to the output file
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(full_html)

    # Write the modules details dictionary to a JSON file
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(modules_details, f, indent=4)

    print(f"HTML documentation generated successfully at '{output_html}'")
    print(f"Module details dictionary saved at '{output_json}'")

if __name__ == '__main__':
    main()