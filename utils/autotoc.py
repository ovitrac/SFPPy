#!/usr/bin/env python3
"""
autotoc.py

Adds a left side panel with a Table of Contents (TOC) to a given HTML file (for instance,
a Jupyter-exported HTML page). The TOC is built from heading tags (h1 up to hN) where N is a
configurable maximum level (default: 4, but can go up to 6 if needed).

The output HTML file uses a layout and styling similar to the SFPPy Wiki project.

Usage:
    python autotoc.py input.html output.html [max_level]

References:
    - Project: SFPPy Wiki Pages
    - Maintained by: INRAE\olivier.vitrac@agroparistech.fr
    - Version: Retrieved from ../utils/VERSION.txt (run from SFPPy/utils/)

    
"""

import os
import re
import sys
import html
from bs4 import BeautifulSoup

# Constants for SFPPy project info
PROJECT_NAME = "🍏⏩🍎 SFPPy Notebooks"
CONTACT = "INRAE\\olivier.vitrac@agroparistech.fr"
REPO = "https://github.com/ovitrac/SFPPy"
BADGE = "https://img.shields.io/badge/GitHub-SFPPy-4CAF50?style=for-the-badge&logo=github"
EMAIL = "olivier.vitrac@agmail.com"


def get_version():
    """
    Extract the version number of SFPPy from the VERSION.txt file located in ../utils/VERSION.txt.
    """
    base_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
    version_file = os.path.join(base_dir, "utils", "VERSION.txt")
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

VERSION = get_version()
SFPPy_VERSION = f"SFPPy v.{VERSION}"

# SFPPy infos
miniheader = f"""
  <div style="display: flex; align-items: center; gap: 12px;">
    <a href="{REPO}" target="_blank">
      <img src="{BADGE}"
           alt="GitHub SFPPy" style="border-radius: 8px;">
    </a>
    <div style="display: flex; align-items: center; font-size: 14px; font-weight: bold;">
      <span style="color: #4CAF50;">SFPPy v{VERSION}</span>
      <a href="mailto:{EMAIL}" title="E-mail the author: Olivier Vitrac" style="margin-left: 8px; font-size: 20px;">📩</a>
    </div>
  </div>
"""

def slugify(text, existing_ids):
    """
    Generate a URL-friendly slug for a given text.
    If the slug already exists in existing_ids, append a suffix.
    """
    slug = re.sub(r'\s+', '-', text.lower())
    slug = re.sub(r'[^a-z0-9\-]', '', slug)
    orig_slug = slug
    count = 1
    while slug in existing_ids:
        slug = f"{orig_slug}-{count}"
        count += 1
    existing_ids.add(slug)
    return slug

def generate_nested_toc(headers):
    """
    Generate nested TOC HTML from a list of headers.
    headers is a list of tuples: (level, text, id)
    """
    toc_html = ""
    current_level = 0
    for level, text, hid in headers:
        while current_level < level:
            toc_html += "<ul>\n"
            current_level += 1
        while current_level > level:
            toc_html += "</ul>\n"
            current_level -= 1
        toc_html += f"<li><a href='#{html.escape(hid)}'>{html.escape(text)}</a></li>\n"
    while current_level > 0:
        toc_html += "</ul>\n"
        current_level -= 1
    return toc_html

def add_toc_and_layout(html_content, max_level=4):
    """
    Parses the HTML content, generates a Table of Contents (TOC) from heading tags (h1 to h{max_level}),
    and injects a left side panel with the TOC, project info, and layout styling.
    """
    soup = BeautifulSoup(html_content, "html.parser")

    # Cap max_level to 6 if user passes a higher value.
    max_level = min(max_level, 6)
    header_tags = [f"h{i}" for i in range(1, max_level+1)]
    headers = []
    existing_ids = set()
    for header in soup.find_all(header_tags):
        # Get header text and remove trailing "¶" if present.
        header_text = header.get_text().strip()
        if header_text.endswith("¶"):
            header_text = header_text[:-1].strip()
        try:
            level = int(header.name[1])
        except ValueError:
            continue
        if level > max_level:
            continue
        # Ensure header has an id for linking.
        if not header.has_attr("id") or not header["id"]:
            header_id = slugify(header_text, existing_ids)
            header["id"] = header_id
        else:
            header_id = header["id"]
            existing_ids.add(header_id)
        headers.append((level, header_text, header_id))

    toc_html = "<div id='toc'><h2>Table of Contents</h2>\n" + generate_nested_toc(headers) + "</div>\n"

    # Create a header element with a toggle button and project title.
    header_elem = soup.new_tag("header")
    toggle_btn = soup.new_tag("button", id="toggleSidebar", **{
        "class": "toggle-btn",
        "aria-label": "Toggle Sidebar",
        "aria-expanded": "false"
    })
    toggle_btn.string = "\u2630"  # Hamburger icon
    header_elem.append(toggle_btn)
    h1_title = soup.new_tag("h1")
    h1_title.string = PROJECT_NAME
    header_elem.append(h1_title)

    # Create nav (left side panel) element and insert version and contact info.
    nav_elem = soup.new_tag("div", id="nav")
    info_html_old = (
        f"<p><strong>Version:</strong> {html.escape(SFPPy_VERSION)}</p>\n"
        f"<p><strong>Maintained by:</strong> {html.escape(CONTACT)}</p>\n<hr>\n"
    )
    info_html = miniheader
    nav_elem.append(BeautifulSoup(info_html, "html.parser"))
    nav_elem.append(BeautifulSoup(toc_html, "html.parser"))

    # Create main content container and move the original body content into it.
    main_elem = soup.new_tag("div", id="main")
    if soup.body:
        original_body_contents = list(soup.body.contents)
        for element in original_body_contents:
            main_elem.append(element.extract())
    else:
        main_elem.append(soup)

    # Create a container to hold the nav and main panels.
    content_elem = soup.new_tag("div", id="content")
    content_elem.append(nav_elem)
    content_elem.append(main_elem)

    # Clear the body and add the header and content container.
    if not soup.body:
        body = soup.new_tag("body")
        soup.append(body)
    else:
        soup.body.clear()
    soup.body.append(header_elem)
    soup.body.append(content_elem)

    # Insert CSS style into the <head>.
    css_style = """
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
    padding: 20px;
    overflow-y: auto;
    box-sizing: border-box;
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
#toc ul {
    list-style-type: none;
    padding-left: 15px;
    margin: 0;
}
#toc li {
    margin: 5px 0;
}
#toc a {
    text-decoration: none;
    color: #007BFF;
}
#toc a:hover {
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
    style_tag = soup.new_tag("style", type="text/css")
    style_tag.string = css_style
    if soup.head:
        soup.head.append(style_tag)
    else:
        head_tag = soup.new_tag("head")
        head_tag.append(style_tag)
        soup.insert(0, head_tag)

    # Insert JavaScript for toggling the sidebar.
    js_script = """
document.addEventListener("DOMContentLoaded", function(){
    var toggleButton = document.getElementById("toggleSidebar");
    var nav = document.getElementById("nav");
    toggleButton.addEventListener("click", function(){
        nav.classList.toggle("collapsed");
        if(nav.classList.contains("collapsed")){
            toggleButton.innerHTML = '<kbd>&#9776;</kbd>';
            toggleButton.setAttribute('aria-expanded', 'false');
        } else {
            toggleButton.innerHTML = '<kbd>&#10005;</kbd>';
            toggleButton.setAttribute('aria-expanded', 'true');
        }
    });
});
"""
    script_tag = soup.new_tag("script")
    script_tag.string = js_script
    if soup.body:
        soup.body.append(script_tag)
    else:
        soup.append(script_tag)

    return str(soup)

def main():
    if len(sys.argv) < 3:
        print("Usage: python autotoc.py input.html output.html [max_level]")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    max_level = 4
    if len(sys.argv) >= 4:
        try:
            max_level = int(sys.argv[3])
        except ValueError:
            print("Invalid max_level value, using default 4.")
            max_level = 4

    with open(input_file, "r", encoding="utf-8") as f:
        html_content = f.read()

    new_html = add_toc_and_layout(html_content, max_level)
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(new_html)
    print(f"TOC and layout added. Output written to {output_file}")

if __name__ == "__main__":
    main()
