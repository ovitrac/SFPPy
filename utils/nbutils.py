#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collection of utilities to manage interactive notebooks

    Author: INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2025-03-17
"""

# %% Dependencies
import os, re, fnmatch
import ipywidgets as widgets
from IPython.display import display, HTML
from IPython import get_ipython


# %% Constants
author = "Olivier Vitrac"
repo = "https://github.com/ovitrac/SFPPy"
web = "https://ovitrac.github.io/SFPPy/"
email = "olivier.vitrac@gmail.com"
badge = "https://img.shields.io/badge/GitHub-SFPPy-4CAF50?style=for-the-badge&logo=github"


# %% static HTML functions

# SFPPy dynamic version number
def get_version():
    """Extract the version number of SFPPy from VERSION.txt."""
    version_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "VERSION.txt"))
    if not  os.path.exists(version_file):
        raise FileExistsError(f"Error: {version_file} not found. Please create VERSION.txt with content: version=\"X.Y.Z\"\n")
    with open(version_file, "r") as f:
        for line in f:
            match = re.match(r'^version\s*=\s*"(.*?)"$', line.strip())
            if match:
                return match.group(1)
    raise ValueError(f"Error version keyword missing in {version_file}") 

# alert
def create_alert(text=None, fontsize=12, color="#FF4D4D"):
    """return a HTML alert"""
    if not text:
        text = "Do not forget to press all green buttons and refresh interfaces with <kbd>Ctrl+enter</kbd>"
    alert = f"""
<div style="border-left: 4px solid {color}; padding: 10px; background: transparent; color: {color}; 
            font-weight: bold; font-size: {fontsize}px; text-align: left;">
    ⚠️ {text}.
</div>
"""
    return HTML(alert)

# subtitle
def create_subtitle(text=None, fontsize=20, color="#4CAF50"):
    """returns a HTML subtitle"""
    if not text:
        text = "Python Framework for Food Contact Compliance and Risk Assessment 🍏⏩🍎"
    subtitle = f"""
<div style="border-left: 4px solid {color}; padding: 10px; background: transparent; 
            color: {color}; font-weight: bold; font-size: 18px; text-align: left;">
    <span style="font-size: {fontsize}px;">{text}</span> 
</div>
"""
    return HTML(subtitle)

# SFPPy logo
def create_logo():
    """returns SFPPy logo in HTML"""
    version = get_version()
    logo = f"""
<div style="display: flex; justify-content: space-between; align-items: flex-end; font-family: monospace; 
            white-space: pre-wrap; overflow: hidden; font-size: 14px; line-height: 1.3; margin-left: 1cm; max-width: 100%;">
    
<!-- Left: Emoji Block -->
<div style="text-align: left;">
🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️
🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍎🍎🍎🍎🍎🍽️🍽️🍏🍏🍏🍏🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️
🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍏🍽️🍽️🍎🍽️🍽️🍽️🍎🍽️🍽️🐍🍽️🍽️🍽️🐍🍽️
🍽️🍽️🍎🍎🍎🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍏🍏🍏🍏🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍽️🐍🍽️🐍🍽️🍽️
🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🐍🍽️🍽️🍽️
🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🐍🍽️🍽️🍽️
🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️
</div>

<!-- Center: GitHub Badge -->
<div style="margin-left: 15px;">
    <a href="{repo}" target="_blank">
        <img src="{badge}" 
             alt="GitHub SFPPy" style="border-radius: 8px;">
    </a>
</div>

<!-- Right: Version & Email (aligned at the bottom) -->
<div style="margin-left: auto; font-weight: bold; font-size: 14px; color: #FF4D4D; text-align: left;">
    <div style="margin-bottom: 5px;">v{version}</div>
    <a href="mailto:{email}" title="E-mail the author">📩</a>
</div>

</div>
</div>
"""
    return HTML(logo)

# synopsis
def create_synopsis(text=None,color="#4CAF50"):
    """returns a HTML synopsis"""
    if not text:
        text = """
        This template illustrates how to evaluate the migration of substances from a polymeric sleeve into a packaged food simulant using 
        SFPPy (<em>Safety of Food Packaging in Python</em>). Automating key tasks—such as retrieving chemical properties, 
        specifying package geometries, applying polymer parameters, and running mass transfer models ensures transparency
        and reproducibility in compliance testing.
        """
    synopsis = f"""
<div style="padding: 10px; border-left: 4px solid {color};">
    <h3 style="color: {color}; margin-top: 0;">Synopsis</h3>
    <p style="margin-bottom: 1em; color: {color};">
{text}
    </p>
</div>
"""
    return HTML(synopsis)

# disclaimer
def create_disclaimer(fontsize=12):
    """returns a HTML disclaimer"""
    disclaimer = f"""
<div style="border-left: 4px solid #FF4D4D; padding: 10px; background: transparent; color: #FF4D4D; 
            font-weight: bold; font-size: {fontsize}px; text-align: left;">
    ⚠️ DISCLAIMER: This material is provided “<b>AS IS</b>” solely for demonstration and training purposes. No warranty, express or implied, is given regarding its accuracy, completeness, or fitness for a particular use. 📌 Users are solely responsible for evaluating its suitability and for ensuring compliance with all applicable regulations. 🔬 The illustrative example highlights the risks of misinterpreting mass transfer phenomena when "migration modeling" is treated as a "black box". 🚫 Neither the authors nor their organizations accept any liability arising from reliance on or use of this material.
</div>
"""
    return HTML(disclaimer)


# create header with version, footer and separator for notebooks
def create_header_footer(what="head", title="SFPPy - Notebook Index 📑",height=4):
    """
    Create an HTML header or footer block for SFPPy notebooks.

    This function generates a styled HTML block to be used as either a header
    or a footer in SFPPy-related notebooks. The header includes the notebook
    title, a GitHub badge linking to the repository, and version and contact
    information. The footer provides a tagline, contact details, and links to
    the SFPPy website and documentation.

    Parameters:
        what (str): Specifies which block to generate.
                    - If it starts with "head" (e.g., "head" or "header"), the function returns the header.
                    - If it starts with "foot" (e.g., "foot" or "footer"), the function returns the footer.
                    - If it starts with "both", header and footer are returned as a tuple (header,footer)
                    - If it starts with "all", a line separator is added as (header,footer,separator)
        title (str): The title text to display in the header. Defaults to "SFPPy - Notebook Index 📑".
        height (int,float): the <hr> height (default=4)

    Returns:
        IPython.display.HTML: An HTML object containing the header or footer design.
        (header,footer),  (header,footer,separator)

    Raises:
        ValueError: If the 'what' parameter does not start with "head" or "foot".
    """

    version = get_version()

    header =  f"""
  <div style="border-radius: 8px; padding: 12px; background: linear-gradient(to right, #4CAF50, #FF4D4D);
              color: white; font-size: 28px; font-weight: bold; display: flex; align-items: center; justify-content: center; position: relative;">
  {title}
  <a href="{repo}" target="_blank" 
    style="position: absolute; right: 12px; top: 10%; transform: translateY(-10%);">
      <img src="{badge}" 
          alt="GitHub SFPPy" style="border-radius: 8px;">
  </a>
  <div style="position: absolute; right: 48px; top: 82%; transform: translateY(-82%); font-size: 14px; font-weight: bold;">
      <span style="color: white;">v{version}</span> 
      <a href="mailto:{email}" title="E-mail the author" style="margin-left: 8px;">📩</a>
      </div>
  </div>
"""

    footer = f"""
<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:{email}" style="color: #fff; text-decoration: underline;">{author}</a> for questions |
  <a href="{repo}" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="{web}" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>
"""
    
    separator = f"""
<hr style="height: {height}px; background-color: #4CAF50; box-shadow: 2px 2px 4px gray; border: none;">
"""

    if what.startswith("head"):
        return HTML(header)
    if what.startswith("foot"):
        return HTML(footer)
    if what.startswith("both"):
        return (HTML(header),HTML(footer))
    if what.startswith("all"):
        return (HTML(header),HTML(footer),HTML(separator))
    raise ValueError(f'{what} is not recognized, use "head", "foot" , both' or 'all')



# %% Widgets
# create a dropdown widget for files and their execution
def create_files_widget(root="/content/SFPPy/", 
                        folder="notebook", 
                        pattern="*.ipynb", 
                        excluded="index*", 
                        actions=["linkcolab", "linklocal", "run"]):
    """
    Create a dropdown widget with file names and one launch button per specified action.
    
    Parameters:
      - root (str): Full installation path.
      - folder (str or list of str): Folder name(s) (relative to root) to search in.
      - pattern (str or list of str): File pattern(s) to include.
      - excluded (str or list of str): Pattern(s) of files to exclude.
      - actions (list of str): List of actions to create buttons for. Supported actions:
             "linkcolab": Button creates an HTML link to open the file on Colab.
             "linklocal": Button creates an HTML link to open the file locally.
             "run":       Button runs a script file (only applicable for .py files).
    
    Returns:
      A tuple: (dropdown_widget, buttons_dict, output_widget)
        - dropdown_widget: an ipywidgets.Dropdown with the list of found files.
        - buttons_dict: a dictionary mapping each action (as key) to its button widget.
        - output_widget: an ipywidgets.Output widget used to capture output (only if "run" is specified), otherwise None.
        
    The user is expected to display these widgets (e.g. via display()).

    Example usage:
        dropdown_widget, btns, output_widget = create_files_widget(
            root="/content/SFPPy/",
            folder="notebook",      # or a list, e.g. ["notebook", "examples"]
            pattern="*.ipynb",       # or list, e.g. ["*.ipynb", "*.py"]
            excluded="index*",       # or list, e.g. ["index*"]
            actions=["linkcolab", "linklocal", "run"]
        )
        display(dropdown_widget)
        for btn in btns.values():
            display(btn)
        if output_widget:
            display(output_widget)
            
    """
    
    # Ensure folder, pattern, and excluded are lists.
    if not isinstance(folder, list):
        folder = [folder]
    if not isinstance(pattern, list):
        pattern = [pattern]
    if not isinstance(excluded, list):
        excluded = [excluded]
    
    # Search for files in each folder.
    file_list = []
    for fld in folder:
        search_path = os.path.join(root, fld)
        if os.path.exists(search_path):
            for f in os.listdir(search_path):
                # Check if f matches any pattern and does not match any excluded pattern.
                if any(fnmatch.fnmatch(f, pat) for pat in pattern) and not any(fnmatch.fnmatch(f, ex) for ex in excluded):
                    # Save the relative path (i.e., folder/file)
                    file_list.append(os.path.join(fld, f))
        else:
            # If the folder does not exist, assume files are directly under root.
            for f in os.listdir(root):
                if any(fnmatch.fnmatch(f, pat) for pat in pattern) and not any(fnmatch.fnmatch(f, ex) for ex in excluded):
                    file_list.append(f)
    
    # Sort the list in ascending order.
    file_list = sorted(file_list)
    
    # Create the dropdown widget.
    dropdown = widgets.Dropdown(
        options=file_list,
        description="Files:",
    )
    
    # Create an output widget for "run" action only.
    run_out = widgets.Output() if "run" in actions else None
    
    # Create an HTML widget to display links if any link action is specified.
    link_out = widgets.HTML(value="") if any(a in actions for a in ["linkcolab", "linklocal"]) else None
    
    # Dictionary to hold buttons for each action.
    buttons = {}
    
    # Action: "linkcolab"
    if "linkcolab" in actions:
        btn_colab = widgets.Button(description="Open on Colab")
        def on_click_colab(b):
            selected_file = dropdown.value
            # Construct the Colab URL (assumes GitHub hosting on the main branch).
            colab_url = f"https://colab.research.google.com/github/ovitrac/SFPPy/blob/main/{selected_file}"
            # Instead of using display(), update the HTML widget.
            if link_out is not None:
                link_out.value = f'<a href="{colab_url}" target="_blank">Click here to open {selected_file} on Colab</a>'
            else:
                display(HTML(f'<a href="{colab_url}" target="_blank">Click here to open {selected_file} on Colab</a>'))
        btn_colab.on_click(on_click_colab)
        buttons["linkcolab"] = btn_colab
    
    # Action: "linklocal"
    if "linklocal" in actions:
        btn_local = widgets.Button(description="Open Locally")
        def on_click_local(b):
            selected_file = dropdown.value
            # For local opening, use a relative link.
            if link_out is not None:
                link_out.value = f'<a href="{selected_file}" target="_blank">Click here to open {selected_file} locally</a>'
            else:
                display(HTML(f'<a href="{selected_file}" target="_blank">Click here to open {selected_file} locally</a>'))
        btn_local.on_click(on_click_local)
        buttons["linklocal"] = btn_local
        
    # Action: "run"
    if "run" in actions:
        btn_run = widgets.Button(description="Run Script")
        def on_click_run(b):
            selected_file = dropdown.value
            if selected_file.endswith(".py"):
                run_out.clear_output()  # Clear previous output.
                with run_out:
                    print(f"Running {selected_file} ...\n")
                    ip = get_ipython()
                    # Execute the file using the %run magic.
                    ip.run_line_magic("run", selected_file)
            else:
                with run_out:
                    print("The selected file is not a Python script (.py).")
        btn_run.on_click(on_click_run)
        buttons["run"] = btn_run
        
    return dropdown, buttons, run_out, link_out

if __name__ == '__main__':
    from utils.nbutils import create_files_widget
    nbdropdown_widget, nbbtns,_ = create_files_widget(root='~/natacha/python/',
                                                     folder="notebooks", 
                                                     pattern="*.ipynb", 
                                                     excluded="index*", 
                                                     actions=["linkcolab", "linklocal"])
    display(nbdropdown_widget)
    for btn in nbbtns.values(): display(btn)