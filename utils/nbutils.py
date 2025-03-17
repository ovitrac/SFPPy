"""
Collection of utilities to manage interactive notebooks

    Author: INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2025-03-17
"""

import os
import fnmatch
import ipywidgets as widgets
from IPython.display import display, HTML
from IPython import get_ipython

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