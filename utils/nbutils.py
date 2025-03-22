#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collection of utilities to manage interactive notebooks

    Author: INRAE\\Olivier Vitrac
    Email: olivier.vitrac@agroparistech.fr
    Last Revised: 2025-03-22
"""

# %% Dependencies
import os, sys, re, fnmatch, datetime, zipfile, shutil, subprocess
import ipywidgets as widgets
from IPython.display import display, HTML, Javascript
from IPython import get_ipython
from pathlib import Path
import nbformat

# %% Constants
author = "Olivier Vitrac"
repo = "https://github.com/ovitrac/SFPPy"
web = "https://ovitrac.github.io/SFPPy/"
email = "olivier.vitrac@gmail.com"
badge = "https://img.shields.io/badge/GitHub-SFPPy-4CAF50?style=for-the-badge&logo=github"
sfppy_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..')) # SFPPy folder

# %% global configurators and exporters

def search_notebook_candidates(filename, search_root=".", max_depth=5):
    """
    Search for matching notebook files downward from `search_root` up to `max_depth`.

    Parameters:
    -----------
    filename : str
        Base name of the notebook (with or without .ipynb).
    search_root : str or Path
        Starting directory for search.
    max_depth : int
        Maximum depth for recursive search.

    Returns:
    --------
    list of Path: All matching .ipynb files with matching base name
    """
    from pathlib import Path
    filename = Path(filename).stem  # strip extension if given
    root = Path(search_root).resolve()
    matches = []

    def _search(path, depth):
        if depth > max_depth:
            return
        for entry in path.iterdir():
            if entry.is_dir():
                _search(entry, depth + 1)
            elif entry.is_file() and entry.name.endswith(".ipynb") and entry.stem == filename:
                matches.append(entry.resolve())

    _search(root, 0)
    return matches


# robust exporter for notebooks (several fallbacks for colab)
def export_notebook(filename=None, add_username=False, fallback_html=True, verbose=True, display_link=True, save_as_zip=True, outputfolder="reports", keep_files=False, colab_download=True):
    """
    Export the current notebook as PDF (Jupyter) or HTML (Colab fallback),
    and also include the original .ipynb and HTML files.
    Adds a timestamp and optional username@host to the file name.
    Optionally creates a .zip archive in the specified output folder.

    If SVG output is used in matplotlib, only HTML export is allowed (PDF would fail).

    Parameters:
    -----------
    filename : str or Path (default: environment variable JPY_SESSION_NAME)
        Name of the notebook to export (with or without .ipynb extension).
        If full path is given, it determines the destination folder.
        If not provided and cannot be inferred, an error is raised.

    add_username : bool (default: False)
        Include the system username@hostname in the exported filename.

    fallback_html : bool (default: True)
        In Colab, fallback to HTML export instead of PDF.

    verbose : bool (default: True)
        Print status messages.

    display_link : bool (default: True)
        Show a download/open link after export (in Colab or Jupyter).

    save_as_zip : bool (default: True)
        If True, create a .zip archive of the exported files.

    outputfolder : str or Path (default: "reports")
        Folder to store the export if filename is not a full path.

    keep_files : bool (default: False)
        If False, remove .ipynb and .html/.pdf after zipping. Ignored if save_as_zip is False.

    colab_download bool (default: True)
        If True, the downloader widget of colab is triggered to download the report file.

    Returns:
    --------
    str : Path to the exported file (PDF or HTML, or ZIP if zipped)
    """
    import matplotlib as mpl

    IN_COLAB = 'google.colab' in sys.modules

    if filename is None:
        filename = os.getenv("JPY_SESSION_NAME")
        if not filename:
            raise RuntimeError("❌ Cannot determine notebook name. Please provide it using the 'filename' argument.")

    path = Path(filename).expanduser().resolve() if Path(filename).is_absolute() else Path(os.getcwd()) / filename
    notebook_path = path if path.suffix == ".ipynb" else path.with_suffix(".ipynb")

    if not notebook_path.exists():
        candidates = search_notebook_candidates(filename, search_root=".", max_depth=3)
        if len(candidates) == 1:
            notebook_path = candidates[0]
            if verbose:
                print(f"✅ Notebook found automatically: {notebook_path}")
        elif len(candidates) > 1:
            print("❌ Multiple candidate notebooks found:")
            for i, match in enumerate(candidates, 1):
                print(f"  {i}. {match}")
            print("👉 Please change to the correct folder using:")
            print(f"   %cd {candidates[0].parent}")
            raise FileNotFoundError("❌ Multiple matching notebooks found. Please disambiguate manually.")
        else:
            raise FileNotFoundError(f"❌ Could not find notebook: {notebook_path}. Please provide the correct filename.")


    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    user = os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"
    host = os.uname().nodename if hasattr(os, 'uname') else os.environ.get("COMPUTERNAME", "unknown-host")
    userhost = f"{user}@{host}" if add_username else ""
    suffix = f"{timestamp}_{userhost}".strip("_")

    using_svg = mpl.rcParams.get("figure.format", None) == "svg"
    if using_svg:
        if verbose:
            print("⚠️ Matplotlib is using SVG output. PDF export will be disabled.")
            print("👉 Please add `from IPython.display import set_matplotlib_formats; set_matplotlib_formats('retina')` to enable PDF export.")
        fallback_html = True

    output_ext = "html" if IN_COLAB or fallback_html else "pdf"
    base_filename = f"{notebook_path.stem}_{suffix}"
    output_dir = notebook_path.parent if Path(filename).is_absolute() else Path(outputfolder).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{base_filename}.{output_ext}"
    output_ipynb_path = output_dir / f"{base_filename}.ipynb"
    zip_path = output_dir / f"{base_filename}.zip"

    command = [
        "jupyter", "nbconvert",
        "--to", output_ext,
        "--output", output_path.name,
        "--output-dir", str(output_dir),
        str(notebook_path.name)
    ]

    try:
        subprocess.run(command, check=True, cwd=notebook_path.parent, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        raise RuntimeError(f"❌ Failed to export notebook to {output_ext.upper()}") from e

    try:
        shutil.copy(notebook_path, output_ipynb_path)
        if verbose:
            print(f"📚 Copied .ipynb as: {output_ipynb_path}")
    except Exception as e:
        raise RuntimeError("❌ Failed to copy .ipynb file.") from e

    if verbose:
        print(f"📤 Exported to: {output_path}")

    if save_as_zip:
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(output_path, arcname=output_path.name)
            zipf.write(output_ipynb_path, arcname=output_ipynb_path.name)
        if not keep_files:
            os.remove(output_path)
            os.remove(output_ipynb_path)
        if verbose:
            print(f"📦 Zipped export to: {zip_path}")

    if display_link:
        try:
            target_path = zip_path if save_as_zip else output_path
            if IN_COLAB:
                from google.colab import files
                files.download(str(target_path))  # trigger download
            else:
                rel_path = os.path.relpath(target_path, os.getcwd())
                label = f"📄 Download {target_path.name}"
                display(HTML(f"<a href='{rel_path}' target='_blank'>{label}</a>"))
        except Exception as e:
            if verbose:
                print(f"⚠️ Could not generate display link.\n Reason: {e}")

    return str(zip_path if save_as_zip else output_path)




def set_figure_format(fmt="png", verbose=True):
    """
    Set the inline figure rendering format for Jupyter or Colab notebooks.

    Parameters
    ----------
    fmt : str
        One of: 'png', 'jpg', 'jpeg', 'svg', 'pdf', 'retina', 'png2x'
        - 'retina' and 'png2x' are aliases for high-DPI PNG (2x resolution)
        - 'pdf' and 'svg' are vector formats (great for publication)
        - 'png' is the default (fastest and most compatible)

    verbose : bool
        If True, prints a confirmation message.

    Raises
    ------
    ValueError
        If an unsupported format is provided.
    """
    VALID_FORMATS = ['jpg', 'jpeg', 'png', 'svg', 'pdf', 'retina', 'png2x']
    fmt = fmt.lower()
    if fmt not in VALID_FORMATS:
        raise ValueError(f"Unsupported figure format: '{fmt}'. Must be one of {VALID_FORMATS}.")

    ipy = get_ipython()
    if ipy is None:
        raise RuntimeError("This function must be run inside a Jupyter or Colab notebook.")

    # Use inline backend with appropriate format
    ipy.run_line_magic("matplotlib", "inline")

    # Retina and png2x are aliases for high-DPI PNG
    if fmt in ["retina", "png2x"]:
        fmt = "retina"
    elif fmt in ["jpg", "jpeg"]:
        # InlineBackend does not support jpg; fallback to PNG and show a warning
        if verbose:
            print("⚠️ Inline display of 'jpg/jpeg' is not supported. Falling back to PNG.")
        fmt = "png"

    ipy.run_line_magic("config", f"InlineBackend.figure_format = '{fmt}'")
    if verbose:
        print(f"📊 Matplotlib inline figure format set to: {fmt}")


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
    logo_old = f"""
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
"""
    logo = f"""
<hr style="height: 4px; background-color: #4CAF50; box-shadow: 2px 2px 4px gray; border: none;">
<div style="display: flex; align-items: flex-start; font-family: monospace; font-size: 14px; line-height: 1.3; margin: 0.5em 0; max-width: 100%;">
  <!-- Left: Emoji “Logo” (Single line, scrollable if too wide) -->
  <div style="white-space: nowrap; overflow-x: auto; margin-right: 15px;">
    🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️<br>
    🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍎🍎🍎🍎🍎🍽️🍽️🍏🍏🍏🍏🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️<br>
    🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍏🍽️🍽️🍎🍽️🍽️🍽️🍎🍽️🍽️🐍🍽️🍽️🍽️🐍🍽️<br>
    🍽️🍽️🍎🍎🍎🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍏🍏🍏🍏🍽️🍽️🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍽️🐍🍽️🐍🍽️🍽️<br>
    🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🐍🍽️🍽️🍽️<br>
    🍽️🍎🍎🍎🍎🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍏🍽️🍽️🍽️🍽️🍽️🍽️🍎🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🐍🍽️🍽️🍽️<br>
    🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️🍽️<br>
  </div>
  <!-- Center: GitHub Badge + Version (column layout) -->
  <div style="display: flex; flex-direction: column; align-items: center; margin-right: 15px;">
    <a href="{repo}" target="_blank">
      <img src="{badge}" alt="GitHub SFPPy" style="border-radius: 8px;">
    </a>
    <div style="margin-top: 4px; font-weight: bold; color: #FF4D4D;">v{version} <a href="mailto:{email}" title="E-mail the author">📩</a></div>
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
    disclaimer_old = f"""
<div style="border-left: 4px solid #FF4D4D; padding: 10px; background: transparent; color: #FF4D4D;
            font-weight: bold; font-size: {fontsize}px; text-align: left;">
    ⚠️ DISCLAIMER: This material is provided “<b>AS IS</b>” solely for demonstration and training purposes. No warranty, express or implied, is given regarding its accuracy, completeness, or fitness for a particular use. 📌 Users are solely responsible for evaluating its suitability and for ensuring compliance with all applicable regulations. 🔬 The illustrative example highlights the risks of misinterpreting mass transfer phenomena when "migration modeling" is treated as a "black box". 🚫 Neither the authors nor their organizations accept any liability arising from reliance on or use of this material.
</div>
"""
    disclaimer = f"""
<div style="border-left: 4px solid #FF4D4D; padding: 10px; background: transparent; color: #FF4D4D;
        font-weight: bold; font-size: {fontsize}px; text-align: left;">
⚠️ <strong>DISCLAIMER</strong><br>
This material is provided “AS IS” solely for demonstration and training. No warranty—express or implied—is given regarding its accuracy, completeness, or suitability for any particular purpose. 📌 Users are fully responsible for assessing its relevance and ensuring compliance with all applicable regulations. 🔬 The illustrative example underscores the risks of treating “migration modeling” as a mere “black box,” potentially leading to misinterpretation of mass transfer phenomena. 🚫 Neither the authors nor their organizations accept any liability for reliance on or use of this material.
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
def create_files_widget(root=sfppy_folder,
                        folder="",
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


# %% Notebook explorer widget and dependencies
def clean_markdown(text):
    """
    Clean a text string by removing extraneous Markdown formatting markers.

    This function removes Markdown heading markers (e.g. "#", "##", etc.) that
    typically appear at the beginning of lines and also removes emphasis markers
    (single or double asterisks) used for italic or bold formatting. Multiplication
    operators (i.e. asterisks surrounded by digits) are preserved. Additionally,
    consecutive empty lines are collapsed into a single empty line.

    Parameters:
        text (str): The input text string, which may include Markdown formatting.

    Returns:
        str: The cleaned text with unnecessary Markdown markers removed.
    # Example usage:
    md_text = '''
    ## Heading Example

    *Italic text* and **bold text** should be cleaned,
    but multiplication 3 * 4 must remain unchanged.

    Another paragraph.


    Extra blank lines should be collapsed.
    '''
    cleaned = clean_markdown(md_text)
    print(cleaned)
    """
    # Remove Markdown headings at the start of lines (e.g. "# ", "## ", etc.)
    text = re.sub(r'^\s*#+\s+', '', text, flags=re.MULTILINE)
    # Remove Markdown bold formatting: replace **text** with text.
    text = re.sub(r'\*\*(.*?)\*\*', r'\1', text)
    # Remove Markdown italic formatting: replace *text* with text.
    text = re.sub(r'\*(.*?)\*', r'\1', text)
    # Remove any remaining asterisks that are not part of a multiplication expression.
    # A multiplication expression is assumed to have a digit before and after the asterisk.
    text = re.sub(r'(?<!\d)\*+(?!\d)', '', text)
    # Collapse multiple empty lines to no more than one empty line.
    text = re.sub(r'\n\s*\n+', '\n\n', text)
    return text



# copy to clipboard
def copy_to_clipboard(text):
    """copy text to clipboard, but very challenging"""
    js = f"""
    (function() {{
      // If the secure clipboard API is available, use it.
      if (navigator.clipboard && window.isSecureContext) {{
          navigator.clipboard.writeText({repr(text)}).then(function() {{
              console.log('Copying text was successful.');
          }}, function(err) {{
              console.error('Failed to copy text: ', err);
          }});
      }} else {{
          // Fallback to execCommand method
          var textArea = document.createElement("textarea");
          textArea.value = {repr(text)};
          // Avoid scrolling to bottom
          textArea.style.position = "fixed";
          textArea.style.top = "-9999px";
          document.body.appendChild(textArea);
          textArea.focus();
          textArea.select();
          try {{
              var successful = document.execCommand('copy');
              console.log('Copying text command was ' + (successful ? 'successful' : 'unsuccessful'));
          }} catch (err) {{
              console.error('Unable to copy', err);
          }}
          document.body.removeChild(textArea);
      }}
    }})();
    """
    print(text)
    display(Javascript(js))

# nb code extraction
def extract_code_segments(nb_path):
    """
    Extract code segments and their associated comments from a Jupyter Notebook.

    Reads the notebook at nb_path and returns a list of tuples (comment, code).
    For each code cell, if the immediately preceding cell is markdown, that cell's
    content is used as the comment; otherwise, the comment is an empty string.

    Parameters:
        nb_path (str): Path to the notebook (.ipynb).

    Returns:
        List[Tuple[str, str]]: A list of (comment, code) pairs.

    Example usage:
    segments = extract_code_segments("path/to/notebook.ipynb")
    for comment, code in segments:
        print("Comment:")
        print(comment)
        print("Code:")
        print(code)
        print("="*40)
    """
    segments = []
    nb = nbformat.read(nb_path, as_version=4)
    cells = nb.cells
    for i, cell in enumerate(cells):
        if cell.cell_type == "code":
            comment = ""
            if i > 0 and cells[i-1].cell_type == "markdown":
                comment = cells[i-1].source.strip()
            code = cell.source.strip()
            segments.append((clean_markdown(comment), code))
    return segments

#nb selector
def create_notebook_explorer(folder=""):
    """
    Create a widget to select a notebook file from the 'notebooks' directory.

    It searches in the directory:
        os.path.join(os.path.dirname(__file__), '..', 'notebooks', folder)
    for all files matching *.ipynb, and returns a dropdown populated with the filenames.
    A button labeled "Open Notebook" is also returned. When pressed, it returns the full
    path to the selected notebook.

    Parameters:
        folder (str): Subfolder within the 'notebooks' directory (default: "").

    Returns:
        tuple: (dropdown_widget, open_button)
            dropdown_widget: ipywidgets.Dropdown containing the notebook filenames.
            open_button: ipywidgets.Button that, when clicked, triggers further actions.
    """
    base_dir = os.path.join(os.path.dirname(__file__), '..', 'notebooks', folder)
    # List all .ipynb files
    notebooks = sorted([f for f in os.listdir(base_dir) if f.endswith('.ipynb')])

    dropdown = widgets.Dropdown(options=notebooks, description="Notebook:")
    open_button = widgets.Button(description="Open Notebook")
    output_area = widgets.Output()

    # For demonstration, we attach a click handler that prints the full path.
    def on_open(b):
        full_path = os.path.abspath(os.path.join(base_dir, dropdown.value))
        print("Selected notebook:", full_path)
        with output_area:
            output_area.clear_output()  # Clear previous output
            navigator = create_notebook_navigator(full_path)
            display(navigator)
    open_button.on_click(on_open)
    display(dropdown, open_button, output_area)



def create_notebook_navigator(nb_path=None):
    """
    Create a navigation widget for a notebook's code segments.

    This widget displays a horizontal navigation bar with buttons:
      [first] [<--] [-->] [last]    "Cell #/N"   [copycode] [copycomment] [close]
    Below it, the segment's associated comment (limited to three visible lines, with
    scrolling if longer) and the full code are shown.

    Users can navigate among segments, copy either the code or the comment to the clipboard,
    and close the navigator. This widget helps overcome Colab's limitation of only opening
    one notebook at a time by letting users view and copy parts of a notebook.

    Parameters:
        nb_path (str): Full path to the notebook (.ipynb) from which to extract segments.

    Returns:
        ipywidgets.Widget: A container widget with the navigation interface.

    Example Usage:
        Notebook Selector widget
        nb_dropdown, open_btn = create_notebook_selector(folder="")
        display(nb_dropdown, open_btn)

        To launch the navigator, you might attach a callback to open_btn that creates and displays the navigator:

        def open_notebook_navigator(b):
            base_dir = os.path.join(os.path.dirname(__file__), '..', 'notebooks', "")
            nb_path = os.path.abspath(os.path.join(base_dir, nb_dropdown.value))
            navigator = create_notebook_navigator(nb_path)
            display(navigator)

        open_btn.on_click(open_notebook_navigator)
    """
    segments = extract_code_segments(nb_path)
    total = len(segments)
    if total == 0:
        return widgets.HTML("<b>No code segments found.</b>")

    current_index = 0  # closure variable for the current segment

    # Navigation buttons
    first_btn = widgets.Button(description="First")
    prev_btn  = widgets.Button(description="<--")
    next_btn  = widgets.Button(description="-->")
    last_btn  = widgets.Button(description="Last")
    pos_label = widgets.Label(value=f"Cell 1/{total}")

    # Copy buttons
    copycode_btn    = widgets.Button(description="Copy Code")
    copycomment_btn = widgets.Button(description="Copy Comment")
    close_btn       = widgets.Button(description="Close")

    # Display areas for comment and code
    comment_area = widgets.Textarea(
        value="",
        description="Comment:",
        layout=widgets.Layout(width="100%", height="80px"),
        disabled=True
    )
    code_area = widgets.Textarea(
        value="",
        description="Code:",
        layout=widgets.Layout(width="100%", height="200px"),
        disabled=True
    )

    # Update function to refresh display based on current_index
    def update_display():
        comment, code = segments[current_index]
        comment_area.value = comment
        code_area.value = code
        pos_label.value = f"Cell {current_index+1}/{total}"

    # Navigation callbacks
    def on_first(b):
        nonlocal current_index
        current_index = 0
        update_display()
    def on_prev(b):
        nonlocal current_index
        if current_index > 0:
            current_index -= 1
            update_display()
    def on_next(b):
        nonlocal current_index
        if current_index < total - 1:
            current_index += 1
            update_display()
    def on_last(b):
        nonlocal current_index
        current_index = total - 1
        update_display()

    first_btn.on_click(on_first)
    prev_btn.on_click(on_prev)
    next_btn.on_click(on_next)
    last_btn.on_click(on_last)

    # Copy callbacks using the helper function
    def on_copycode(b):
        copy_to_clipboard(segments[current_index][1])
    def on_copycomment(b):
        copy_to_clipboard(segments[current_index][0])
    copycode_btn.on_click(on_copycode)
    copycomment_btn.on_click(on_copycomment)

    # Close callback: hides the navigator container.
    def on_close(b):
        container.layout.display = "none"
    close_btn.on_click(on_close)

    # Assemble the top navigation bar
    nav_bar = widgets.HBox([
        first_btn, prev_btn, next_btn, last_btn,
        pos_label,
        copycode_btn, copycomment_btn, close_btn
    ])

    # Create a container for the whole navigator
    container = widgets.VBox([
        nav_bar,
        comment_area,
        widgets.HTML("<hr>"),
        code_area
    ])

    # Initialize display
    update_display()
    return container


# %% for debugging
if __name__ == '__main__':
    nbdropdown_widget, nbbtns,_ = create_files_widget(root='~/natacha/python/',
                                                     folder="notebooks",
                                                     pattern="*.ipynb",
                                                     excluded="index*",
                                                     actions=["linkcolab", "linklocal"])
    display(nbdropdown_widget)
    for btn in nbbtns.values(): display(btn)
