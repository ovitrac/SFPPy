"""
Bootstrap Function to Initialize SFPPy Environment (Enhanced from initialize_sfppy.py)

This function is designed for use in local notebooks, Google Colab, or any temporary or read-only environment.
It automatically ensures that SFPPy is available, installed, and imported correctly. It supports offline/local
mode and optionally clones the GitHub repository if needed.

**Two Usages**

___ Locally use it as ___
from utils.bootstrap import bootstrap_sfppy
bootstrap_sfppy(verbose=True, clone="ifcolab")

__ Remotely use it as ___
import urllib.request
exec(urllib.request.urlopen(
    "https://raw.githubusercontent.com/ovitrac/SFPPy/main/utils/bootstrap_sfppy.py"
).read())
bootstrap_sfppy(verbose=True, clone="ifcolab")


**INRAE\Olivier Vitrac**  
Email: [olivier.vitrac@agroparistech.fr](mailto:olivier.vitrac@agroparistech.fr)
"""

import os
import sys
import urllib.request
from IPython import get_ipython


def bootstrap_sfppy(verbose=True, clone="ifcolab"):
    """
    Initializes the SFPPy environment in Colab, Jupyter, or local notebooks.

    Parameters:
    -----------
    verbose : bool
        If True, prints status and actions.

    clone : str
        Controls whether to clone SFPPy if not found.
        Accepted values:
        - "always": always clone the repository if not found.
        - "ifcolab": clone only if running in Google Colab.
        - "never"  : do not clone, raise error if not found.

    Returns:
    --------
    bool : True if initialization succeeded, False otherwise.
    """
    IN_COLAB = 'google.colab' in sys.modules

    def is_sfppy_available():
        try:
            import patankar.loadpubchem
            sfppy_folder = os.path.abspath(os.path.join(os.path.dirname(patankar.loadpubchem.__file__), ".."))
            if verbose:
                print("📁 SFPPy installation folder:", sfppy_folder)
            return True
        except ImportError:
            return False

    def run_initialization_script(script_path):
        if verbose:
            print(f"🔧 Running initialization script: {script_path}")
        get_ipython().run_line_magic("run", f"-i {script_path} --no-verbose --no-cd")

    if is_sfppy_available():
        if verbose:
            print("✅ SFPPy is already available in the environment.")
        return True

    if verbose:
        print("⚠️  SFPPy not detected. Attempting initialization...")

    base_dirs = ["./SFPPy/utils", "../utils", "./utils", ".", "../../utils", "../../../utils"]
    candidates = [os.path.abspath(os.path.join(d, "initialize_sfppy.py"))
                  for d in base_dirs if os.path.exists(os.path.join(d, "initialize_sfppy.py"))]

    if candidates:
        run_initialization_script(candidates[0])
        if is_sfppy_available():
            if verbose:
                print("✅ Success")
            return True
        else:
            if verbose:
                print("❌ Failure.\n❓ Is your installation corrupted?")
            return False

    # No local script found — try cloning (if allowed)
    if clone == "never" or (clone == "ifcolab" and not IN_COLAB):
        if verbose:
            print("❌ No initialization script found and cloning is disabled.")
        return False

    if verbose:
        print("📥 Cloning SFPPy repository from GitHub...")
    try:
        if not os.path.exists("SFPPy"):
            os.system("git clone https://github.com/ovitrac/SFPPy.git") # we clone the first time
        else:
            os.system("git -C SFPPy pull") # we update
        init_script = os.path.abspath(os.path.join("SFPPy", "utils", "initialize_sfppy.py"))
        run_initialization_script(init_script)
    except Exception as e:
        if verbose:
            print("❌ Cloning or initialization failed:", e)
        return False

    if is_sfppy_available():
        if verbose:
            print("✅ SFPPy successfully initialized after cloning.")
        return True
    else:
        if verbose:
            print("❌ Initialization failed. patankar.loadpubchem still not importable.")
        return False
