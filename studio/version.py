"""
SFPPy Studio version information.
"""

__version__ = "0.3.0"
__version_info__ = (0, 3, 0)
__author__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@gmail.com"
__organization__ = "INRAE + Generative Simulation"

# Import SFPPy version
try:
    import sys
    from pathlib import Path
    # Add parent to path for patankar import
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from patankar.migration import __version__ as SFPPy_version
except ImportError:
    SFPPy_version = "unknown"


def get_version_info() -> dict:
    """Get complete version information for traceability."""
    import platform
    import socket

    return {
        "studio_version": __version__,
        "sfppy_version": SFPPy_version,
        "python_version": platform.python_version(),
        "machine_name": socket.gethostname(),
        "machine_os": platform.system(),
        "machine_platform": platform.platform(),
    }
