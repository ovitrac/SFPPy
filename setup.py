"""
setup.py — Legacy SFPPy Installation Script
============================================

DEPRECATION NOTICE:
    This setup.py is provided for backward compatibility only.
    For new installations, please use pyproject.toml instead:

        pip install .              # Modern installation
        pip install -e .           # Editable/development mode

    This file will be removed in a future version (SFPPy >= 2.0).

@project: SFPPy — Safe Food Packaging in Python
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@agroparistech.fr
@license: MIT
"""

from setuptools import setup, find_packages
import sys
import os
import warnings

# Issue deprecation warning when setup.py is invoked directly
if __name__ == "__main__" or "setup.py" in sys.argv[0]:
    warnings.warn(
        "\n"
        "=" * 70 + "\n"
        "DEPRECATION WARNING: setup.py is deprecated.\n"
        "\n"
        "SFPPy now uses pyproject.toml for package configuration.\n"
        "This setup.py is maintained for backward compatibility only\n"
        "and will be removed in SFPPy >= 2.0.\n"
        "\n"
        "For installation, use:\n"
        "    pip install .         # Standard installation\n"
        "    pip install -e .      # Development/editable mode\n"
        "\n"
        "=" * 70,
        DeprecationWarning,
        stacklevel=2
    )

# Ensure SFPPy root is in sys.path
sfppy_root = os.path.dirname(os.path.abspath(__file__))
if sfppy_root not in sys.path:
    sys.path.insert(0, sfppy_root)

setup(
    name="SFPPy",
    version="1.6",
    description="Software Simulating Mass Transfer from Food Packaging",
    author="Olivier Vitrac",
    author_email="olivier.vitrac@agroparistech.fr",
    url="https://github.com/ovitrac/SFPPy",
    packages=find_packages(include=['patankar', 'patankar.*']),
    install_requires=[
        "numpy>=1.21.0", "matplotlib>=3.4.0", "scipy>=1.7.0", "pandas>=1.3.0", "openpyxl>=3.0.10", "pillow>=8.0.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    include_package_data=True,
    zip_safe=True,
    entry_points={  # CLI command for running scripts from anywhere
        "console_scripts": [
            "sfppy=cli:main",
        ],
    },
)
