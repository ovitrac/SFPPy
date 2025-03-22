from setuptools import setup, find_packages
import sys
import os

# Ensure SFPPy root is in sys.path
sfppy_root = os.path.dirname(os.path.abspath(__file__))
if sfppy_root not in sys.path:
    sys.path.insert(0, sfppy_root)

setup(
    name="SFPPy",
    version="1.37",
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
