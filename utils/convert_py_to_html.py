#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    convert_py_to_html.py
    =====================
    Script to Convert Python Files to HTML with Syntax Highlighting Using Pygments.

    This script reads a Python file and generates an HTML file with syntax highlighting.
    The output HTML includes:
        - Syntax-highlighted Python code, styled using Pygments.
        - A red-bordered notice at the top indicating that the raw code is displayed,
        along with the date and time of the conversion.
        - Fully embedded CSS for standalone viewing.

    The script is useful as a fallback when automated documentation generation fails,
    allowing users to view raw Python code in a visually formatted HTML page.

    Usage:
        python convert_py_to_html.py <input_file.py> <output_file.html>

    Positional Arguments:
        input_file (str): Path to the input Python file to be converted.
        output_file (str): Path where the resulting HTML file will be saved.

    Requirements:
        - Python 3.x
        - Pygments library (install via `pip install pygments`)

    Behavior:
        - Reads the content of the specified Python file.
        - Generates HTML with syntax highlighting and a notice message.
        - Saves the resulting HTML file to the specified output path.

    Example:
        python convert_py_to_html.py example.py example.html

        This will read `example.py` and create an HTML file named `example.html`
        in the current directory, containing the highlighted Python code.

    Error Handling:
        - If the input file does not exist, the script will display an error message.
        - Any issues during file read/write operations will also be logged.

    Note:
        Ensure that the input file path is correct, and you have write permissions
        to the output directory.

    Author:
        - **INRAE\Olivier Vitrac**
        - **Email:** olivier.vitrac@agroparistech.fr
        - **Last Revised:** 2024-12-26
        
"""

import os
import datetime
from pygments import highlight
from pygments.lexers import PythonLexer
from pygments.formatters import HtmlFormatter

def convert_python_to_html(input_file, output_file):
    """
    Converts a Python file to an HTML file with syntax highlighting using Pygments.
    Adds a frame at the top with a message about raw code display and the conversion date.
    """
    try:
        # Read the Python file
        with open(input_file, 'r') as f:
            code = f.read()
        
        # Generate the HTML with Pygments
        formatter = HtmlFormatter(full=True, noclasses=True)  # Embed CSS
        highlighted_code = highlight(code, PythonLexer(), formatter)
        
        # Add a frame with a message and date
        message = f"""
        <div style="border: 2px solid red; padding: 10px; margin-bottom: 20px;">
            <h3 style="color: red;">Notice</h3>
            <p>This is a raw view of the Python source code due to an error in generating the documentation.</p>
            <p><strong>Date of Conversion:</strong> {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        """
        
        # Combine the frame and the highlighted code
        final_html = highlighted_code.replace("<body>", f"<body>{message}", 1)
        
        # Write the final HTML to the output file
        with open(output_file, 'w') as f:
            f.write(final_html)
        
        print(f"Converted {input_file} to {output_file}")
    
    except Exception as e:
        print(f"Failed to convert {input_file} to HTML: {e}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python convert_py_to_html.py <input_file.py> <output_file.html>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' does not exist.")
        else:
            convert_python_to_html(input_file, output_file)
