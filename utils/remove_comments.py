#!/usr/bin/env python3
"""
Script to remove comments and docstrings from a Python file while keeping shebang and UTF-8 encoding declarations.

Usage:
    python remove_comments.py filename.py [suffix] [destination_folder]

Parameters:
    filename.py:       The Python file to process (mandatory).
    suffix:            Optional suffix for the output file (default: .nocomment).
    destination_folder: Optional folder for saving the processed file (default: same as input file).
"""

import sys
import os
import ast
import tokenize
from io import StringIO
import re

# Function to check if a file is Python code
def is_python_file(file_path):
    """Check if the given file is a Python file."""
    _, ext = os.path.splitext(file_path)
    return ext == ".py"

def remove_docstrings(source):
    """
    Remove docstrings from the source code using the ast module.

    Args:
        source (str): The original Python source code.

    Returns:
        str: The source code without docstrings.
    """
    try:
        # Parse the source code into an AST
        parsed_ast = ast.parse(source)
    except SyntaxError as e:
        print(f"Syntax error while parsing the file: {e}")
        sys.exit(1)

    # Define a NodeTransformer to remove docstrings
    class DocstringRemover(ast.NodeTransformer):
        def visit_FunctionDef(self, node):
            self.generic_visit(node)
            if ast.get_docstring(node):
                node.body = node.body[1:]
            return node

        def visit_ClassDef(self, node):
            self.generic_visit(node)
            if ast.get_docstring(node):
                node.body = node.body[1:]
            return node

        def visit_Module(self, node):
            self.generic_visit(node)
            if ast.get_docstring(node):
                node.body = node.body[1:]
            return node

    # Remove docstrings
    remover = DocstringRemover()
    cleaned_ast = remover.visit(parsed_ast)
    ast.fix_missing_locations(cleaned_ast)

    # Convert the AST back to source code
    try:
        cleaned_source = ast.unparse(cleaned_ast)
    except AttributeError:
        # For Python versions < 3.9 where ast.unparse is not available
        try:
            import astor
            cleaned_source = astor.to_source(cleaned_ast)
        except ImportError:
            print("For Python versions below 3.9, please install the 'astor' library: pip install astor")
            sys.exit(1)

    return cleaned_source

def remove_comments(source):
    """
    Remove comments from the source code using the tokenize module.

    Args:
        source (str): The Python source code.

    Returns:
        str: The source code without comments.
    """
    result = []
    g = tokenize.generate_tokens(StringIO(source).readline)
    try:
        for toknum, tokval, start, end, line in g:
            if toknum == tokenize.COMMENT:
                continue  # Skip comments
            elif toknum == tokenize.NL:
                # Preserve standalone newlines
                result.append((toknum, tokval))
            elif toknum == tokenize.NEWLINE:
                result.append((toknum, tokval))
            else:
                result.append((toknum, tokval))
    except tokenize.TokenError as e:
        print(f"Tokenization error: {e}")
        sys.exit(1)
    cleaned_code = tokenize.untokenize(result)
    return cleaned_code

def preserve_preserved_lines(source):
    """
    Preserve shebang and encoding declarations.

    Args:
        source (str): The original Python source code.

    Returns:
        tuple: (preserved_lines, remaining_source)
    """
    preserved_lines = []
    remaining_source = source
    lines = source.splitlines(keepends=True)
    idx = 0
    encoding_re = re.compile(r'coding[:=]\s*([-\w.]+)')

    # Preserve shebang and encoding declarations at the top of the file
    while idx < len(lines):
        line = lines[idx]
        if line.startswith('#!'):
            preserved_lines.append(line)
            idx += 1
        elif line.startswith('#') and encoding_re.search(line):
            preserved_lines.append(line)
            idx += 1
        elif line.strip() == '':
            preserved_lines.append(line)
            idx += 1
        else:
            break

    remaining_source = ''.join(lines[idx:])
    return preserved_lines, remaining_source

def minimize_blank_lines(source):
    """
    Collapse multiple consecutive blank lines into a single blank line.

    Args:
        source (str): The Python source code.

    Returns:
        str: The source code with minimized blank lines.
    """
    # Replace multiple blank lines with two newlines
    source = re.sub(r'\n\s*\n+', '\n\n', source)
    # Optionally, strip trailing whitespace on each line
    source = '\n'.join(line.rstrip() for line in source.splitlines())
    return source.strip() + '\n'  # Ensure the file ends with a single newline

def remove_comments_and_docstrings(source):
    """
    Remove both docstrings and comments from the source code.

    Args:
        source (str): The original Python source code.

    Returns:
        str: The cleaned source code.
    """
    preserved_lines, remaining_source = preserve_preserved_lines(source)
    # Remove docstrings
    no_docstrings = remove_docstrings(remaining_source)
    # Remove comments
    no_comments = remove_comments(no_docstrings)
    # Minimize blank lines
    cleaned_code = minimize_blank_lines(no_comments)
    # Combine preserved lines with cleaned code
    final_content = ''.join(preserved_lines) + cleaned_code
    return final_content

# Main function
def main():
    if len(sys.argv) < 2:
        print("Usage: python remove_comments.py filename.py [suffix] [destination_folder]")
        sys.exit(1)

    input_file = sys.argv[1]
    suffix = sys.argv[2] if len(sys.argv) > 2 else ".nocomment"
    destination_folder = sys.argv[3] if len(sys.argv) > 3 else None

    if not os.path.isfile(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    if not is_python_file(input_file):
        print(f"Error: File '{input_file}' is not a Python file.")
        sys.exit(1)

    try:
        with open(input_file, "r", encoding="utf-8") as f:
            file_content = f.read()

        # Remove comments and docstrings
        cleaned_content = remove_comments_and_docstrings(file_content)

        # Generate output filename
        base, ext = os.path.splitext(os.path.basename(input_file))
        output_file = base + suffix + ext

        if destination_folder:
            # If the destination folder is relative, make it relative to the input file
            if not os.path.isabs(destination_folder):
                destination_folder = os.path.join(os.path.dirname(input_file), destination_folder)
            os.makedirs(destination_folder, exist_ok=True)
            output_file = os.path.join(destination_folder, output_file)
        else:
            output_file = os.path.join(os.path.dirname(input_file), output_file)

        # Write the cleaned content to the output file
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(cleaned_content)

        print(f"Processed file saved to: {output_file}")

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
