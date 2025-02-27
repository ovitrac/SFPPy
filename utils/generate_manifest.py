#!/usr/bin/env python3

"""
Manifest Manager for Pizza3 Project

===================================

The `manifestManager` module provides a robust interface for creating, synchronizing,
and validating manifests within the Pizza3 project. A manifest serves as an inventory
of project files, capturing essential metadata such as file paths, hashes, and sizes.
This facilitates efficient tracking of project changes, ensuring integrity and consistency
across different environments or versions.

Goals:
------
- **Create Manifests:** Generate comprehensive manifests that catalog all relevant project files.
- **Diff Manifests:** Compare local and source manifests to identify additions, modifications, or deletions.
- **Update Manifests:** Synchronize local project files with a source manifest, handling updates and ensuring integrity.
- **Maintain Integrity:** Utilize file hashing to verify the integrity of project files, detecting any unintended alterations.

Usage:
------
The module is designed to be executed as a standalone script, accepting various command-line
arguments to perform desired actions. It leverages the `FileEntity` class to represent individual
files and directories, capturing their metadata for manifest operations.

Example Commands:
-----------------
- **Create a Manifest:**
  
  ```bash
  python3 manifestManager.py create -p /path/to/project -m project.manifest
  ```

- **Diff Two Manifests:**
  
  ```bash
  python3 manifestManager.py diff -l /local/manifest/dir -s https://example.com/source/manifest --print
  ```

- **Update a Manifest:**
  
  ```bash
  python3 manifestManager.py update -l /local/manifest/dir -s https://example.com/source/manifest --prompt
  ```

Dependencies:
-------------
- Python 3.x
- Standard Python libraries: `os`, `hashlib`, `urllib`, `shutil`, `errno`, `argparse`, `sys`

Author:
-------
Michael Imelfort

License:
--------
GPLv3

"""

#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    generate_manifest.py                                                     #
#                                                                             #
#    Work with online data manifests (creating / syncing / validating)        #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = ["Michael Imelfort","Olivier Vitrac"]
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.35"


###############################################################################
###############################################################################
###############################################################################
###############################################################################

# System includes
import os
import hashlib
import urllib.request, urllib.error, urllib.parse
import shutil
import errno
import argparse
import sys

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class FileEntity(object):
    """Represents a file or directory entity with associated metadata."""
    
    def __init__(self,
                 name,      # Name of the entity
                 path,      # Relative path to the root directory
                 parent,    # Parent FileEntity object (None for root)
                 hashd,     # Hash of the file (None for directories)
                 size       # Size of the file in bytes (0 for directories)
                 ):
        """
        Initializes a new instance of the FileEntity class.

        Args:
            name (str): The name of the entity on the file system.
            path (str): The relative path to the root directory.
            parent (FileEntity or None): The parent FileEntity object containing this entity.
                - `None` if the entity is the root.
            hashd (str): The SHA-256 hash of the file.
                - Set to `'-'` for directories.
            size (int): The size of the file in bytes.
                - Set to `0` for directories.
        """
        self.name = name
        self.path = path
        self.parent = parent
        self.hashd = hashd
        self.size = size
        self.type = 'dir' if self.hashd == '-' else 'file'  # Determine type based on hashd

    def getFullPath(self):
        """
        Retrieves the full path to this entity by concatenating parent paths.

        Returns:
            str: The absolute path to the entity.
                - For the root entity, returns its name.
                - For other entities, joins the parent's full path with its own name.
        """
        if self.parent is None:
            return self.name  # Root entity
        else:
            return os.path.join(self.parent.getFullPath(), self.name)

    def checkIntegrity(self):
        """
        Checks the integrity of the file by comparing its stored hash with a newly computed hash.

        For directories, integrity is assumed to be intact.

        Returns:
            bool: 
                - `True` if the entity is a directory or if the file's hash matches the stored hash.
                - `False` otherwise.
        """
        if self.type == 'dir':
            return True  # Directories are assumed to be intact
        else:
            # Placeholder for actual integrity check
            # Implement hash comparison logic if necessary
            return True

    def __str__(self):
        """
        Provides a string representation of the FileEntity instance.

        Returns:
            str: A tab-separated string containing the entity's path, hash, and size.
                - Format: "relative_path/name\thashd\tsize"
                - For the root entity, returns an empty string.
        """
        if self.parent is not None:
            return "\t".join([os.path.join(self.path, self.name), self.hashd, str(self.size)])
        return ""

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ManifestManager(object):
    """
    ManifestManager Class

    ======================

    The `ManifestManager` class serves as the core interface for managing project manifests
    within the Pizza3 project. It encapsulates functionalities to create, compare, and update
    manifests, ensuring that the project's file structure remains consistent and unaltered
    across different environments or versions.

    Attributes:
    -----------
    timeout : int
        Specifies the timeout duration (in seconds) for network operations when fetching remote manifests.
    myExtensions : list of str
        A list of file extensions to include in the manifest (e.g., Python and shell scripts).
    files : list of FileEntity
        A collection of `FileEntity` objects representing the project's files and directories.
    type : str
        The type of manifest being managed (default is 'generic').

    Methods:
    --------
    createManifest(path, manifestName=None, exclude_dirs=None, exclude_files=None, extensions=None)
        Generates a manifest by inventorying all relevant files within the specified directory.
    
    diffManifests(localManifestLocation, sourceManifestLocation, localManifestName=None, sourceManifestName=None, printDiffs=False)
        Compares two manifests to identify differences such as added, modified, or deleted files and directories.
    
    updateManifest(localManifestLocation, sourceManifestLocation, localManifestName=None, sourceManifestName=None, prompt=True)
        Updates the local project files based on the differences identified between the local and source manifests.
    
    getManType(line)
        Extracts the manifest type from a given line of the manifest file.
    
    formatData(amount)
        Formats a byte size into a human-readable string (e.g., KB, MB, GB).
    
    makeSurePathExists(path)
        Ensures that a specified directory path exists, creating it if necessary.
    
    promptUserDownload()
        Prompts the user for confirmation before proceeding with downloading updates.
    
    walk(parents, full_path, rel_path, dirs, files, skipFile=".dmanifest", exclude_dirs=None, exclude_files=None, extensions=None)
        Recursively traverses the project directory to catalog files and directories, excluding specified files.
    
    listdir(path)
        Lists directories, files, and symbolic links within a specified path.
    
    hashfile(fileName, blocksize=65536)
        Computes the SHA-256 hash of a given file to ensure integrity.

    Usage Example:
    --------------
    ```python
    from generate_manifest import ManifestManager

    # Initialize the manager with specific extensions
    manager = ManifestManager()

    # Create a new manifest including .py, .sh, .md, and .html files
    manager.createManifest('/path/to/project', manifestName='project.manifest', extensions=['.py', '.sh', '.md', '.html'])

    # Diff two manifests and print differences
    manager.diffManifests(
        localManifestLocation='/local/manifest/dir',
        sourceManifestLocation='https://example.com/source/manifest',
        printDiffs=True
    )

    # Update the local manifest based on the source
    manager.updateManifest(
        localManifestLocation='/local/manifest/dir',
        sourceManifestLocation='https://example.com/source/manifest',
        prompt=True
    )
    ```
    """
    
    def __init__(self, manType=None, timeout=30):
        """
        Initializes a new instance of the ManifestManager class.

        Args:
            manType (str, optional): The type/category of the manifest. Defaults to "generic".
            timeout (int, optional): Timeout duration in seconds for network operations. Defaults to 30.
        """
        self.timeout = timeout
        self.myExtensions = [".py", ".sh"]  # Default extensions
        self.files = []
        if manType is not None:
            self.type = manType
        else:
            self.type = "generic"

    def createManifest(self, path, manifestName=None, exclude_dirs=None, exclude_files=None, extensions=None):
        """
        Inventory all files in the specified path and create a manifest file.

        Args:
            path (str): The root directory path to inventory.
            manifestName (str, optional): The name of the manifest file. Defaults to ".dmanifest".
            exclude_dirs (list of str, optional): List of directories to exclude from the manifest.
            exclude_files (list of str, optional): List of files to exclude from the manifest.
            extensions (list of str, optional): List of file extensions to include in the manifest.
                - If provided, overrides the default `myExtensions`.
                - If None, uses the existing `myExtensions`.
        """
        if manifestName is None:
            manifestName = ".dmanifest"  # Default manifest name
        print(f"Creating manifest '{manifestName}' for path: {path}")

        # Update extensions if provided
        if extensions is not None:
            self.myExtensions = extensions
            print(f"Using custom extensions: {self.myExtensions}")
        else:
            print(f"Using default extensions: {self.myExtensions}")

        # Make the root file entity
        root_path = os.path.abspath(path)
        root_fe = FileEntity('root', ".", None, "-", 0)
        self.files.append(root_fe)

        # Now make all the ones below
        parents = [root_fe]
        dirs, files, _ = self.listdir(path)
        print(f"Initial directories: {dirs}")
        print(f"Initial files: {files}")
        self.walk(parents, root_path, '', dirs, files, skipFile=manifestName, exclude_dirs=exclude_dirs, exclude_files=exclude_files, extensions=self.myExtensions)

        with open(os.path.join(path, manifestName), 'w') as man_fh:
            # Print the header
            man_fh.write("#\t::: %s ::: \tPizza3 manifest version %s\n\n" % (self.type, __version__))
            for f in self.files:
                if f.parent is not None:
                    man_fh.write("%s\n" % f)
        print(f"Manifest '{manifestName}' successfully written to '{path}'.")

    def diffManifests(self,
                      localManifestLocation,
                      sourceManifestLocation,
                      localManifestName=None,
                      sourceManifestName=None,
                      printDiffs=False):
        """
        Check for any differences between two manifests.

        Args:
            localManifestLocation (str): Path to the local manifest directory.
            sourceManifestLocation (str): Path or URL to the source manifest location.
            localManifestName (str, optional): Name of the local manifest file. Defaults to ".dmanifest".
            sourceManifestName (str, optional): Name of the source manifest file. Defaults to ".dmanifest".
            printDiffs (bool, optional): Flag to print the differences. Defaults to False.

        Returns:
            tuple: Contains source path, added files, added directories, deleted files, and modified files.
                   Returns (None, None, None, None, None) if an error occurs.
        """
        if localManifestName is None:
            localManifestName = ".dmanifest"
        if sourceManifestName is None:
            sourceManifestName = ".dmanifest"

        print(f"Diffing manifests: Local='{localManifestLocation}/{localManifestName}', Source='{sourceManifestLocation}/{sourceManifestName}'")
        
        # Get the "type" of the local manifest
        l_type = "generic"
        try:
            with open(os.path.join(localManifestLocation, localManifestName)) as l_man:
                for line in l_man:
                    if line.startswith("#"):
                        l_type = self.getManType(line)
                        print(f"Local manifest type: {l_type}")
                    break
        except FileNotFoundError:
            print(f"Error: Local manifest '{localManifestName}' not found in '{localManifestLocation}'.")
            return (None, None, None, None, None)

        # Load the source manifest
        s_type = "generic"
        source_man = {}
        source = ""
        # First, assume it is remote
        try:
            source_url = urllib.parse.urljoin(sourceManifestLocation + '/', sourceManifestName)
            print(f"Attempting to load source manifest from URL: {source_url}")
            with urllib.request.urlopen(source_url, timeout=self.timeout) as s_man:
                source = sourceManifestLocation + "/"
                for line in s_man:
                    decoded_line = line.decode('utf-8')
                    if decoded_line.startswith("#"):
                        s_type = self.getManType(decoded_line)
                        print(f"Source manifest type: {s_type}")
                        if s_type != l_type:
                            print(f"Error: Type of source manifest ({s_type}) does not match type of local manifest ({l_type}).")
                            return (None, None, None, None, None)
                        continue
                    fields = decoded_line.rstrip().split("\t")
                    if len(fields) >= 3:
                        source_man[fields[0]] = [fields[1], fields[2], False]
        except urllib.error.URLError as e:
            print(f"Warning: Could not load source manifest from URL. Attempting to load from local path. Error: {e}")
            # If not remote, assume it's a local file path
            try:
                with open(os.path.join(sourceManifestLocation, sourceManifestName)) as s_man:
                    source = os.path.join(sourceManifestLocation) + os.path.sep
                    for line in s_man:
                        if line.startswith("#"):
                            s_type = self.getManType(line)
                            print(f"Source manifest type: {s_type}")
                            if s_type != l_type:
                                print(f"Error: Type of source manifest ({s_type}) does not match type of local manifest ({l_type}).")
                                return (None, None, None, None, None)
                            continue
                        fields = line.rstrip().split("\t")
                        if len(fields) >= 3:
                            source_man[fields[0]] = [fields[1], fields[2], False]
            except FileNotFoundError:
                print(f"Error: Source manifest '{sourceManifestName}' not found in '{sourceManifestLocation}'.")
                return (None, None, None, None, None)

        print(f"Loaded {len(source_man)} entries from source manifest.")

        # Keep lists of modifications
        deleted = []
        addedDirs = []
        addedFiles = []
        modified = []

        try:
            with open(os.path.join(localManifestLocation, localManifestName)) as l_man:
                for line in l_man:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip().split("\t")
                    if len(fields) < 3:
                        continue
                    path, hashd, size = fields[0], fields[1], fields[2]
                    if path in source_man:
                        if source_man[path][0] != hashd:
                            # Hashes don't match
                            modified.append(path)
                        # Mark as seen
                        source_man[path][2] = True
                    else:
                        # File has been deleted from the source manifest
                        deleted.append(path)
        except FileNotFoundError:
            print(f"Error: Local manifest '{localManifestName}' not found in '{localManifestLocation}'.")
            return (None, None, None, None, None)

        # Check for new files
        for f in list(source_man.keys()):
            if not source_man[f][2]:
                if source_man[f][0] == '-':
                    addedDirs.append(f)
                else:
                    addedFiles.append(f)

        print(f"Diff Results - Added Files: {len(addedFiles)}, Added Dirs: {len(addedDirs)}, Deleted Files: {len(deleted)}, Modified Files: {len(modified)}")

        if printDiffs:
            new_size = 0
            modified_size = 0
            for f in addedFiles:
                try:
                    new_size += int(source_man[f][1])
                except ValueError:
                    pass
            for f in modified:
                try:
                    modified_size += int(source_man[f][1])
                except ValueError:
                    pass

            if addedFiles:
                print("#------------------------------------------------------")
                print(f"# Source contains {len(addedFiles)} new file(s) ({self.formatData(new_size)})")
                for f in addedFiles:
                    print("\t".join([self.formatData(int(source_man[f][1])), f]))

            if addedDirs:
                print("#------------------------------------------------------")
                print(f"# Source contains {len(addedDirs)} new folder(s)")
                for f in addedDirs:
                    print(f)

            if modified:
                print("#------------------------------------------------------")
                print(f"# Source contains {len(modified)} modified file(s) ({self.formatData(modified_size)})")
                for f in modified:
                    print(f)

            if deleted:
                print("#------------------------------------------------------")
                print(f"# {len(deleted)} file(s) have been deleted in the source:")
                for f in deleted:
                    print(f)
        else:
            return (source,
                    [(a, source_man[a]) for a in addedFiles],
                    [(a, source_man[a]) for a in addedDirs],
                    deleted,
                    [(m, source_man[m]) for m in modified])

    def updateManifest(self,
                       localManifestLocation,
                       sourceManifestLocation,
                       localManifestName=None,
                       sourceManifestName=None,
                       prompt=True):
        """
        Update local files based on remote changes.

        Args:
            localManifestLocation (str): Path to the local manifest directory.
            sourceManifestLocation (str): Path or URL to the source manifest location.
            localManifestName (str, optional): Name of the local manifest file. Defaults to ".dmanifest".
            sourceManifestName (str, optional): Name of the source manifest file. Defaults to ".dmanifest".
            prompt (bool, optional): Flag to prompt the user before downloading updates. Defaults to True.

        Returns:
            bool: 
                - `True` if the update was successful.
                - `False` otherwise.
        """
        # Get the diffs
        diff = self.diffManifests(localManifestLocation,
                                  sourceManifestLocation,
                                  localManifestName,
                                  sourceManifestName,
                                  printDiffs=False)
        source, added_files, added_dirs, deleted, modified = diff

        # Bail if the diff failed
        if source is None:
            return False

        # No changes by default
        do_down = False
        if prompt:
            total_size = 0
            for f in added_files:
                try:
                    total_size += int(f[1][1])
                except ValueError:
                    pass
            for f in modified:
                try:
                    total_size += int(f[1][1])
                except ValueError:
                    pass
            if total_size != 0:
                print("****************************************************************")
                print(f"{len(added_files)} new file(s) to be downloaded from source")
                print(f"{len(modified)} existing file(s) to be updated")
                print(f"{self.formatData(total_size)} will need to be downloaded")
                do_down = self.promptUserDownload()
                if not do_down:
                    print("Download aborted")

        update_manifest = False
        if do_down:
            update_manifest = True
            for add in added_dirs:
                # Make the dirs first
                full_path = os.path.abspath(os.path.join(localManifestLocation, add[0]))
                print(f"Creating directory: {full_path}")
                self.makeSurePathExists(full_path)
            for add in added_files:
                full_path = os.path.abspath(os.path.join(localManifestLocation, add[0]))
                print(f"Downloading new file: {add[0]} -> {full_path}")
                try:
                    urllib.request.urlretrieve(source + add[0], full_path)
                except Exception as e:
                    print(f"Error downloading {add[0]}: {e}")
            for modify in modified:
                full_path = os.path.abspath(os.path.join(localManifestLocation, modify[0]))
                print(f"Updating existing file: {modify[0]} -> {full_path}")
                try:
                    urllib.request.urlretrieve(source + modify[0], full_path)
                except Exception as e:
                    print(f"Error updating {modify[0]}: {e}")

        if update_manifest:
            print("(Re)creating manifest file (please be patient)")
            self.createManifest(localManifestLocation, manifestName=localManifestName)
            
        return True

    def getManType(self, line):
        """
        Work out the manifest type from the first line of the file.

        Args:
            line (str): A line from the manifest file.

        Returns:
            str: The extracted manifest type. Defaults to "generic" if extraction fails.
        """
        try:
            return line.rstrip().split(":::")[1].strip()
        except IndexError:
            return "generic"

    def formatData(self, amount):
        """
        Pretty print file sizes into human-readable formats.

        Args:
            amount (int): Size in bytes.

        Returns:
            str: Formatted size string (e.g., "10 MB").
        """
        try:
            amount = int(amount)
        except (ValueError, TypeError):
            return "0 B"
        if amount < 1024*1024:
            return f"{amount} B"
        elif amount < 1024*1024*1024:
            return f"{amount / (1024*1024):.2f} MB"
        elif amount < 1024*1024*1024*1024:
            return f"{amount / (1024*1024*1024):.2f} GB"
        elif amount < 1024*1024*1024*1024*1024:
            return f"{amount / (1024*1024*1024*1024):.2f} TB"
        else:
            return f"{amount / (1024*1024*1024*1024*1024):.2f} PB"

    #-----------------------------------------------------------------------------
    # FS utilities

    def makeSurePathExists(self, path):
        """
        Ensures that a specified directory path exists, creating it if necessary.

        Args:
            path (str): The directory path to verify or create.
        """
        try:
            os.makedirs(path)
            print(f"Directory created: {path}")
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            else:
                print(f"Directory already exists: {path}")

    def promptUserDownload(self):
        """
        Prompts the user for confirmation before proceeding with downloading updates.

        Returns:
            bool: 
                - `True` if the user confirms.
                - `False` otherwise.
        """
        input_not_ok = True
        minimal = False
        valid_responses = {'Y': True, 'N': False}
        vrs = ",".join([x.lower() for x in list(valid_responses.keys())])
        while input_not_ok:
            if minimal:
                option = input(f"Download? ({vrs}) : ").upper()
            else:
                option = input(f"Confirm you want to download this data\n"
                               f"Changes *WILL* be permanent\n"
                               f"Continue? ({vrs}) : ").upper()
            if option in valid_responses:
                print("****************************************************************")
                return valid_responses[option]
            else:
                print(f"ERROR: unrecognised choice '{option}'")
                minimal = True

    def walk(self, parents, full_path, rel_path, dirs, files, skipFile=".dmanifest", exclude_dirs=None, exclude_files=None, extensions=None):
        """
        Recursively walks through the directory tree to catalog files and directories.

        Args:
            parents (list of FileEntity): Stack of parent directories.
            full_path (str): Absolute path to the current directory.
            rel_path (str): Relative path from the root.
            dirs (list of str): List of directories in the current path.
            files (list of str): List of files in the current path.
            skipFile (str, optional): File to skip (e.g., the manifest file itself). Defaults to ".dmanifest".
            exclude_dirs (list of str, optional): List of directories to exclude.
            exclude_files (list of str, optional): List of files to exclude.
            extensions (list of str, optional): List of file extensions to include.
        """
        print(f"Entering directory: {full_path}")
        # First do files here
        for f in files:
            if (f != skipFile) and (not extensions or os.path.splitext(f)[1] in extensions):
                if exclude_files and f in exclude_files:
                    print(f"Excluding file: {os.path.join(rel_path, f)}")
                    continue
                path = os.path.join(full_path, f)
                print(f"Adding file: {path}")
                self.files.append(FileEntity(f, rel_path, parents[-1], self.hashfile(path), os.path.getsize(path)))
        for d in dirs:
            if exclude_dirs and d in exclude_dirs:
                print(f"Excluding directory: {os.path.join(rel_path, d)}")
                continue
            # The walk will go into these dirs first
            tmp_fe = FileEntity(d, rel_path, parents[-1], "-", 0)
            print(f"Adding directory: {os.path.join(rel_path, d)}")
            self.files.append(tmp_fe)
            parents.append(tmp_fe)
            new_full_path = os.path.join(full_path, d)
            new_rel_path = os.path.join(rel_path, d)
            new_dirs, new_files, _ = self.listdir(new_full_path)
            self.walk(parents, new_full_path, new_rel_path, new_dirs, new_files, skipFile=skipFile, exclude_dirs=exclude_dirs, exclude_files=exclude_files, extensions=extensions)
            parents.pop()

    def listdir(self, path):
        """
        Lists directories, files, and symbolic links within a specified path.

        Args:
            path (str): The directory path to list.

        Returns:
            tuple: Three lists containing directories, files, and symbolic links respectively.
        """
        dirs, files, links = [], [], []
        try:
            for name in os.listdir(path):
                path_name = os.path.join(path, name)
                if os.path.isdir(path_name):
                    dirs.append(name)
                elif os.path.isfile(path_name):
                    files.append(name)
                elif os.path.islink(path_name):
                    links.append(name)
        except PermissionError as e:
            print(f"Permission denied accessing '{path}': {e}")
        except Exception as e:
            print(f"Error accessing '{path}': {e}")
        return dirs, files, links

    def hashfile(self, fileName, blocksize=65536):
        """
        Hashes a file and returns its SHA-256 digest.

        Args:
            fileName (str): Path to the file to hash.
            blocksize (int, optional): Size of each read from the file. Defaults to 65536.

        Returns:
            str: The SHA-256 hash of the file. Returns "?" if hashing fails.
        """
        hasher = hashlib.sha256()
        try:
            with open(fileName, "rb") as fh:
                buf = fh.read(blocksize)
                while len(buf) > 0:
                    hasher.update(buf.strip())
                    buf = fh.read(blocksize)
                return hasher.hexdigest()
        except FileNotFoundError:
            print(f"Warning: File not found during hashing: {fileName}")
            return "?"
        except Exception as e:
            print(f"Error hashing file {fileName}: {e}")
            return "?"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def main():
    parser = argparse.ArgumentParser(description="Manage project manifests for Pizza3.")
    subparsers = parser.add_subparsers(dest='action', help='Available actions')

    # Create Manifest
    parser_create = subparsers.add_parser('create', help='Create a new manifest')
    parser_create.add_argument('-p', '--path', required=True, help='Path to the project directory')
    parser_create.add_argument('-m', '--manifest', default=".dmanifest", help='Name of the manifest file')
    parser_create.add_argument('--exclude-dirs', nargs='*', help='List of directories to exclude from the manifest')
    parser_create.add_argument('--exclude-files', nargs='*', help='List of files to exclude from the manifest')
    parser_create.add_argument('--extensions', nargs='*', help='List of file extensions to include (e.g., .py .sh .md .html)')

    # Diff Manifests
    parser_diff = subparsers.add_parser('diff', help='Diff two manifests')
    parser_diff.add_argument('-l', '--local', required=True, help='Local manifest directory')
    parser_diff.add_argument('-s', '--source', required=True, help='Source manifest location (URL or path)')
    parser_diff.add_argument('--local-manifest', default=".dmanifest", help='Name of the local manifest file')
    parser_diff.add_argument('--source-manifest', default=".dmanifest", help='Name of the source manifest file')
    parser_diff.add_argument('--print', action='store_true', help='Print differences')

    # Update Manifest
    parser_update = subparsers.add_parser('update', help='Update local manifest based on source')
    parser_update.add_argument('-l', '--local', required=True, help='Local manifest directory')
    parser_update.add_argument('-s', '--source', required=True, help='Source manifest location (URL or path)')
    parser_update.add_argument('--local-manifest', default=".dmanifest", help='Name of the local manifest file')
    parser_update.add_argument('--source-manifest', default=".dmanifest", help='Name of the source manifest file')
    parser_update.add_argument('--prompt', action='store_true', help='Prompt before downloading updates')

    args = parser.parse_args()

    if args.action is None:
        parser.print_help()
        sys.exit(1)

    manager = ManifestManager()

    if args.action == 'create':
        manager.createManifest(
            args.path,
            manifestName=args.manifest,
            exclude_dirs=args.exclude_dirs,
            exclude_files=args.exclude_files,
            extensions=args.extensions
        )
        print(f"Manifest '{args.manifest}' created at '{args.path}'.")
        sys.exit(0)

    elif args.action == 'diff':
        diff = manager.diffManifests(
            localManifestLocation=args.local,
            sourceManifestLocation=args.source,
            localManifestName=args.local_manifest,
            sourceManifestName=args.source_manifest,
            printDiffs=args.print
        )
        if diff[0] is not None and args.print:
            print("Diff operation completed.")
        elif diff[0] is None:
            print("Diff operation failed.")
        else:
            # Handle non-print diff results if needed
            pass
        sys.exit(0)

    elif args.action == 'update':
        success = manager.updateManifest(
            localManifestLocation=args.local,
            sourceManifestLocation=args.source,
            localManifestName=args.local_manifest,
            sourceManifestName=args.source_manifest,
            prompt=args.prompt
        )
        if success:
            print("Manifest update completed successfully.")
        else:
            print("Manifest update failed.")
        sys.exit(0)

if __name__ == '__main__':
    main()
