# __init__.py
# Patankar library for Python 3.x
#   It is the porting of a collection of tools known as *pantankar.m tools
#   from the migration toolbox
#
#   INRAE\olivier.vitrac@agroparistech.fr

# For the end-user:
#    1) add patankar/ to your PYTHONPATH
#    export PYTHONPATH="${PYTHONPATH}:/home/me/python/patankar"
#    2) in your Python code
#    import patankar *

# list of public classes: data, dump, raster, struct

# $ last revision - 2022-02-14 $

from patankar.layer import layer
from patankar.private.struct import struct