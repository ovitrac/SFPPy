#!/usr/bin/env bash

# backup all codes
# INRAE\Olivier Vitrac - olivier.vitrac@agroparistech.fr

find . -type f      \
    -iname '*.m'    \
-o  -iname '*.asv'  \
-o  -iname '*.m~'   \
-o  -iname '*.pynb' \
-o  -iname '*.py'   \
-o  -iname '*.sh'   \
-o  -iname '*.txt'  \
| zip -rTgp -9 "${PWD##*/}_backup_`hostname`_`date +"%Y_%m_%d__%H-%M"`.zip" -@
