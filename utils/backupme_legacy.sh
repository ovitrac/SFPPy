#!/usr/bin/env bash

# backup all codes
# INRAE\Olivier Vitrac - olivier.vitrac@agroparistech.fr

#-o  -iname 'dump*frames'  \

mainfolder="./.."
backfolder=$mainfolder/history
mkdir -p $backfolder

find $mainfolder -type f      \
    -iname '*.m'    \
-o  -iname '*.asv'  \
-o  -iname '*.m~'   \
-o  -iname '*.pynb' \
-o  -iname '*.py'   \
-o  -iname '*.sh'   \
-o  -iname '*.txt'  \
| zip -rTgp -9 "$backfolder/${PWD##*/}_backup_`whoami`@`hostname`_`date +"%Y_%m_%d__%H-%M"`.zip" -@