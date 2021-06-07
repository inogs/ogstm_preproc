#! /bin/bash
INITDIR=/gpfs/scratch/userexternal/plazzari/eas_v12/INIT
for I in `ls $INITDIR/*nc`; do ln -s  $I ; done
