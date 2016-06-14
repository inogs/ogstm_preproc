#! /bin/bash

HERE=$PWD
ORIG_DIR=/pico/scratch/userexternal/plazzari/REANALISI-FREE-SURFACE/FORCINGS/UNZIPPED
cd $ORIG_DIR
DATE=1999
#DATE=20150[56]
ls -1 $DATE*U.nc > $HERE/nomefile_U
ls -1 $DATE*V.nc > $HERE/nomefile_V
ls -1 $DATE*W.nc > $HERE/nomefile_W
ls -1 $DATE*T.nc > $HERE/nomefile_T

