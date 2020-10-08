#! /bin/bash
# Generation of a climatology at 1/24


CLIM_1_Km=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_Climatology_1km.nc
MASKFILE=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc
CLIM_1_24=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_Climatology_24.nc

# inside a serial job

python interp_climatology.py -i $CLIM_1_Km -m $MASKFILE -o $CLIM_1_24

# ... should be continued using scripts in bit.sea/Sat/KD490 