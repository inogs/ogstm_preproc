#! /bin/bash

MASKFILE=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask.nc
INPUTDIR=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/inputs
OUTPUTDIR=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/BC/OUTPUT

python atlantic_generator.py -i $INPUTDIR -o $OUTPUTDIR -m $MASKFILE 
