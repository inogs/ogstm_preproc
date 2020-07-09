#! /bin/bash

ORIGDIR=/gpfs/scratch/userexternal/gcoidess/REA_24/REA_24_SAT/KD490/DAILY/ORIG
RAW_DATA=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/raw_data
CLIM_DEC=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/CLIM_DEC
MASKFILE=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_meshmask.nc
CLIMFILE=/gpfs/scratch/userexternal/gbolzon0/REA24/SAT/Climatology_CHL.nc

# step 1 - inside a serial job
python mask_and_counter.py -i $ORIGDIR  -v KD490 -o KD_map_occurency.npy

# step 2 - a few seconds
python mask_gen.py -i KD_map_occurency.npy -m $MASKFILE -o /gpfs/scratch/userexternal/gbolzon0/REA24/SAT/KD490/KD490_occurrency.nc


# from here ahead domain decomposition by 5x5 is hardcoded
# Step 3 -- 1.5 h on 25 cores
mpirun -np 25 python setup_rawdata.py -i $ORIGDIR -o $RAW_DATA -m $MASKFILE -v KD490 --valid_max 10.0

# Step 4 -- 3 h on 25 cores
mpirun -np 25 python climatology.py  -i $ORIGDIR -r $RAW_DATA -o $CLIM_DEC -m $MASKFILE -v KD490

# Step 5 -- a few seconds
python output_writer.py -i $CLIM_DEC -m $MASKFILE -o $CLIMFILE -v KD490
