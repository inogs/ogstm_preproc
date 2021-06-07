#!/bin/bash

#SBATCH -N 6 -n 48
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=115gb
#SBATCH --account=IscrB_MED21BIO
#SBATCH --partition=bdw_usr_prod


cd $SLURM_SUBMIT_DIR

source /marconi_scratch/userexternal/gbolzon0/2017_RA16/CODE/ogstm/compilers/machine_modules/marconi.intel


cp /marconi_scratch/userexternal/gbolzon0/2017_RA16/CODE/ogstm/bin/ogstm.xx .
ulimit -s unlimited
date

export OMP_NUM_THREADS=4
python writeStartEnd.py
mpirun -np 48 ./ogstm.xx  -genvall

date


